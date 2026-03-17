import io
import json
import time
from datetime import datetime
from pathlib import Path

import requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class GalaxyMAFFTClient:
    """Run MAFFT alignments on a Galaxy server via API."""

    def __init__(self, url: str, api_key: str, timeout: int = 900):
        self.base_url = url.rstrip("/")
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"x-api-key": api_key})

    def _get(self, path: str, **kwargs):
        response = self.session.get(f"{self.base_url}{path}", timeout=60, **kwargs)
        response.raise_for_status()
        return response

    def _post(self, path: str, payload: dict):
        response = self.session.post(
            f"{self.base_url}{path}", json=payload, timeout=60
        )
        response.raise_for_status()
        return response

    def _find_mafft_tool_id(self) -> str:
        response = self._get("/api/tools", params={"q": "mafft", "qv": "mafft"})
        tools = response.json()
        if not tools:
            raise RuntimeError("No MAFFT tool found on this Galaxy server.")

        for tool in tools:
            name = (tool.get("name") or "").lower()
            tool_id = tool.get("id", "")
            if "mafft" in name or "mafft" in tool_id.lower():
                return tool_id

        raise RuntimeError("Could not identify a MAFFT tool ID from Galaxy API response.")

    def _first_data_input_name(self, tool_id: str, history_id: str) -> str:
        response = self._get(f"/api/tools/{tool_id}/build", params={"history_id": history_id})
        inputs = response.json().get("inputs", [])
        for inp in inputs:
            if inp.get("type") == "data":
                return inp["name"]
        raise RuntimeError("MAFFT tool does not expose a data input field.")

    def _wait_dataset_ok(self, dataset_id: str):
        deadline = time.time() + self.timeout
        while time.time() < deadline:
            ds = self._get(f"/api/datasets/{dataset_id}").json()
            state = ds.get("state")
            if state == "ok":
                return
            if state in {"error", "failed", "discarded"}:
                raise RuntimeError(f"Galaxy dataset {dataset_id} failed with state '{state}'.")
            time.sleep(4)
        raise TimeoutError(f"Timed out waiting for dataset {dataset_id} to complete.")

    def align_fasta(self, fasta_text: str) -> str:
        history = self._post("/api/histories", {"name": f"phylogen_mafft_{int(time.time())}"}).json()
        history_id = history["id"]

        upload_payload = {
            "history_id": history_id,
            "tool_id": "upload1",
            "inputs": {
                "files_0|type": "upload_dataset",
                "files_0|url_paste": fasta_text,
                "file_type": "fasta",
            },
        }
        upload_result = self._post("/api/tools", upload_payload).json()
        input_dataset_id = upload_result["outputs"][0]["id"]
        self._wait_dataset_ok(input_dataset_id)

        mafft_tool_id = self._find_mafft_tool_id()
        input_name = self._first_data_input_name(mafft_tool_id, history_id)

        run_payload = {
            "history_id": history_id,
            "tool_id": mafft_tool_id,
            "inputs": {input_name: {"src": "hda", "id": input_dataset_id}},
        }
        run_result = self._post("/api/tools", run_payload).json()
        if not run_result.get("outputs"):
            raise RuntimeError("Galaxy MAFFT run produced no outputs.")

        out_dataset_id = run_result["outputs"][0]["id"]
        self._wait_dataset_ok(out_dataset_id)

        data = self._get(f"/api/datasets/{out_dataset_id}/display").text
        self._post(f"/api/histories/{history_id}", {"deleted": True, "purged": True})
        return data


def fetch_ids(taxon, api_key):
    handle = Entrez.esearch(db="nuccore", term=f'"{taxon}"[Organism]', retmax=10000, api_key=api_key)
    record = Entrez.read(handle)
    handle.close()
    uids = record.get("IdList", [])
    if not uids:
        return []

    accessions = []
    for i in range(0, len(uids), 200):
        chunk = uids[i : i + 200]
        h = Entrez.esummary(db="nuccore", id=",".join(chunk), api_key=api_key)
        summaries = Entrez.read(h)
        h.close()
        for summary in summaries:
            acc = summary.get("AccessionVersion") or summary.get("Caption")
            if acc:
                accessions.append(acc)
    return accessions


def fetch_metadata(ids, api_key):
    metadata = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i : i + 200]
        handle = Entrez.esummary(
            db="nuccore", id=",".join(chunk), retmode="json", api_key=api_key
        )
        js = json.load(handle)
        handle.close()

        for uid in js.get("result", {}).get("uids", []):
            info = js["result"].get(uid, {})
            meta = {"release_date": info.get("createdate", "")}
            subtype = info.get("subtype", "").split("|")
            subname = info.get("subname", "").split("|")
            for k, v in zip(subtype, subname):
                meta[k] = v
            acc = info.get("accessionversion", uid)
            metadata[acc] = meta

    for i in range(0, len(ids), 100):
        chunk = ids[i : i + 100]
        h = Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="gb", retmode="text", api_key=api_key)
        for record in SeqIO.parse(h, "genbank"):
            meta = metadata.setdefault(record.id, {})
            meta.setdefault("release_date", record.annotations.get("date", ""))
            for feat in record.features:
                if feat.type != "source":
                    continue
                q = feat.qualifiers
                if "isolate" in q and not meta.get("isolate"):
                    meta["isolate"] = q["isolate"][0]
                if "host" in q and not meta.get("host"):
                    meta["host"] = q["host"][0]
                if "country" in q and not meta.get("country"):
                    meta["country"] = q["country"][0]
                if "geo_loc_name" in q and not meta.get("country"):
                    meta["country"] = q["geo_loc_name"][0]
                if "collection_date" in q and not meta.get("collection_date"):
                    meta["collection_date"] = q["collection_date"][0]
                break
        h.close()

    return metadata


def fetch_fasta(ids, api_key):
    parts = []
    for i in range(0, len(ids), 500):
        chunk = ids[i : i + 500]
        h = Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="fasta", retmode="text", api_key=api_key)
        parts.append(h.read())
        h.close()
    return "".join(parts)


def fetch_labeled_peptides(ids, api_key):
    pep_map = {}
    labels = set()
    for i in range(0, len(ids), 50):
        chunk = ids[i : i + 50]
        h = Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="gb", retmode="text", api_key=api_key)
        for record in SeqIO.parse(h, "genbank"):
            rec_labels = {}
            for feat in record.features:
                if feat.type != "mat_peptide":
                    continue
                seq = feat.extract(record.seq)
                codon_start = int(feat.qualifiers.get("codon_start", ["1"])[0])
                if codon_start > 1:
                    seq = seq[codon_start - 1 :]
                if len(seq) % 3:
                    seq = seq + Seq("N" * (3 - len(seq) % 3))
                prot_id = feat.qualifiers.get("protein_id", [""])[0]
                gene = feat.qualifiers.get("gene", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                label = (gene or product or prot_id or "protein").replace(" ", "_")
                rec_labels[label] = SeqRecord(seq, id=record.id)
                labels.add(label)
            if not rec_labels:
                rec_labels["genome"] = SeqRecord(record.seq, id=record.id)
                labels.add("genome")
            pep_map[record.id] = rec_labels
        h.close()
    return pep_map, sorted(labels)


def prompt_canonical_labels(labels):
    print("Detected mat_peptide labels:")
    for label in labels:
        print(f" - {label}")
    canonical = {}
    ordered = []
    for label in labels:
        alias = input(f"Canonical name for '{label}' (empty=keep): ").strip() or label
        canonical[label] = alias
        if alias not in ordered:
            ordered.append(alias)
    return canonical, ordered


def align_peptides_with_galaxy(ids, api_key, galaxy_client):
    pep_map, labels = fetch_labeled_peptides(ids, api_key)
    canonical, ordered_groups = prompt_canonical_labels(labels)

    concatenated = {sid: [] for sid in ids}
    per_outputs = []

    for group in ordered_groups:
        records = []
        for sid in ids:
            match = None
            for label, rec in pep_map.get(sid, {}).items():
                if canonical[label] == group:
                    match = rec
                    break
            if match:
                records.append(match)

        if not records:
            continue

        in_buf = io.StringIO()
        SeqIO.write(records, in_buf, "fasta")
        aligned_text = galaxy_client.align_fasta(in_buf.getvalue())
        aligned = {r.id: r for r in SeqIO.parse(io.StringIO(aligned_text), "fasta")}
        part_len = len(next(iter(aligned.values())).seq) if aligned else 0

        out_records = []
        for sid in ids:
            if sid in aligned:
                seq = aligned[sid].seq
            else:
                seq = Seq("-" * part_len)
            concatenated[sid].append(str(seq))
            out_records.append(SeqRecord(seq, id=sid))

        out = io.StringIO()
        SeqIO.write(out_records, out, "fasta")
        per_outputs.append((group, out.getvalue()))

    final = io.StringIO()
    final_records = [
        SeqRecord(Seq("".join(parts)), id=sid) for sid, parts in concatenated.items() if parts
    ]
    SeqIO.write(final_records, final, "fasta")
    return final.getvalue(), per_outputs


def main():
    print("NCBI Virus Fetcher + Galaxy MAFFT")
    Entrez.email = input("Email for NCBI Entrez: ").strip() or "anonymous@example.com"
    ncbi_key = input("NCBI API Key: ").strip()
    galaxy_url = input("Galaxy URL [https://usegalaxy.org]: ").strip() or "https://usegalaxy.org"
    galaxy_key = input("Galaxy API Key: ").strip()
    taxon = input("Taxon Name [Potyvirus]: ").strip() or "Potyvirus"

    base = taxon.replace(" ", "_")
    out_dir = Path.home() / "projects" / "phylogen" / f"outputs_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out_dir}")

    ids = fetch_ids(taxon, ncbi_key)
    if not ids:
        print("No sequences found.")
        return
    print(f"Found {len(ids)} sequences.")

    metadata = fetch_metadata(ids, ncbi_key)
    fasta_data = fetch_fasta(ids, ncbi_key)

    fasta_file = out_dir / f"{base}.fasta"
    with open(fasta_file, "w") as f:
        for rec in SeqIO.parse(io.StringIO(fasta_data), "fasta"):
            meta = metadata.get(rec.id, {})
            header = (
                f"{rec.id} | isolate={meta.get('isolate', '')} | host={meta.get('host', '')} | "
                f"country={meta.get('country', '')} | collection_date={meta.get('collection_date', '')} | "
                f"release_date={meta.get('release_date', '')}"
            )
            f.write(f">{header}\n{rec.seq}\n")
    print(f"Sequences written to {fasta_file}")

    if input("Align mat_peptides via Galaxy MAFFT? [y/N]: ").strip().lower() != "y":
        return

    galaxy = GalaxyMAFFTClient(galaxy_url, galaxy_key)
    alignment, per_peptide = align_peptides_with_galaxy(ids, ncbi_key, galaxy)

    aln_file = out_dir / f"{base}_mat_peptide_alignment.fasta"
    with open(aln_file, "w") as f:
        f.write(alignment)
    print(f"Combined alignment written to {aln_file}")

    for label, aln in per_peptide:
        name = label.replace("|", "_")
        part_file = out_dir / f"{base}_{name}_alignment.fasta"
        with open(part_file, "w") as f:
            f.write(aln)
        print(f"{label} alignment written to {part_file}")


if __name__ == "__main__":
    main()
