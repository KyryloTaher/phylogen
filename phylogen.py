from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import io
import os
import tempfile
import subprocess
import shutil
from datetime import datetime
from pathlib import Path

def fetch_ids(taxon, api_key):
    """Return accession.version identifiers for the given taxon."""

    search_term = f"\"{taxon}\"[Organism]"
    handle = Entrez.esearch(
        db="nuccore", term=search_term, retmax=10000, api_key=api_key
    )
    record = Entrez.read(handle)
    handle.close()
    uid_list = record.get("IdList", [])
    if not uid_list:
        return []

    accessions = []
    for chunk in [uid_list[i:i + 200] for i in range(0, len(uid_list), 200)]:
        handle = Entrez.esummary(db="nuccore", id=",".join(chunk), api_key=api_key)
        summaries = Entrez.read(handle)
        handle.close()
        for summary in summaries:
            acc = summary.get("AccessionVersion") or summary.get("Caption")
            if acc:
                accessions.append(acc)

    return accessions

def fetch_metadata(ids, api_key):
    """Return metadata for each accession in *ids*.

    The function first queries ``esummary`` to retrieve the basic mapping from
    UID to accession and any available metadata.  Some records lack host,
    country or collection date in ``esummary`` so we additionally fetch the
    GenBank record and parse the ``source`` feature to fill in any missing
    fields.
    """

    import json

    metadata = {}

    # Fetch summary metadata for quick access to release dates and mapping to
    # accession.version identifiers.
    for chunk in [ids[i:i + 200] for i in range(0, len(ids), 200)]:
        handle = Entrez.esummary(
            db="nuccore", id=",".join(chunk), retmode="json", api_key=api_key
        )
        js = json.load(handle)
        handle.close()

        for uid in js.get("result", {}).get("uids", []):
            info = js["result"].get(uid, {})
            meta = {}
            subtype = info.get("subtype", "").split("|")
            subname = info.get("subname", "").split("|")
            for k, v in zip(subtype, subname):
                meta[k] = v
            meta["release_date"] = info.get("createdate", "")
            acc = info.get("accessionversion", uid)
            metadata[acc] = meta

    # Retrieve GenBank records to capture missing metadata such as host or
    # collection_date.
    for chunk in [list(metadata.keys())[i:i + 100] for i in range(0, len(metadata), 100)]:
        handle = Entrez.efetch(
            db="nuccore", id=",".join(chunk), rettype="gb", retmode="text", api_key=api_key
        )
        for record in SeqIO.parse(handle, "genbank"):
            meta = metadata.setdefault(record.id, {})
            # Fill release date from record annotations if missing.
            meta.setdefault("release_date", record.annotations.get("date", ""))
            for feat in record.features:
                if feat.type == "source":
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
        handle.close()

    return metadata

def fetch_fasta(ids, api_key):
    sequences = []
    for chunk in [ids[i:i+500] for i in range(0, len(ids), 500)]:
        handle = Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="fasta", retmode="text", api_key=api_key)
        data = handle.read()
        handle.close()
        sequences.append(data)
    return "".join(sequences)

def fetch_all_features(ids, api_key):
    """Return features for all sequences in *ids* as a single string."""
    lines = []
    for chunk in [ids[i:i + 100] for i in range(0, len(ids), 100)]:
        handle = Entrez.efetch(
            db="nuccore", id=",".join(chunk), rettype="gb", retmode="text", api_key=api_key
        )
        for record in SeqIO.parse(handle, "genbank"):
            lines.append(f">{record.id}")
            for feat in record.features:
                lines.append(f"{feat.type}: {feat.location}")
        handle.close()
    return "\n".join(lines)

def fetch_refseq(taxon, api_key):
    """Return RefSeq record information including proteins."""

    term = f"\"{taxon}\"[Organism] AND srcdb_refseq[PROP]"
    handle = Entrez.esearch(db="nuccore", term=term, retmax=1, api_key=api_key)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return None, None, None, []

    refseq_id = ids[0]

    handle = Entrez.efetch(
        db="nuccore", id=refseq_id, rettype="gb", retmode="text", api_key=api_key
    )
    gb = handle.read()
    handle.close()

    fasta = None
    features = []
    proteins = []
    for record in SeqIO.parse(io.StringIO(gb), "genbank"):
        fasta = f">{record.id}\n{record.seq}\n"
        for feat in record.features:
            features.append(f"{feat.type}: {feat.location}")
            if feat.type == "mat_peptide":
                seq = feat.extract(record.seq)
                codon_start = int(feat.qualifiers.get("codon_start", ["1"])[0])
                if codon_start > 1:
                    seq = seq[codon_start - 1:]
                if len(seq) % 3:
                    pad = 3 - len(seq) % 3
                    seq = seq + Seq("N" * pad)
                prot_seq = seq.translate(to_stop=False)
                prot_id = feat.qualifiers.get("protein_id", [""])[0]
                gene = feat.qualifiers.get("gene", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                label = gene or product or prot_id or f"protein_{len(proteins)+1}"
                label = label.replace(" ", "_")
                proteins.append(SeqRecord(prot_seq, id=label))

    return refseq_id, fasta, "\n".join(features), proteins


def _extract_mat_peptide_features(record):
    """Return mat_peptide sequences, codon_start values and coordinates."""

    pep_list = []
    for feat in record.features:
        if feat.type == "mat_peptide":
            seq = feat.extract(record.seq)
            codon_start = int(feat.qualifiers.get("codon_start", ["1"])[0])
            start = int(feat.location.nofuzzy_start) + 1
            end = int(feat.location.nofuzzy_end)
            coords = f"{start}-{end}"
            pep_list.append((seq, codon_start, coords))
    if not pep_list:
        coords = f"1-{len(record.seq)}"
        pep_list.append((record.seq, 1, coords))
    return pep_list


def _run_mafft(records):
    """Align a list of SeqRecord objects with MAFFT L-INS-i."""
    if not records:
        return []
    mafft_path = shutil.which("mafft")
    if mafft_path is None:
        raise FileNotFoundError(
            "MAFFT executable not found. Please install MAFFT and ensure it is in your PATH."
        )
    with tempfile.TemporaryDirectory() as tmpdir:
        in_f = os.path.join(tmpdir, "in.fasta")
        SeqIO.write(records, in_f, "fasta")
        cmd = [mafft_path, "--localpair", "--maxiterate", "1000", in_f]
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
        return list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))


def _align_translate_back_with_ref(cds_records, ref_protein=None):
    """Translate mat_peptide regions, include an optional reference protein and align."""

    aa_records = []
    codons = {}
    ref_id = None

    if ref_protein is not None:
        ref_id = ref_protein.id
        aa_records.append(ref_protein)

    for rec in cds_records:
        seq = rec.seq
        codon_start = int(rec.annotations.get("codon_start", 1))
        if codon_start > 1:
            seq = seq[codon_start - 1:]
        if len(seq) % 3:
            pad = 3 - len(seq) % 3
            seq = seq + Seq("N" * pad)
        codon_list = [str(seq[i:i + 3]) for i in range(0, len(seq), 3)]
        codons[rec.id] = codon_list
        aa_records.append(SeqRecord(seq.translate(to_stop=False), id=rec.id))

    aligned_aa = _run_mafft(aa_records)

    aligned_nt = {}
    ref_aligned = None
    for rec in aligned_aa:
        if ref_id is not None and rec.id == ref_id:
            ref_aligned = rec.seq
            continue
        codon_list = codons[rec.id]
        idx = 0
        nt_frag = []
        for aa in str(rec.seq):
            if aa == "-":
                nt_frag.append("---")
            else:
                nt_frag.append(codon_list[idx])
                idx += 1
        aligned_nt[rec.id] = SeqRecord(Seq("".join(nt_frag)), id=rec.id)

    return aligned_nt, ref_aligned


def align_mat_peptides(fasta_file, ids, api_key, ref_proteins=None):
    """Align every mat_peptide to each RefSeq mat_peptide individually."""

    records = [rec for rec in SeqIO.parse(fasta_file, "fasta") if rec.id in ids]

    pep_parts = {}
    for chunk in [ids[i:i + 50] for i in range(0, len(ids), 50)]:
        handle = Entrez.efetch(
            db="nuccore", id=",".join(chunk), rettype="gb", retmode="text", api_key=api_key
        )
        for record in SeqIO.parse(handle, "genbank"):
            pep_parts[record.id] = _extract_mat_peptide_features(record)
        handle.close()

    outputs = {}

    if not ref_proteins:
        return outputs

    all_peptides = []
    for sid, parts in pep_parts.items():
        for seq, codon_start, coords in parts:
            label = f"{sid}/_{coords}"
            all_peptides.append(
                SeqRecord(seq, id=label, annotations={"codon_start": codon_start})
            )

    for ref_rec in ref_proteins:
        aligned, ref_aln = _align_translate_back_with_ref(all_peptides, ref_rec)
        final_records = []

        if ref_aln is not None:
            final_records.append(SeqRecord(ref_aln, id=f"RefSeq|{ref_rec.id}"))

        for pep in all_peptides:
            if pep.id in aligned:
                final_records.append(SeqRecord(aligned[pep.id].seq, id=pep.id))

        output = io.StringIO()
        SeqIO.write(final_records, output, "fasta")
        outputs[ref_rec.id] = output.getvalue()

    return outputs

def main():
    print("NCBI Virus Fetcher")
    api_key = input("NCBI API Key: ").strip()
    taxon = input("Taxon Name [Potyvirus]: ").strip() or "Potyvirus"
    base = taxon.replace(" ", "_")

    # Determine output directory within the user's home/projects/phylogen folder
    base_dir = Path.home() / "projects" / "phylogen"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = base_dir / f"outputs_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Outputs will be stored in {output_dir}")
    if not api_key or not taxon:
        print("Please provide both API key and taxon name")
        return

    Entrez.email = "example@example.com"

    print("Fetching data from NCBI...")
    ids = fetch_ids(taxon, api_key)
    if not ids:
        print("No sequences found")
        return
    print(f"Found {len(ids)} sequences")

    metadata = fetch_metadata(ids, api_key)
    fasta_data = fetch_fasta(ids, api_key)
    features_data = fetch_all_features(ids, api_key)

    fasta_file = output_dir / f"{base}.fasta"
    with open(fasta_file, "w") as f:
        for record in SeqIO.parse(io.StringIO(fasta_data), "fasta"):
            meta = metadata.get(record.id, {})
            header = (
                f"{record.id} | isolate={meta.get('isolate','')} | "
                f"host={meta.get('host','')} | "
                f"country={meta.get('country','')} | "
                f"collection_date={meta.get('collection_date','')} | "
                f"release_date={meta.get('release_date','')}"
            )
            f.write(f">{header}\n{record.seq}\n")
    print(f"Sequences written to {fasta_file}")

    seq_feat_file = output_dir / f"{base}_sequences_features.txt"
    with open(seq_feat_file, "w") as f:
        f.write(features_data)
    print(f"Sequence features written to {seq_feat_file}")

    ref_id, ref_fasta, features, ref_proteins = fetch_refseq(taxon, api_key)
    if ref_id:
        ref_file = output_dir / f"{base}_refseq.fasta"
        with open(ref_file, "w") as f:
            f.write(ref_fasta)
        feat_file = output_dir / f"{base}_refseq_features.txt"
        with open(feat_file, "w") as f:
            f.write(features)
        print(f"RefSeq written to {ref_file}")
        print(f"RefSeq features written to {feat_file}")
    else:
        print("No refseq found for this taxon")


    choice = input("Align mat_peptide sequences only? [y/N]: ").strip().lower()
    if choice == "y" and ref_proteins:
        ids_no_ref = [i for i in ids if i != ref_id]
        alignments = align_mat_peptides(fasta_file, ids_no_ref, api_key, ref_proteins)
        for rname, data in alignments.items():
            safe_name = rname.replace(" ", "_")
            align_file = output_dir / f"{base}_{safe_name}_alignment.fasta"
            with open(align_file, "w") as af:
                af.write(data)
            print(f"Alignment for {rname} written to {align_file}")

if __name__ == "__main__":
    main()

