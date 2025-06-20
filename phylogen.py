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
    """Return RefSeq record information including proteins and mat_peptide positions."""

    term = f"\"{taxon}\"[Organism] AND srcdb_refseq[PROP]"
    handle = Entrez.esearch(db="nuccore", term=term, retmax=1, api_key=api_key)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return None, None, None, []

    refseq_uid = ids[0]

    handle = Entrez.efetch(
        db="nuccore", id=refseq_uid, rettype="gb", retmode="text", api_key=api_key
    )
    gb = handle.read()
    handle.close()

    fasta = None
    features = []
    proteins = []
    pep_positions = []
    refseq_acc = None
    for record in SeqIO.parse(io.StringIO(gb), "genbank"):
        refseq_acc = record.id
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
                pep_positions.append(
                    {
                        "label": label,
                        "start": int(feat.location.start),
                        "end": int(feat.location.end),
                        "codon_start": codon_start,
                    }
                )

    return refseq_acc, fasta, "\n".join(features), proteins, pep_positions

def _extract_mat_peptide_features(record):
    """Return mat_peptide sequences and their codon_start values."""

    pep_list = []
    for feat in record.features:
        if feat.type == "mat_peptide":
            seq = feat.extract(record.seq)
            codon_start = int(feat.qualifiers.get("codon_start", ["1"])[0])
            pep_list.append((seq, codon_start))
    if not pep_list:
        pep_list.append((record.seq, 1))
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


def _choose_frame_no_stop(seq):
    """Return codon_start (1-based) for frame without stop codons."""
    s = str(seq).replace("-", "")
    for frame in range(3):
        sub = s[frame:]
        sub = sub[: len(sub) // 3 * 3]
        if not sub:
            return frame + 1
        aa = Seq(sub).translate(to_stop=False)
        if "*" not in str(aa):
            return frame + 1
    return 1


def align_mat_peptides_two_step(
    fasta_file, ids, api_key, ref_id, ref_pep_positions, ref_proteins=None
):
    """Two-step alignment of mat_peptides based on RefSeq coordinates.

    Returns a tuple `(concatenated, per_peptide)` where `concatenated` is a FASTA-formatted string with all peptides joined per accession and `per_peptide` is a list of `(label, fasta_string)` for each individual RefSeq mat_peptide alignment.
    """
    rec_dict = {
        rec.id: rec
        for rec in SeqIO.parse(fasta_file, "fasta")
        if rec.id in ids or rec.id == ref_id
    }
    if ref_id not in rec_dict:
        # fetch reference sequence if not present in fasta_file
        handle = Entrez.efetch(
            db="nuccore",
            id=ref_id,
            rettype="fasta",
            retmode="text",
            api_key=api_key,
        )
        rec = SeqIO.read(handle, "fasta")
        handle.close()
        rec_dict[ref_id] = rec
    order = [i for i in ids if i != ref_id]

    # first alignment of full sequences including refseq
    full_aln = _run_mafft([rec_dict[ref_id]] + [rec_dict[i] for i in order])
    aln_dict = {r.id: r.seq for r in full_aln}
    ref_aln = aln_dict[ref_id]

    # map reference positions to alignment columns
    pos_to_col = {}
    pos = 0
    for idx, nt in enumerate(str(ref_aln)):
        if nt != "-":
            pos_to_col[pos] = idx
            pos += 1

    concatenated = {sid: [] for sid in order}
    per_peptide_outputs = []

    for idx, feat in enumerate(ref_pep_positions):
        start = feat["start"]
        end = feat["end"]
        start_col = pos_to_col.get(start)
        end_col = pos_to_col.get(end - 1)
        if start_col is None or end_col is None:
            continue
        end_col += 1

        pep_records = []
        for sid in order:
            sub = aln_dict[sid][start_col:end_col]
            raw = Seq(str(sub).replace("-", ""))
            frame = _choose_frame_no_stop(raw)
            pep_records.append(
                SeqRecord(raw, id=sid, annotations={"codon_start": frame})
            )

        ref_rec = None
        if ref_proteins and idx < len(ref_proteins):
            ref_rec = ref_proteins[idx]

        aligned_nt, ref_aln_seq = _align_translate_back_with_ref(pep_records, ref_rec)

        part_len = 0
        if ref_aln_seq is not None:
            part_len = len(ref_aln_seq) * 3
        elif aligned_nt:
            part_len = len(next(iter(aligned_nt.values())).seq)

        per_records = []
        if ref_aln_seq is not None:
            per_records.append(SeqRecord(ref_aln_seq, id=f"RefSeq|{ref_rec.id}"))
        for sid in order:
            if sid in aligned_nt:
                seq_str = str(aligned_nt[sid].seq)
                concatenated[sid].append(seq_str)
                per_records.append(SeqRecord(aligned_nt[sid].seq, id=sid))
            else:
                gap = "-" * part_len
                concatenated[sid].append(gap)
                per_records.append(SeqRecord(Seq(gap), id=sid))
        out = io.StringIO()
        SeqIO.write(per_records, out, "fasta")
        per_peptide_outputs.append((feat["label"], out.getvalue()))

    final_records = [SeqRecord(Seq("".join(v)), id=k) for k, v in concatenated.items()]

    output = io.StringIO()
    SeqIO.write(final_records, output, "fasta")
    return output.getvalue(), per_peptide_outputs


def align_mat_peptides(fasta_file, ids, api_key, ref_proteins=None):
    """Retrieve mat_peptide regions, align each separately and include RefSeq proteins."""

    records = [rec for rec in SeqIO.parse(fasta_file, "fasta") if rec.id in ids]
    order = [rec.id for rec in records]

    pep_parts = {}
    for chunk in [ids[i:i + 50] for i in range(0, len(ids), 50)]:
        handle = Entrez.efetch(
            db="nuccore", id=",".join(chunk), rettype="gb", retmode="text", api_key=api_key
        )
        for record in SeqIO.parse(handle, "genbank"):
            pep_parts[record.id] = _extract_mat_peptide_features(record)
        handle.close()

    max_parts = max(len(v) for v in pep_parts.values()) if pep_parts else 0
    n_parts = max(max_parts, len(ref_proteins) if ref_proteins else 0)

    final_records = []

    for idx in range(n_parts):
        pep_recs = []
        for sid in order:
            if idx < len(pep_parts.get(sid, [])):
                seq, codon_start = pep_parts[sid][idx]
                pep_recs.append(
                    SeqRecord(seq, id=sid, annotations={"codon_start": codon_start})
                )

        ref_rec = None
        if ref_proteins and idx < len(ref_proteins):
            ref_rec = ref_proteins[idx]

        if not pep_recs and ref_rec is None:
            continue

        aligned, ref_aln = _align_translate_back_with_ref(pep_recs, ref_rec)

        part_len = None
        if ref_aln is not None:
            part_len = len(ref_aln) * 3
        elif aligned:
            part_len = len(next(iter(aligned.values())).seq)

        if ref_aln is not None:
            final_records.append(
                SeqRecord(ref_aln, id=f"RefSeq|{ref_rec.id}")
            )

        for sid in order:
            label = f"{sid}|MP{idx + 1}"
            if ref_rec is not None:
                label = f"{sid}|{ref_rec.id}"
            if sid in aligned:
                final_records.append(SeqRecord(aligned[sid].seq, id=label))
            else:
                if part_len is not None:
                    final_records.append(SeqRecord(Seq("-" * part_len), id=label))

    output = io.StringIO()
    SeqIO.write(final_records, output, "fasta")
    return output.getvalue()

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

    ref_id, ref_fasta, features, ref_proteins, pep_positions = fetch_refseq(taxon, api_key)
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
    if choice == "y":
        ids_no_ref = [i for i in ids if i != ref_id]
        aligned, per_pep = align_mat_peptides_two_step(
            fasta_file,
            ids_no_ref,
            api_key,
            ref_id,
            pep_positions,
            ref_proteins,
        )
        align_file = output_dir / f"{base}_mat_peptide_alignment.fasta"
        with open(align_file, "w") as af:
            af.write(aligned)
        print(f"mat_peptide alignment written to {align_file}")

        for label, aln in per_pep:
            fname = label.replace("|", "_")
            part_file = output_dir / f"{base}_{fname}_alignment.fasta"
            with open(part_file, "w") as pf:
                pf.write(aln)
            print(f"mat_peptide {label} alignment written to {part_file}")

if __name__ == "__main__":
    main()

