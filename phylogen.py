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

def fetch_ids(taxon, api_key, refseq_only=False):
    """Return accession.version identifiers for the given taxon.

    If *refseq_only* is True, restrict the search to RefSeq records.
    """

    search_term = f"\"{taxon}\"[Organism]"
    if refseq_only:
        search_term += " AND srcdb_refseq[PROP]"
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
    """Return RefSeq record information including mat_peptide nucleotide segments."""

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
                prot_id = feat.qualifiers.get("protein_id", [""])[0]
                gene = feat.qualifiers.get("gene", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                label = gene or product or prot_id or f"protein_{len(proteins)+1}"
                label = label.replace(" ", "_")
                proteins.append(
                    SeqRecord(seq, id=label, annotations={"codon_start": codon_start})
                )
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


def _extract_labeled_mat_peptides(record):
    """Return mapping of label -> SeqRecord for mat_peptides."""

    result = {}
    for feat in record.features:
        if feat.type == "mat_peptide":
            seq = feat.extract(record.seq)
            codon_start = int(feat.qualifiers.get("codon_start", ["1"])[0])
            if codon_start > 1:
                seq = seq[codon_start - 1:]
            if len(seq) % 3:
                pad = 3 - len(seq) % 3
                seq = seq + Seq("N" * pad)
            prot_id = feat.qualifiers.get("protein_id", [""])[0]
            gene = feat.qualifiers.get("gene", [""])[0]
            product = feat.qualifiers.get("product", [""])[0]
            label = gene or product or prot_id or f"protein_{len(result)+1}"
            label = label.replace(" ", "_")
            result[label] = SeqRecord(seq, id=label, annotations={"codon_start": codon_start})
    if not result:
        result["genome"] = SeqRecord(record.seq, id="genome", annotations={"codon_start": 1})
    return result


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


def _align_translate_back_with_ref(cds_records, ref_record=None):
    """Translate mat_peptide regions, include an optional reference segment and align."""

    aa_records = []
    codons = {}
    ref_id = None

    if ref_record is not None:
        ref_id = ref_record.id
        seq = ref_record.seq
        codon_start = int(ref_record.annotations.get("codon_start", 1))
        if codon_start > 1:
            seq = seq[codon_start - 1:]
        if len(seq) % 3:
            pad = 3 - len(seq) % 3
            seq = seq + Seq("N" * pad)
        codons[ref_record.id] = [str(seq[i:i + 3]) for i in range(0, len(seq), 3)]
        aa_records.append(SeqRecord(seq.translate(to_stop=False), id=ref_record.id))

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
        codon_list = codons[rec.id]
        idx = 0
        nt_frag = []
        for aa in str(rec.seq):
            if aa == "-":
                nt_frag.append("---")
            else:
                nt_frag.append(codon_list[idx])
                idx += 1
        seq_nt = Seq("".join(nt_frag))
        if ref_id is not None and rec.id == ref_id:
            ref_aligned = SeqRecord(seq_nt, id=rec.id)
        else:
            aligned_nt[rec.id] = SeqRecord(seq_nt, id=rec.id)


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


def _compute_identity_map(aln_str):
    """Return mapping of sequence ID to identity with the reference.

    The identity is calculated using only the region from the first to last
    non-gap nucleotide in each sequence. Gaps are ignored in the calculation.
    """
    recs = list(SeqIO.parse(io.StringIO(aln_str), "fasta"))
    ref = None
    for r in recs:
        if r.id.startswith("RefSeq|"):
            ref = r
            break
    if ref is None:
        return {}

    ref_seq = str(ref.seq)
    result = {}
    for rec in recs:
        if rec.id == ref.id:
            continue
        seq = str(rec.seq)
        try:
            first = next(i for i, c in enumerate(seq) if c != "-")
            last = len(seq) - 1 - next(i for i, c in enumerate(reversed(seq)) if c != "-")
        except StopIteration:
            result[rec.id] = 0.0
            continue

        sub = seq[first : last + 1]
        ref_sub = ref_seq[first : last + 1]
        matches = 0
        valid = 0
        for a, b in zip(sub, ref_sub):
            if a != "-" and b != "-":
                valid += 1
                if a == b:
                    matches += 1
        result[rec.id] = matches / valid if valid else 0.0
    return result


def _filter_alignment(aln_str, keep_ids):
    """Return FASTA alignment string with only the specified IDs kept."""
    recs = [
        r
        for r in SeqIO.parse(io.StringIO(aln_str), "fasta")
        if r.id.startswith("RefSeq|") or r.id in keep_ids
    ]
    out = io.StringIO()
    SeqIO.write(recs, out, "fasta")
    return out.getvalue()


def _prompt_label_groups(labels):
    """Return mapping of label -> canonical label using user input.

    The user is first shown all distinct labels and then asked to provide a
    canonical name for each one. Pressing enter keeps the original label.
    """

    print("Mat_peptide labels detected:")
    for lab in sorted(labels):
        print(f" - {lab}")

    canonical = {}
    order = []
    for label in sorted(labels):
        ans = input(f"Canonical name for '{label}' (leave empty to keep): ").strip()
        if not ans:
            ans = label
        canonical[label] = ans
        if ans not in order:
            order.append(ans)
    return canonical, order


def align_mat_peptides_by_name(ids, api_key):
    """Align mat_peptides by their feature names across multiple sequences.

    The user will be asked to provide canonical names for each distinct
    mat_peptide label encountered so that peptides with synonymous names can be
    aligned together.
    """
    """Return mapping of label -> canonical label using user input."""

    canonical = {}
    groups = []
    for label in sorted(labels):
        assigned = False
        for canon in groups:
            ans = input(f'Does "{label}" mean the same as "{canon}"? [y/N]: ').strip().lower()
            if ans == "y":
                canonical[label] = canon
                assigned = True
                break
        if not assigned:
            canonical[label] = label
            groups.append(label)
    return canonical, groups


def align_mat_peptides_by_name(ids, api_key):
    """Align mat_peptides by their feature names across multiple sequences."""
    pep_map = {}
    all_labels = set()
    for chunk in [ids[i:i + 50] for i in range(0, len(ids), 50)]:
        handle = Entrez.efetch(
            db="nuccore", id=",".join(chunk), rettype="gb", retmode="text", api_key=api_key
        )
        for record in SeqIO.parse(handle, "genbank"):
            labeled = _extract_labeled_mat_peptides(record)
            pep_map[record.id] = labeled
            all_labels.update(labeled.keys())
        handle.close()

    canonical, order = _prompt_label_groups(all_labels)

    concatenated = {sid: [] for sid in ids}
    completeness = {sid: {} for sid in ids}
    per_outputs = []

    for label in order:
        pep_recs = []
        for sid in ids:
            rec_map = pep_map.get(sid, {})
            found = None
            for lab, rec in rec_map.items():
                if canonical[lab] == label:
                    found = rec
                    break
            if found is not None:
                pep_recs.append(SeqRecord(found.seq, id=sid, annotations=found.annotations))

        aligned, _ = _align_translate_back_with_ref(pep_recs)

        part_len = len(next(iter(aligned.values())).seq) if aligned else 0

        per_records = []
        for sid in ids:
            if sid in aligned:
                seq_str = str(aligned[sid].seq)
            else:
                seq_str = "-" * part_len
            codons = [seq_str[i:i+3] for i in range(0, len(seq_str), 3)]
            comp = False
            if codons and codons[0] != "---" and codons[-1] != "---":
                comp = True
            completeness[sid][label] = comp
            concatenated[sid].append(seq_str)
            per_records.append(SeqRecord(Seq(seq_str), id=sid, description="complete" if comp else "partial"))
        out = io.StringIO()
        SeqIO.write(per_records, out, "fasta")
        per_outputs.append((label, out.getvalue()))

    final_records = []
    for sid, parts in concatenated.items():
        overall = all(completeness[sid].get(lbl, False) for lbl in order)
        final_records.append(
            SeqRecord(Seq("".join(parts)), id=sid, description="complete" if overall else "partial")
        )

    out = io.StringIO()
    SeqIO.write(final_records, out, "fasta")
    return out.getvalue(), per_outputs, completeness


def align_mat_peptides_two_step(
    fasta_file, ids, api_key, ref_id, ref_pep_positions, ref_proteins=None
):
    """Two-step alignment of mat_peptides based on RefSeq coordinates.

    Returns a tuple `(concatenated, per_peptide, completeness)` where
    `concatenated` is a FASTA-formatted string with all peptides joined per
    accession, `per_peptide` is a list of `(label, fasta_string)` for each
    individual RefSeq mat_peptide alignment and `completeness` is a mapping of
    accession -> peptide label -> boolean indicating whether the first and last
    reference amino acids are present.
    """
    """Two-step alignment of mat_peptides based on RefSeq coordinates."""

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
    completeness = {sid: {} for sid in order}
    labels = []


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
            part_len = len(ref_aln_seq.seq)
        elif aligned_nt:
            part_len = len(next(iter(aligned_nt.values())).seq)

        # Determine first and last reference codon positions
        first_idx = last_idx = None
        if ref_aln_seq is not None:
            codons = [str(ref_aln_seq.seq)[i:i+3] for i in range(0, len(ref_aln_seq.seq), 3)]
            for i, c in enumerate(codons):
                if c != '---':
                    first_idx = i
                    break
            for i, c in enumerate(reversed(codons)):
                if c != '---':
                    last_idx = len(codons) - 1 - i
                    break


        per_records = []
        if ref_aln_seq is not None:
            per_records.append(SeqRecord(ref_aln_seq.seq, id=f"RefSeq|{ref_rec.id}"))
        label = feat["label"]
        labels.append(label)
        for sid in order:
            if sid in aligned_nt:
                seq_str = str(aligned_nt[sid].seq)
            else:
                seq_str = "-" * part_len
            codons = [seq_str[i:i+3] for i in range(0, len(seq_str), 3)]
            comp = False
            if first_idx is not None and last_idx is not None and len(codons) > last_idx:
                comp = codons[first_idx] != '---' and codons[last_idx] != '---'
            completeness[sid][label] = comp
            concatenated[sid].append(seq_str)
            per_records.append(SeqRecord(Seq(seq_str), id=sid, description="complete" if comp else "partial"))
        out = io.StringIO()
        SeqIO.write(per_records, out, "fasta")
        per_peptide_outputs.append((label, out.getvalue()))

    final_records = []
    for sid, parts in concatenated.items():
        overall = all(completeness[sid].get(lbl, False) for lbl in labels)
        final_records.append(
            SeqRecord(Seq("".join(parts)), id=sid, description="complete" if overall else "partial")
        )

    output = io.StringIO()
    SeqIO.write(final_records, output, "fasta")
    return output.getvalue(), per_peptide_outputs, completeness

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
            part_len = len(ref_aln.seq)
        elif aligned:
            part_len = len(next(iter(aligned.values())).seq)

        if ref_aln is not None:
            final_records.append(
                SeqRecord(ref_aln.seq, id=f"RefSeq|{ref_rec.id}")
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
    higher = input("Is this genus level or higher? [y/N]: ").strip().lower() == "y"

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
    ids = fetch_ids(taxon, api_key, refseq_only=higher)
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
        filter_opt = input("Filter by completeness? [all/complete]: ").strip().lower() or "all"
        if higher:
            aligned, per_pep, completeness = align_mat_peptides_by_name(
                ids,
                api_key,
            )
        else:
            ids_no_ref = [i for i in ids if i != ref_id]
            aligned, per_pep, completeness = align_mat_peptides_two_step(
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

        if higher:
            id_pass = set(ids)
        else:
            identity_map = _compute_identity_map(aligned)
            id_pass = {sid for sid, val in identity_map.items() if val >= 0.5}
            if id_pass:
                id_file = output_dir / f"{base}_mat_peptide_alignment_id50.fasta"
                with open(id_file, "w") as af:
                    af.write(_filter_alignment(aligned, id_pass))
                print(
                    f"mat_peptide alignment (>=50% identity) written to {id_file}"
                )

        if filter_opt == "complete":
            comp_ids = {
                sid
                for sid, vals in completeness.items()
                if all(vals.get(lbl, False) for lbl in vals)
            }
            keep_ids = comp_ids & id_pass
            filt_records = [
                rec
                for rec in SeqIO.parse(io.StringIO(aligned), "fasta")
                if rec.id in keep_ids or rec.id.startswith("RefSeq|")
            ]
            out = io.StringIO()
            SeqIO.write(filt_records, out, "fasta")
            comp_file = (
                output_dir
                / f"{base}_mat_peptide_alignment_complete.fasta"
            )
            with open(comp_file, "w") as cf:
                cf.write(out.getvalue())
            print(
                f"Complete mat_peptide alignment written to {comp_file}"
            )
        if filter_opt == "complete":
            comp_ids = [sid for sid, vals in completeness.items() if all(vals.get(l, False) for l in vals)]
            filt_records = [
                rec
                for rec in SeqIO.parse(io.StringIO(aligned), "fasta")
                if rec.id in comp_ids or rec.id.startswith("RefSeq|")
            ]
            out = io.StringIO()
            SeqIO.write(filt_records, out, "fasta")
            comp_file = output_dir / f"{base}_mat_peptide_alignment_complete.fasta"
            with open(comp_file, "w") as cf:
                cf.write(out.getvalue())
            print(f"Complete mat_peptide alignment written to {comp_file}")

        for label, aln in per_pep:
            fname = label.replace("|", "_")
            part_file = output_dir / f"{base}_{fname}_alignment.fasta"
            with open(part_file, "w") as pf:
                pf.write(aln)
            print(f"mat_peptide {label} alignment written to {part_file}")

            if not higher and id_pass:
                part_file_i = output_dir / f"{base}_{fname}_alignment_id50.fasta"
                with open(part_file_i, "w") as pfi:
                    pfi.write(_filter_alignment(aln, id_pass))
                print(
                    f"mat_peptide {label} >=50% identity alignment written to {part_file_i}"
                )

            if filter_opt == "complete":
                filt_records = [
                    rec
                    for rec in SeqIO.parse(io.StringIO(aln), "fasta")
                    if rec.id.startswith("RefSeq|") or completeness.get(rec.id, {}).get(label, False)
                ]
                out = io.StringIO()
                SeqIO.write(filt_records, out, "fasta")
                part_file_c = output_dir / f"{base}_{fname}_alignment_complete.fasta"
                with open(part_file_c, "w") as pf2:
                    pf2.write(out.getvalue())
                print(f"mat_peptide {label} complete-only alignment written to {part_file_c}")

if __name__ == "__main__":
    main()

