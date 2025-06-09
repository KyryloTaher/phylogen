import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import io
import os
import zipfile
import tempfile
import subprocess

def fetch_ids(taxon, api_key):
    search_term = f"\"{taxon}\"[Organism]"
    handle = Entrez.esearch(db="nuccore", term=search_term, retmax=10000, api_key=api_key)
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])

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
    term = f"\"{taxon}\"[Organism] AND srcdb_refseq[PROP]"
    handle = Entrez.esearch(db="nuccore", term=term, retmax=1, api_key=api_key)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return None, None, None
    refseq_id = ids[0]
    handle = Entrez.efetch(db="nuccore", id=refseq_id, rettype="fasta", retmode="text", api_key=api_key)
    fasta = handle.read()
    handle.close()
    handle = Entrez.efetch(db="nuccore", id=refseq_id, rettype="gb", retmode="text", api_key=api_key)
    gb = handle.read()
    handle.close()
    features = []
    for record in SeqIO.parse(io.StringIO(gb), "genbank"):
        for feat in record.features:
            features.append(f"{feat.type}: {feat.location}")
    return refseq_id, fasta, "\n".join(features)


def _parse_partitions(record):
    """Return 5'UTR, CDS and 3'UTR sequences from a GenBank record."""
    utr5 = Seq("")
    utr3 = Seq("")
    cds_seq = None
    cds_loc = None
    for feat in record.features:
        if feat.type == "CDS" and cds_seq is None:
            cds_seq = feat.extract(record.seq)
            cds_loc = feat.location
        elif feat.type in {"5'UTR", "five_prime_UTR"}:
            utr5 = feat.extract(record.seq)
        elif feat.type in {"3'UTR", "three_prime_UTR"}:
            utr3 = feat.extract(record.seq)
    if cds_seq is None:
        cds_seq = record.seq
        cds_loc = record.features[0].location if record.features else None
    if cds_loc is not None:
        start = int(cds_loc.start)
        end = int(cds_loc.end)
        if not utr5:
            utr5 = record.seq[:start]
        if not utr3:
            utr3 = record.seq[end:]
    return utr5, cds_seq, utr3


def _run_mafft(records):
    """Align a list of SeqRecord objects with MAFFT L-INS-i."""
    if not records:
        return []
    with tempfile.TemporaryDirectory() as tmpdir:
        in_f = os.path.join(tmpdir, "in.fasta")
        SeqIO.write(records, in_f, "fasta")
        result = subprocess.run(
            ["mafft", "--localpair", "--maxiterate", "1000", in_f],
            capture_output=True, text=True, check=True
        )
        return list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))


def _align_translate_back(cds_records):
    """Translate CDS, align amino acids and back-translate to nucleotides."""
    aa_records = []
    codons = {}
    for rec in cds_records:
        codon_list = [str(rec.seq[i:i+3]) for i in range(0, len(rec.seq), 3)]
        codons[rec.id] = codon_list
        aa_records.append(SeqRecord(rec.seq.translate(to_stop=False), id=rec.id))
    aligned_aa = _run_mafft(aa_records)
    aligned_nt = {}
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
        aligned_nt[rec.id] = SeqRecord(Seq("".join(nt_frag)), id=rec.id)
    return aligned_nt


def partition_and_align(fasta_file, ids, api_key):
    """Partition sequences into UTRs and CDS and align each part."""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    order = [rec.id for rec in records]
    partitions = {}
    for chunk in [ids[i:i + 50] for i in range(0, len(ids), 50)]:
        handle = Entrez.efetch(
            db="nuccore", id=",".join(chunk), rettype="gb", retmode="text", api_key=api_key
        )
        for record in SeqIO.parse(handle, "genbank"):
            utr5, cds, utr3 = _parse_partitions(record)
            partitions[record.id] = {"utr5": utr5, "cds": cds, "utr3": utr3}
        handle.close()
    utr5_recs = [SeqRecord(partitions[i]["utr5"], id=i) for i in order]
    cds_recs = [SeqRecord(partitions[i]["cds"], id=i) for i in order]
    utr3_recs = [SeqRecord(partitions[i]["utr3"], id=i) for i in order]

    aligned_utr5 = {r.id: r for r in _run_mafft(utr5_recs)}
    aligned_cds = _align_translate_back(cds_recs)
    aligned_utr3 = {r.id: r for r in _run_mafft(utr3_recs)}

    final_records = []
    for i in order:
        seq = (
            aligned_utr5.get(i, SeqRecord(Seq(""))).seq
            + aligned_cds.get(i, SeqRecord(Seq(""))).seq
            + aligned_utr3.get(i, SeqRecord(Seq(""))).seq
        )
        final_records.append(SeqRecord(seq, id=i))

    output = io.StringIO()
    SeqIO.write(final_records, output, "fasta")
    return output.getvalue()

def main():
    st.title("NCBI Virus Fetcher")
    api_key = st.text_input("NCBI API Key")
    taxon = st.text_input("Taxon Name", "Potyvirus")
    base = taxon.replace(" ", "_")
    if st.button("Fetch Data"):
        if not api_key or not taxon:
            st.error("Please provide both API key and taxon name")
            return
        Entrez.email = "example@example.com"
        ids = fetch_ids(taxon, api_key)
        if not ids:
            st.error("No sequences found")
            return
        st.write(f"Found {len(ids)} sequences")
        metadata = fetch_metadata(ids, api_key)
        fasta_data = fetch_fasta(ids, api_key)
        features_data = fetch_all_features(ids, api_key)
        fasta_file = f"{base}.fasta"
        with open(fasta_file, "w") as f:
            for record in SeqIO.parse(io.StringIO(fasta_data), "fasta"):
                meta = metadata.get(record.id, {})
                # Use only the accession ID in the header to avoid long descriptions
                header = (
                    f"{record.id} | isolate={meta.get('isolate','')} | "
                    f"host={meta.get('host','')} | "
                    f"country={meta.get('country','')} | "
                    f"collection_date={meta.get('collection_date','')} | "
                    f"release_date={meta.get('release_date','')}"
                )
                f.write(f">{header}\n{record.seq}\n")
        st.download_button("Download Sequences", open(fasta_file, "rb"), file_name=fasta_file)
        st.session_state['fasta_file'] = fasta_file
        st.session_state['ids'] = ids
        st.session_state['base'] = base
        seq_feat_file = f"{base}_sequences_features.txt"
        with open(seq_feat_file, "w") as f:
            f.write(features_data)
        output_files = [fasta_file, seq_feat_file]
        ref_id, ref_fasta, features = fetch_refseq(taxon, api_key)
        if ref_id:
            ref_file = f"{base}_refseq.fasta"
            with open(ref_file, "w") as f:
                f.write(ref_fasta)
            feat_file = f"{base}_features.txt"
            with open(feat_file, "w") as f:
                f.write(features)
            output_files.extend([ref_file, feat_file])
            st.download_button("Download RefSeq", open(ref_file, "rb"), file_name=ref_file)
            st.download_button("Download Features", open(feat_file, "rb"), file_name=feat_file)
        else:
            st.write("No refseq found for this taxon")
        zip_file = f"{base}_outputs.zip"
        with zipfile.ZipFile(zip_file, "w") as zf:
            for fpath in output_files:
                zf.write(fpath, arcname=os.path.basename(fpath))
        st.download_button("Download All Outputs", open(zip_file, "rb"), file_name=zip_file, mime="application/zip")

    if st.button("Partition and Align"):
        fasta_path = st.session_state.get('fasta_file')
        id_list = st.session_state.get('ids')
        base = st.session_state.get('base', base)
        if not fasta_path or not id_list:
            st.error("Please fetch data first")
        else:
            aligned = partition_and_align(fasta_path, id_list, api_key)
            align_file = f"{base}_partitioned_alignment.fasta"
            with open(align_file, "w") as af:
                af.write(aligned)
            st.download_button(
                "Download Partitioned Alignment",
                open(align_file, "rb"),
                file_name=align_file,
            )

if __name__ == "__main__":
    main()
