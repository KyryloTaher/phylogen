import streamlit as st
from Bio import Entrez, SeqIO
import io
import os
import zipfile

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

def main():
    st.title("NCBI Virus Fetcher")
    api_key = st.text_input("NCBI API Key")
    taxon = st.text_input("Taxon Name", "Potyvirus")
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
        base = taxon.replace(" ", "_")
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

if __name__ == "__main__":
    main()
