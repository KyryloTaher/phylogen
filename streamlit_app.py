import streamlit as st
from Bio import Entrez, SeqIO
import io
import os

def fetch_ids(taxon, api_key):
    search_term = f"\"{taxon}\"[Organism]"
    handle = Entrez.esearch(db="nuccore", term=search_term, retmax=10000, api_key=api_key)
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])

def fetch_metadata(ids, api_key):
    metadata = {}
    for chunk in [ids[i:i+200] for i in range(0, len(ids), 200)]:
        handle = Entrez.esummary(db="nuccore", id=",".join(chunk), retmode="json", api_key=api_key)
        data = handle.read()
        handle.close()
        import json
        js = json.loads(data)
        for uid in js.get("result", {}).get("uids", []):
            info = js["result"].get(uid, {})
            meta = {}
            labels = info.get("subtype", "").split("|") if info.get("subtype") else []
            values = info.get("subname", "").split("|") if info.get("subname") else []
            for i, lab in enumerate(labels):
                val = values[i] if i < len(values) else ""
                meta[lab] = val
            subtype = info.get("subtype", "").split("|")
            subname = info.get("subname", "").split("|")
            for k, v in zip(subtype, subname):
                meta[k] = v
            meta["release_date"] = info.get("createdate", "")
            metadata[uid] = meta
    return metadata

def fetch_fasta(ids, api_key):
    sequences = []
    for chunk in [ids[i:i+500] for i in range(0, len(ids), 500)]:
        handle = Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="fasta", retmode="text", api_key=api_key)
        data = handle.read()
        handle.close()
        sequences.append(data)
    return "".join(sequences)

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
            qual = "; ".join(f"{k}={','.join(v)}" for k, v in feat.qualifiers.items())
            features.append(f"{feat.type}\t{feat.location}\t{qual}")
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
        base = taxon.replace(" ", "_")
        fasta_file = f"{base}.fasta"
        with open(fasta_file, "w") as f:
            for record in SeqIO.parse(io.StringIO(fasta_data), "fasta"):
                meta = metadata.get(record.id, {})
                header = (
                    f"{record.id}"
                    f" | isolate={meta.get('isolate','')}"
                    f" | host={meta.get('host','')}"
                    f" | country={meta.get('country','')}"
                    f" | collection_date={meta.get('collection_date','')}"
                    f" | release_date={meta.get('release_date','')}"
                )

                header = f"{record.description} | host={meta.get('host','')} | country={meta.get('country','')} | collection_date={meta.get('collection_date','')} | release_date={meta.get('release_date','')}"
                f.write(f">{header}\n{record.seq}\n")
        st.download_button("Download Sequences", open(fasta_file, "rb"), file_name=fasta_file)
        ref_id, ref_fasta, features = fetch_refseq(taxon, api_key)
        if ref_id:
            ref_file = f"{base}_refseq.fasta"
            with open(ref_file, "w") as f:
                f.write(ref_fasta)
            feat_file = f"{base}_features.txt"
            with open(feat_file, "w") as f:
                f.write(features)
            st.download_button("Download RefSeq", open(ref_file, "rb"), file_name=ref_file)
            st.download_button("Download Features", open(feat_file, "rb"), file_name=feat_file)
        else:
            st.write("No refseq found for this taxon")

if __name__ == "__main__":
    main()
