# Phylogen Virus Fetcher

This repository contains a small Streamlit application that queries NCBI for virus sequences by taxon name. The app uses NCBI's Entrez API to download sequences and metadata as well as a reference sequence with its annotated features.

## Requirements

```
pip install -r requirements.txt
```

## Usage

Run the Streamlit application:

```
streamlit run streamlit_app.py
```

Provide your NCBI API key and a taxon name (e.g. `Potyvirus`). The app will
produce several downloadable files which are also bundled into a single zip
archive:

- `<Taxon_Name>.fasta` – all sequences with isolate, host, country,
  collection date and release date in the headers. Metadata is gathered from
  both `esummary` and GenBank records to ensure these fields are populated when
  available.
- `<Taxon_Name>_sequences_features.txt` – features for every sequence returned
  for the taxon.
- `<Taxon_Name>_refseq.fasta` – the RefSeq sequence for the taxon (if available).
- `<Taxon_Name>_features.txt` – list of features (UTRs, CDS, mat_peptides, etc.) from the RefSeq record.
- `<Taxon_Name>_outputs.zip` – archive containing all of the above files for easy download.
