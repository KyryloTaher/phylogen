# Phylogen Virus Fetcher

This repository contains a small command line tool that queries NCBI for virus sequences by taxon name. It uses NCBI's Entrez API to download sequences and metadata as well as a reference sequence with its annotated features.

## Requirements

```
pip install -r requirements.txt
```

[MAFFT](https://mafft.cbrc.jp/alignment/software/) executable to be
installed and available on your `PATH`.

## Usage

Run the script:

```
python phylogen.py
```

The program will prompt for your NCBI API key and a taxon name (e.g. `Potyvirus`).
It produces several files:

- `<Taxon_Name>.fasta` – all sequences with isolate, host, country,
  collection date and release date in the headers. Metadata is gathered from
  both `esummary` and GenBank records to ensure these fields are populated when
  available.
- `<Taxon_Name>_sequences_features.txt` – features for every sequence returned
  for the taxon.
- `<Taxon_Name>_refseq.fasta` – the RefSeq sequence for the taxon (if available).
- `<Taxon_Name>_refseq_features.txt` – list of features (UTRs, CDS, mat_peptides, etc.) from the RefSeq record.
- RefSeq peptide labels have spaces replaced with underscores so MAFFT does not
  truncate them during alignment.

Optionally, the script can align all mature peptide regions (`mat_peptide` features).
Each peptide is aligned individually with the corresponding RefSeq nucleotide segment
placed at the top of the alignment. The resulting file `<Taxon_Name>_mat_peptide_alignment.fasta`
contains the nucleotide alignments for every peptide with gaps added when a sequence
lacks a particular peptide.
Additionally, separate alignment files `<Taxon_Name>_<Peptide_Label>_alignment.fasta`
are produced for each RefSeq mat-peptide.
