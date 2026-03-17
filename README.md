# Phylogen Virus Fetcher

This command-line tool queries NCBI for virus sequences by taxon name, enriches sequence headers with metadata, and aligns mature peptide regions (`mat_peptide`) using **MAFFT on a Galaxy server** (default: `https://usegalaxy.org`).

## Requirements

```bash
pip install -r requirements.txt
```

You need:
- NCBI API key
- Galaxy API key (for a Galaxy instance exposing a MAFFT tool)

## Usage

```bash
python phylogen.py
```

The script prompts for:
- Entrez email
- NCBI API key
- Galaxy URL (default `https://usegalaxy.org`)
- Galaxy API key
- Taxon name

## Outputs

The script writes output under:

`~/projects/phylogen/outputs_<timestamp>/`

Main files include:
- `<Taxon_Name>.fasta` — downloaded NCBI sequences with metadata in headers (`isolate`, `host`, `country`, `collection_date`, `release_date`).
- `<Taxon_Name>_mat_peptide_alignment.fasta` — concatenated alignment across peptide groups.
- `<Taxon_Name>_<Peptide_Label>_alignment.fasta` — per-peptide group alignments.

## Notes

- Local MAFFT is no longer used.
- MAFFT runs are submitted via Galaxy API by creating a temporary history, uploading FASTA, running MAFFT, downloading output, and deleting the temporary history.
