# Ediacara

**Work in progress**

The **Edi**nburgh Genome Foundry's **a**lignment **c**omp**ara**tor.


## Install

```
# pip install ediacara
pip install --upgrade git+https://github.com/Edinburgh-Genome-Foundry/Ediacara.git@main
```
Using Ediacara requires the plasmid reference sequences (Genbank), the alignments in [PAF format](https://lh3.github.io/minimap2/minimap2.html#10), and a TSV file of coverage (depth) counts.

[Minimap2](https://lh3.github.io/minimap2/) can create SAM or PAF alignments from fastq and reference fasta files.
[`paftools.js`](https://github.com/lh3/minimap2/blob/master/misc/README.md) can convert from SAM to PAF.
The TSV file is the output of `samtools depth -aa`.

[Badread](https://github.com/rrwick/Badread) is a long read simulator that generates fastq files from reference sequences.


## Usage

An example for creating alignments for a reference named `ABC_1`:
```
# badread simulate --reference ABC_1.fa --quantity 50x --length 1500,500 > ABC_1.fastq
minimap2 -a ABC_1.fa ABC_1.fastq > ABC_1.sam
samtools sort -O sam -T sample.sort -o ABC_1_sorted.sam ABC_1.sam
samtools depth -aa ABC_1_sorted.sam > ABC_1.tsv
paftools.js sam2paf ABC_1.sam > ABC_1.paf
```

Run Ediacara and plot coverage:
```python
import ediacara as edi
import dnacauldron

tsv_file = "/path/to/ABC_1.tsv"
tsv = edi.ComparatorGroup.load_tsv(tsv_file)
paf_path = "/path/to/ABC_1.paf"
paf = edi.ComparatorGroup.load_paf(paf_path)
reference_paths = ["/path/to/ABC_1.gb"]
records = dnacauldron.biotools.load_records_from_files(files=reference_paths, use_file_names_as_ids=True)

references = {record.id: record for record in records}
alignments = {"paf": paf, "tsv": tsv}
comparator_group = edi.ComparatorGroup(references, alignments)
comparator = comparator_group.create_comparator("ABC_1")

edi.write_pdf_report("report.pdf", comparator)
```


## Versioning

Ediacara uses the [semantic versioning](https://semver.org) scheme.


## Copyright

Copyright 2021 Edinburgh Genome Foundry

Ediacara was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp).
