<p align="center">
<img alt="Ediacara logo" title="Ediacara" src="images/logo.png" width="120">
</p>

# Ediacara

![version](https://img.shields.io/badge/current_version-0.2.0-blue)
[![build](https://github.com/Edinburgh-Genome-Foundry/Ediacara/actions/workflows/build.yml/badge.svg)](https://github.com/Edinburgh-Genome-Foundry/Ediacara/actions/workflows/build.yml)

Ediacara helps interpreting sequencing data of assembled DNA constructs (plasmids). It's used in the [Sequeduct pipeline](https://github.com/Edinburgh-Genome-Foundry/Sequeduct/).

## Install

```bash
# pip install ediacara
pip install --upgrade git+https://github.com/Edinburgh-Genome-Foundry/Ediacara.git@main
```

The plotting feature requires NCBI BLAST+. On Ubuntu, install it with:

```bash
apt-get install ncbi-blast+
```

## Usage

Using Ediacara requires the plasmid reference sequences (Genbank), the alignments in [PAF format](https://lh3.github.io/minimap2/minimap2.html#10), and a TSV file of coverage (depth) counts. Optionally, the user can also specify consensus (or *de novo* assembly) files that will be compared with the references.

[Minimap2](https://lh3.github.io/minimap2/) can create SAM or PAF alignments from fastq and reference fasta files.
[`paftools.js`](https://github.com/lh3/minimap2/blob/master/misc/README.md) can convert from SAM to PAF.
The TSV file is the output of `samtools depth -aa`.

[Canu](https://canu.readthedocs.io) can create an assembly file from filtered FASTQ reads that align best to the reference.

If we don't have sequencing data, [Badread](https://github.com/rrwick/Badread) can be used to simulate fastq files from reference sequences.

An example for creating alignments in the terminal, for a reference named `ABC_1`:

```bash
# badread simulate --reference ABC_1.fa --quantity 100x --length 5000,500 > ABC_1.fastq
minimap2 -a ABC_1.fa ABC_1.fastq > ABC_1.sam
samtools sort -O sam -T sample.sort -o ABC_1_sorted.sam ABC_1.sam
samtools depth -aa ABC_1_sorted.sam > ABC_1.tsv
paftools.js sam2paf ABC_1.sam > ABC_1.paf
```

Create a *de novo* assembly with Canu:

```bash
./canu-2.1.1/bin/canu -p egf -d ABC_1 genomeSize=8k -nanopore filtered_fastq/ABC_1_filtered.fastq
```

Create a PDF report in Python with Ediacara:

```python
import ediacara as edi
import dnacauldron  # for loading sequence files

tsv_file = "/path/to/ABC_1.tsv"
tsv = edi.ComparatorGroup.load_tsv(tsv_file)
paf_path = "/path/to/ABC_1.paf"
paf = edi.ComparatorGroup.load_paf(paf_path)
reference_paths = ["/path/to/ABC_1.gb"]
records = dnacauldron.biotools.load_records_from_files(files=reference_paths, use_file_names_as_ids=True)

references = {record.id: record for record in records}
alignments = {"paf": paf, "tsv": tsv}
comparator_group = edi.ComparatorGroup(references, alignments)

assembly_paths = {"ABC_1": "/path/to/canu_assembly/ABC_1/egf.contigs.fasta"}
comparator_group.perform_all_comparisons(assembly_paths=assembly_paths)
# Write a PDF report on all constructs:
edi.write_comparatorgroup_report("report_ABC.pdf", comparator_group)
```

<p align="center">
<img alt="Ediacara" title="Ediacara" src="images/analysis_screenshot.png" width="360">
</p>

This page of the report shows a missing part ('feature_14'). The presence of unaligned read segments suggest that there is another ~1.6 kbp sequence in its place (this was the part's carrier plasmid in this case). The correct section also contains a point mutation as seen in the red row of the variant call table. Suspected false positives (systematic sequencing errors) are shown in grey.

If we want to investigate which of the parts are present in the construct, we can align the part sequences against the consensus or de novo sequence, and create PDF report:

```python
assembly = edi.Assembly(
    assembly_path="/path/to/canu_assembly/ABC_1/egf.contigs.fasta",
    reference_path="/path/to/ABC_1.gb",
    alignment_path="/path/to/ABC_1_part_alignment.paf",
    assembly_plan=assembly_plan,
)

assemblybatch = edi.AssemblyBatch(assemblies=[assembly], name="EGF review")
assemblybatch.perform_all_interpretations_in_group()

edi.write_assembly_analysis_report("review_report_pdf", assemblybatch)
```

<p align="center">
<img alt="Ediacara" title="Ediacara" src="images/review_screenshot.png" width="360">
</p>

The plot shows alignment regions as annotations: green annotations (genetic parts) are expected in the construct; red annotations show other parts that align to the construct. Grey shows the alignment of the reference expected sequence for the plasmid. In this case, we have no errors as the green segments cover the whole sequence, and the reference aligns fully as well.

## Versioning

Ediacara uses the [semantic versioning](https://semver.org) scheme.

## License = GPLv3+

Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh

Ediacara was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp), and is released under the GPLv3 license.
