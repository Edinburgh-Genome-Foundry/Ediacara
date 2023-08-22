# Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of Ediacara.
#
# Ediacara is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Ediacara is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

import os
import pandas as pd

from Bio import SeqIO
import ediacara as edi

data_dir = os.path.join("tests", "data")

samplesheet_csv = os.path.join(data_dir, "entries.csv")
params_projectname = "EGF test"


def test_write_sequencinggroup_report(tmpdir):
    pdf_file = os.path.join(str(tmpdir), "report.pdf")
    entries = pd.read_csv(samplesheet_csv, header=None)

    entries.columns = [
        "projectname",
        "entry",
        "barcode",
        "sample",
        "fasta",
        "vcf",
        "paf",
        "tsv",
        "consensus_fasta",
    ]
    # have them in order in the pdf
    entries.sort_values(by=["barcode", "sample"], inplace=True)

    comparatorgroups = []
    for index, row in entries.iterrows():
        entry = row["entry"]
        sample = row["sample"]
        vcf = os.path.join(data_dir, row["vcf"])

        reference_gb = os.path.join(data_dir, row["sample"] + ".gb")
        record = SeqIO.read(reference_gb, "genbank")
        record.id = sample
        references = {record.id: record}

        tsv_file = os.path.join(data_dir, row["tsv"])
        paf_path = os.path.join(data_dir, row["paf"])

        tsv = edi.ComparatorGroup.load_tsv(tsv_file)
        paf = edi.ComparatorGroup.load_paf(paf_path)

        assembly_paths = {sample: os.path.join(data_dir, row["consensus_fasta"])}
        vcf_paths = {sample: vcf}

        comparator_group = edi.ComparatorGroup(
            references=references,
            alignments={"paf": paf, "tsv": tsv},
            barcode=row["barcode"],
            assembly_paths=assembly_paths,
            vcf_paths=vcf_paths,
        )

        list_of_constructs = [sample]
        for element in list_of_constructs:
            comparator_group.add_comparator(element)

        comparatorgroups += [comparator_group]

    # Create PDF report
    sequencinggroup = edi.SequencingGroup(comparatorgroups, name=params_projectname)
    sequencinggroup.perform_all_comparisons_in_sequencinggroup()
    edi.write_sequencinggroup_report(target=pdf_file, sequencinggroup=sequencinggroup)


def test_write_assembly_analysis_report(tmpdir):

    samplesheet_csv = os.path.join(data_dir, "review", "entries_review.csv")
    assembly_plan_path = os.path.join(data_dir, "review", "assembly_plan.csv")
    params_projectname = "Test review"
    pdf_file = os.path.join(str(tmpdir), "review_report.pdf")

    entries = pd.read_csv(samplesheet_csv, header=None)
    entries.columns = [
        "project",
        "entry",
        "barcode",
        "sample",
        "result",
        "gb",
        "fa",
        "consensus",
        "paf",
    ]
    # have them in order in the pdf
    entries.sort_values(by=["barcode", "sample"], inplace=True)

    consensus_list = []
    for index, row in entries.iterrows():
        assembly = edi.Assembly(
            assembly_path=os.path.join(data_dir, row["consensus"]),
            reference_path=os.path.join(data_dir, row["gb"]),
            alignment_path=os.path.join(data_dir, "review", row["paf"]),
            assembly_plan=assembly_plan_path,
        )
        consensus_list += [assembly]

    assemblybatch = edi.AssemblyBatch(
        assemblies=consensus_list, name=params_projectname
    )
    assemblybatch.perform_all_interpretations_in_group()

    edi.write_assembly_analysis_report(pdf_file, assemblybatch)
