import pandas as pd

from Bio import SeqFeature
from Bio import SeqIO

import dna_features_viewer

from .Comparator import ComparatorGroup


class AssemblyTranslator(dna_features_viewer.BiopythonTranslator):
    """Custom translator to display only the aligned parts."""

    def compute_feature_color(self, feature):
        ediacara_qualifier = "ediacara"
        if feature.qualifiers[ediacara_qualifier] == "wrong_part":
            return "#f5685e"
        elif feature.qualifiers[ediacara_qualifier] == "correct_part":
            return "#79d300"
        elif feature.qualifiers[ediacara_qualifier] == "unknown_part":
            return "#FFFFFF"
        elif feature.qualifiers[ediacara_qualifier] == "reference":
            return "#d3d3d3"
        else:
            return "#7245dc"  # default dna_features_viewer colour


class AssemblyBatch:
    """Batch of Assembly class instances.
    
    
    **Parameters**
    
    **assemblies**
    > List of `Assembly` instances.

    **name**
    > Name of the assembly analysis project (`str`).
    """

    def __init__(self, assemblies, name="Unnamed"):
        self.assemblies = assemblies
        self.name = name

    def perform_all_interpretations_in_group(self):
        for assembly in self.assemblies:
            assembly.interpret_alignment()


class Assembly:
    """Compare an assembly sequence to parts and simulated reference.


    **Parameters**

    **assembly_path**
    > Path to the *de novo* assembled FASTA sequence file (`str`).

    **reference**
    > Path to reference Genbank file (`str`).
    
    **alignment**
    > Path (`str`) to minimap2 alignment PAF file, created with the `-cx asm5` options.
    In this, the part and reference sequences are aligned against the de novo sequence.

    **assembly_plan**
    > Path of assembly plan CSV file (DNA Cauldron output format, with header line).
    May contain additional entries, the correct line is chosen by reference sequence ID.
    """

    def __init__(
        self, assembly_path, reference_path, alignment_path, assembly_plan=None
    ):
        self.assembly = SeqIO.read(handle=assembly_path, format="fasta")
        self.reference = SeqIO.read(handle=reference_path, format="genbank")
        self.paf = ComparatorGroup.load_paf(alignment_path)
        self.paf.columns = self.paf.columns[:-1].to_list() + ["CIGAR"]  # -1 is last

        if assembly_plan is None:
            self.assembly_plan = None
        else:
            plan_df = pd.read_csv(assembly_plan, skiprows=1, header=None)  # skip header
            self.assembly_plan = plan_df[plan_df[0] == self.reference.id]
            if len(self.assembly_plan) == 0:
                raise ValueError("Error! Assembly plan doesn't contain the reference!")
            if len(self.assembly_plan) > 1:
                raise ValueError(
                    "Error! More than one assembly plan entry matches the reference!"
                )

    def interpret_alignment(self):
        for index, row in self.paf.iterrows():
            # May be useful to exclude reference:
            # if row["query_name"] == self.reference.id:
            #     continue
            part_type = self.evaluate_part_name(row["query_name"])

            location = SeqFeature.FeatureLocation(
                row["target_start"], row["target_end"]
            )
            feature = SeqFeature.SeqFeature(
                location=location,
                type="misc_feature",
                id=row["query_name"],
                qualifiers={"label": row["query_name"], "ediacara": part_type},
            )

            self.assembly.features.append(feature)
        self.assembly_figure = self.plot_assembly()

    def subset_paf(self):
        selected_columns = [
            "query_name",
            "query_length",
            "query_start",
            "query_end",
            "strand",
            "target_start",
            "target_end",
            "mapping_matches",
            "mapping_size",
            "mapping_quality",
        ]
        paf_subset = self.paf[selected_columns]
        new_columnnames = [
            "Name",
            "Length",
            "Start",
            "End",
            "Strand",
            "T Start",
            "T End",
            "Matches",
            "Size",
            "Quality",
        ]
        paf_subset.columns = new_columnnames

        return paf_subset

    def evaluate_part_name(self, name):
        if self.assembly_plan is None:
            return "unknown_part"
        elif name == self.reference.id:
            return "reference"
        elif name in self.assembly_plan.iloc[0].to_list():  # has only one line
            return "correct_part"
        else:
            return "wrong_part"

    def plot_assembly(self):
        graphic_record = AssemblyTranslator().translate_record(self.assembly)
        ax, _ = graphic_record.plot(figure_width=8, strand_in_label_threshold=7)
        return ax
