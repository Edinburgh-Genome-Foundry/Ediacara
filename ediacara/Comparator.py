import pandas


class ComparatorGroup:
    """Analyse alignments to a set of references.


    **Parameters**

    **references**
    > Dictionary of name: Biopython record.

    **alignments**
    > Dictionary of pandas dataframes: {"paf": df_paf, "tsv": df_tsv}. The PAF dataframe
    was created from the output of `paftools.js sam2paf` and the TSV dataframe was
    created from the output of `samtools depth -aa`.
    """
    def __init__(self, references, alignments):
        self.references = references
        self.paf = alignments["paf"]
        self.tsv = alignments["tsv"]
        self.comparators = []

    def create_comparator(self, reference_name):
        subset_paf = self.paf[self.paf["target_name"] == reference_name]
        subset_tsv = self.tsv[self.tsv[0] == reference_name]  # first column is name
        alignment = {"paf": subset_paf, "tsv": subset_tsv}

        comparator = Comparator({reference_name: self.references[reference_name]}, alignment)

        return comparator

    @staticmethod
    def load_tsv(tsv_path):
        """Create a dataframe from a tsv file of depth counts.


        **Parameters**

        **tsv_path**
        > Path to output of `samtools depth -aa` (`str`).
        """
        tsv = pandas.read_csv(tsv_path, sep="\t", header=None)
        tsv.columns = ["name", "position", "depth"]

        return tsv


class Comparator:
    """Analyse alignments to a reference.


    **Parameters**

    **reference**
    > Dictionary of one name: one Biopython record.

    **alignment**
    > Dictionary of pandas dataframes: {"paf": df_paf, "tsv": df_tsv}. The dataframes
    contain entries only for the reference.
    """
    def __init__(self, reference, alignment):
        self.reference = reference
        self.paf = alignment["paf"]
        self.tsv = alignment["tsv"]
