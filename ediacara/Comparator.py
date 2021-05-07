import statistics

import matplotlib.pyplot as plt
import numpy
import pandas

import Bio

import dna_features_viewer


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
        subset_tsv = self.tsv[self.tsv["name"] == reference_name]
        alignment = {"paf": subset_paf, "tsv": subset_tsv}

        comparator = Comparator(
            {reference_name: self.references[reference_name]}, alignment
        )

        return comparator

    @staticmethod
    def load_paf(paf_path):
        """Create a dataframe from a paf file of alignments.


        **Parameters**

        **paf_path**
        > Path (`str`) to PAF file: `paftools.js sam2paf aln.sam > aln.paf`.
        """
        paf = pandas.read_csv(paf_path, sep="\t", header=None)

        # First 12 columns are defined by the format:
        columns_12 = [
            "query_name",
            "query_length",
            "query_start",
            "query_end",
            "strand",
            "target_name",
            "target_length",
            "target_start",
            "target_end",
            "mapping_matches",
            "mapping_size",
            "mapping_quality",
        ]
        # For additional columns optionally created during alignment:
        columns_13plus = [str(index) for index in range(12, len(paf.columns))]
        paf.columns = columns_12 + columns_13plus

        return paf

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
        self.name = next(iter(reference))
        self.record = next(iter(reference.values()))
        self.paf = alignment["paf"]
        self.tsv = alignment["tsv"]

    def plot_coverage(self):
        """Plot the reference with the coverage and weak reads."""

        xx = numpy.arange(len(self.record.seq))  # for the plot x axis
        yy = self.tsv["depth"].to_list()  # for plotting coverage
        if not hasattr(self, "zz"):
            self.get_weak_read_histogram_data(cutoff=0.8)
            # for plotting weak reads (less than 80% of the read maps)

        # Plot
        ######
        fig, (ax1, ax2, ax3) = plt.subplots(
            3, 1, figsize=(12, 6), sharex=True, gridspec_kw={"height_ratios": [4, 1, 1]}
        )
        graphic_record = dna_features_viewer.BiopythonTranslator().translate_record(
            self.record
        )
        graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)

        # Plot coverage
        ax2.fill_between(xx, yy, alpha=0.8)
        median_yy = statistics.median(yy)
        stdev_yy = statistics.stdev(yy)
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel("Coverage [x]")
        ax2.set_xlabel("Base [bp]")
        ax2.axhline(
            y=median_yy, xmin=0, xmax=1, color="grey", linestyle="-", linewidth=0.8
        )
        ax2.axhline(
            y=median_yy - stdev_yy,
            xmin=0,
            xmax=1,
            color="grey",
            linestyle="--",
            linewidth=0.8,
        )
        ax2.axhline(
            y=median_yy + stdev_yy,
            xmin=0,
            xmax=1,
            color="grey",
            linestyle="--",
            linewidth=0.8,
        )

        # Plot low-quality reads:
        ax3.fill_between(xx, self.zz, alpha=0.8, color="red")
        ax3.set_ylim(bottom=0)
        ax3.set_ylabel("Weak reads")
        ax3.set_xlabel("Base [bp]")
        ax3.axhline(
            y=median_yy * 0.15,  # 15% of max coverage as guide
            xmin=0,
            xmax=1,
            color="grey",
            linestyle="--",
            linewidth=0.8,
        )

    def get_weak_read_histogram_data(self, cutoff):
        """Create data for the plotting function.


        **Parameters**

        **cutoff**
        > (`float`). Recommended: 0.8.
        """
        # For each read, calculate the proportion of bases that map, then create an array
        # for each read that is below the cutoff. The sum of these arrays will give the
        # depth of aligning bad reads for each base of the reference.

        self.proportion_maps = (
            self.paf["mapping_matches"] / self.paf["query_length"].to_list()
        )
        filtered_paf = self.paf[self.proportion_maps < cutoff]
        array_length = filtered_paf.iloc[0]["target_length"]
        histogram_array = []

        for index, row in filtered_paf.iterrows():
            read_array = [0] * array_length
            start = row["target_start"]
            end = row["target_end"]
            values = [1] * (end - start)
            read_array[start:end] = values
            histogram_array.append(read_array)

        histogram_numpy_array = numpy.array(histogram_array)

        self.zz = numpy.sum(histogram_numpy_array, axis=0)
