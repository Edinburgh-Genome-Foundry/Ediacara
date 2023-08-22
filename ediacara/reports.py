# Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of Ediacara.
#
# Ediacara is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Ediacara is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

from datetime import datetime
import os

import matplotlib.pyplot as plt

from pdf_reports import (
    add_css_class,
    dataframe_to_html,
    pug_to_html,
    style_table_rows,
    write_report,
)
import pdf_reports.tools as pdf_tools

from .version import __version__
from .Comparator import result_keywords

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
SEQUENCINGGROUP_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "run_simulation_report.pug")
ASSEMBLY_ANALYSIS_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "assembly_report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, "report_style.css")


def end_pug_to_html(template, **context):
    now = datetime.now().strftime("%Y-%m-%d")
    defaults = {
        "sidebar_text": "Generated on %s by Ediacara (version %s)" % (now, __version__),
        "ediacara_logo_url": os.path.join(ASSETS_PATH, "imgs", "logo.png"),
    }
    for k in defaults:
        if k not in context:
            context[k] = defaults[k]
    return pug_to_html(template, **context)


def write_assembly_analysis_report(target, assemblybatch):
    """Write a sequencing run report with a PDF summary.


    **Parameters**

    **target**
    > Path for PDF file.

    **assemblybatch**
    > `AssemblyBatch` instance.
    """
    for assembly in assemblybatch.assemblies:
        assembly.assembly_figure_data = pdf_tools.figure_data(
            assembly.assembly_figure, fmt="svg"
        )

        assembly.report_table = dataframe_to_html(
            assembly.subset_paf(), extra_classes=("definition",)
        )

    html = end_pug_to_html(
        ASSEMBLY_ANALYSIS_REPORT_TEMPLATE, assemblybatch=assemblybatch
    )
    write_report(html, target, extra_stylesheets=(STYLESHEET,))


def write_sequencinggroup_report(target, sequencinggroup):
    """Write a sequencing run report with a PDF summary.


    **Parameters**

    **target**
    > Path for PDF file.

    **sequencinggroup**
    > `SequencingGroup` instance.
    """

    def tr_modifier(tr):
        tds = list(tr.find_all("td"))
        if len(tds) == 0:
            return
        result = tds[1]  # second element of list is the result symbol
        if result.text == result_keywords["good"]:
            add_css_class(tr, "positive")
        elif result.text == result_keywords["error"]:
            add_css_class(tr, "negative")

    if not sequencinggroup.comparisons_performed:
        return "Run perform_all_comparisons_in_sequencinggroup()!"

    for comparatorgroup in sequencinggroup.comparatorgroups:
        comparatorgroup.report_table = dataframe_to_html(
            comparatorgroup.summary_table, extra_classes=("definition",)
        )

        # This colours the summary table:
        comparatorgroup.report_table = style_table_rows(
            comparatorgroup.report_table, tr_modifier
        )

        # Histogram of reads in the report summary
        histogram_height = comparatorgroup.fastq_plot.get_size_inches()[1]
        comparatorgroup.fastq_figure_data = pdf_tools.figure_data(
            comparatorgroup.fastq_plot, fmt="svg", size=[8, histogram_height]
        )
        for comparator in comparatorgroup.comparators:
            comparator.figure_data = pdf_tools.figure_data(comparator.fig, fmt="svg")
            comparator.insert_plot_data = pdf_tools.figure_data(
                comparator.insert_plot, fmt="svg"
            )

            # Keep first few only:
            if comparator.vcf_table.shape[0] > sequencinggroup.vcf_cutoff:  # 0 for rows
                vcf_table = comparator.vcf_table[: sequencinggroup.vcf_cutoff]
                comparator.vcf_table_message = True
            else:
                vcf_table = comparator.vcf_table
                comparator.vcf_table_message = False

            # This column has lists, does not look good in print:
            vcf_table["ALT"] = vcf_table["ALT"].apply(lambda x: ", ".join(x))
            # Convert to string then truncate too long entries and add three dots:
            columnwidth = 15  # wider columns don't display well in PDF
            vcf_table = vcf_table.apply(
                lambda x: x.astype(str).apply(
                    lambda y: y[:columnwidth] + "..." if len(y) > columnwidth else y
                )
            )

            def tr_modifier_for_vcf_table(tr):
                tds = list(tr.find_all("td"))
                if len(tds) == 0:
                    return
                true_positive = tds[-1]  # last column shows call
                if true_positive.text == "1":
                    add_css_class(tr, "negative")
                # else:
                #     add_css_class(tr, "positive")

            comparator.vcf_table_html = dataframe_to_html(
                vcf_table, extra_classes=("definition",)
            )
            comparator.vcf_table_html = style_table_rows(
                comparator.vcf_table_html, tr_modifier_for_vcf_table
            )

            if hasattr(comparator, "is_comparison_successful"):
                if comparator.is_comparison_successful:
                    height = comparator.comparison_figure.figure.get_size_inches()[1]
                    if height > 10:
                        height = 10  # to fit on one page
                    comparator.comparison_figure_data = pdf_tools.figure_data(
                        comparator.comparison_figure, fmt="svg", size=[7, height]
                    )
                else:
                    comparator.comparison_figure_data = None

            if comparator.has_consensus:
                comparator.has_comparison_error = True
                if comparator.geneblocks_outcome == "none":
                    comparator.geneblocks_text = (
                        "Missing consensus file for comparison!"
                    )
                elif comparator.geneblocks_outcome == "incorrect_length":
                    comparator.geneblocks_text = (
                        "Incorrect length! " + comparator.incorrect_length_msg
                    )
                elif comparator.geneblocks_outcome == "too_many_differences":
                    comparator.geneblocks_text = (
                        "Plot not created. "
                        + "Too many differences between consensus and reference."
                    )

                elif comparator.geneblocks_outcome == "geneblocks_error":
                    comparator.geneblocks_text = (
                        "GeneBlocks comparison of consensus and reference failed."
                    )
                elif comparator.geneblocks_outcome == "swapped_diffblocks":
                    comparator.geneblocks_text = (
                        "Note: plotting the reference against the consensus failed."
                        + " Instead, this plot compares the consensus to the reference "
                        + comparator.name
                        + " , therefore there are no annotations."
                    )
                    comparator.has_comparison_error = False
                elif comparator.geneblocks_outcome == "all_good":
                    comparator.geneblocks_text = (
                        "<b>"
                        + comparator.name
                        + "</b> reference vs consensus of reads:"
                    )
                    comparator.has_comparison_error = False

        # plt.close("all")

    html = end_pug_to_html(
        SEQUENCINGGROUP_REPORT_TEMPLATE, sequencinggroup=sequencinggroup
    )
    write_report(html, target, extra_stylesheets=(STYLESHEET,))
