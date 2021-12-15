from datetime import datetime
import os

import matplotlib.pyplot as plt
import pandas

from pdf_reports import (
    add_css_class,
    dataframe_to_html,
    pug_to_html,
    style_table_rows,
    write_report,
)
import pdf_reports.tools as pdf_tools

from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
# REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "simulation_report.pug")
# GROUP_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "group_simulation_report.pug")
SEQUENCINGGROUP_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "run_simulation_report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, "report_style.css")


def end_pug_to_html(template, **context):
    now = datetime.now().strftime("%Y-%m-%d")
    defaults = {
        "sidebar_text": "Generated on %s by EGF's alignment comparator (version %s)"
        % (now, __version__),
        "end_logo_url": os.path.join(ASSETS_PATH, "imgs", "logo.png"),
    }
    for k in defaults:
        if k not in context:
            context[k] = defaults[k]
    return pug_to_html(template, **context)


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
        if result.text == "☑":
            add_css_class(tr, "positive")
        elif result.text == "☒":
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
        comparatorgroup.fastq_figure_data = pdf_tools.figure_data(
            comparatorgroup.fastq_plot, fmt="svg"
        )
        for comparator in comparatorgroup.comparators:
            comparator.figure_data = pdf_tools.figure_data(comparator.fig, fmt="svg")

            # Keep first few only:
            if comparator.vcf_table.shape[0] > sequencinggroup.vcf_cutoff:  # 0 for rows
                vcf_table = comparator.vcf_table[: sequencinggroup.vcf_cutoff]
                comparator.vcf_table_message = True
            else:
                vcf_table = comparator.vcf_table
                comparator.vcf_table_message = False

            # Convert to string then truncate too long entries and add three dots:
            columnwidth = 15  # wider columns don't display well in PDF
            vcf_table = vcf_table.apply(
                lambda x: x.astype(str).apply(
                    lambda y: y[:columnwidth] + "..." if len(y) > columnwidth else y
                )
            )

            comparator.vcf_table_html = dataframe_to_html(
                vcf_table, extra_classes=("definition",)
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

            if comparator.has_de_novo:
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
                        "Note: the plot compares the consensus to the reference "
                        + comparator.name
                        + " therefore there are no annotations."
                    )
                    comparator.has_comparison_error = False
                elif comparator.geneblocks_outcome == "all_good":
                    comparator.geneblocks_text = (
                        "<b>"
                        + comparator.name
                        + "</b> reference vs consensus of reads:"
                    )
                    comparator.has_comparison_error = False

                if (
                    hasattr(comparator, "is_assembly_reverse_complement")
                    and comparator.is_assembly_reverse_complement
                ):
                    comparator.geneblocks_text += (
                        "Note: the consensus is the "
                        "reverse complement of the reference."
                    )
        plt.close("all")

    html = end_pug_to_html(
        SEQUENCINGGROUP_REPORT_TEMPLATE, sequencinggroup=sequencinggroup
    )
    write_report(html, target, extra_stylesheets=(STYLESHEET,))


# def write_comparatorgroup_report(target, comparatorgroup):
#     """Write an alignment report with a PDF summary.


#     **Parameters**

#     **target**
#     > Path for PDF file.

#     **comparatorgroup**
#     > `ComparatorGroup` instance.
#     """
#     if not comparatorgroup.comparisons_performed:
#         return "Run perform_all_comparisons()!"

#     comparatorgroup.report_table = dataframe_to_html(
#         comparatorgroup.summary_table, extra_classes=("definition",)
#     )

#     def tr_modifier(tr):
#         tds = list(tr.find_all("td"))
#         if len(tds) == 0:
#             return
#         result = tds[1]  # second element of list is the result symbol
#         if result.text == "☑":
#             add_css_class(tr, "positive")
#         elif result.text == "☒":
#             add_css_class(tr, "negative")

#     # This colours the summary table:
#     comparatorgroup.report_table = style_table_rows(
#         comparatorgroup.report_table, tr_modifier
#     )

#     # Histogram of reads in the report summary
#     comparatorgroup.fastq_figure_data = pdf_tools.figure_data(
#         comparatorgroup.fastq_plot, fmt="svg"
#     )

#     for comparator in comparatorgroup.comparators:
#         comparator.figure_data = pdf_tools.figure_data(comparator.fig, fmt="svg")

#         comparator.vcf_table_html = dataframe_to_html(
#             comparator.vcf_table, extra_classes=("definition",)
#         )

#         if hasattr(comparator, "is_comparison_successful"):
#             if comparator.is_comparison_successful:
#                 height = comparator.comparison_figure.figure.get_size_inches()[1]
#                 if height > 10:
#                     height = 10  # to fit on one page
#                 comparator.comparison_figure_data = pdf_tools.figure_data(
#                     comparator.comparison_figure, fmt="svg", size=[7, height]
#                 )
#             else:
#                 comparator.comparison_figure_data = None

#             if comparator.has_de_novo:
#                 comparator.has_comparison_error = True
#                 if comparator.geneblocks_outcome == "none":
#                     comparator.geneblocks_text = (
#                         "Missing consensus file for comparison!"
#                     )
#                 elif comparator.geneblocks_outcome == "incorrect_length":
#                     comparator.geneblocks_text = (
#                         "Incorrect length! " + comparator.incorrect_length_msg
#                     )
#                 elif comparator.geneblocks_outcome == "geneblocks_error":
#                     comparator.geneblocks_text = (
#                         "GeneBlocks comparison of consensus and reference failed."
#                     )
#                 elif comparator.geneblocks_outcome == "swapped_diffblocks":
#                     comparator.geneblocks_text = (
#                         "Note: the plot compares the consensus to the reference "
#                         + comparator.name
#                         + " therefore there are no annotations."
#                     )
#                     comparator.has_comparison_error = False
#                 elif comparator.geneblocks_outcome == "all_good":
#                     comparator.geneblocks_text = (
#                         "<b>"
#                         + comparator.name
#                         + "</b> reference vs consensus of reads:"
#                     )
#                     comparator.has_comparison_error = False

#             if (
#                 hasattr(comparator, "is_assembly_reverse_complement")
#                 and comparator.is_assembly_reverse_complement
#             ):
#                 comparator.geneblocks_text += (
#                     "Note: the consensus is the reverse complement of the reference."
#                 )

#     html = end_pug_to_html(GROUP_REPORT_TEMPLATE, comparatorgroup=comparatorgroup)
#     write_report(html, target, extra_stylesheets=(STYLESHEET,))


# def write_pdf_report(
#     target, comparator, assembly_path, csv_path=None,
# ):
#     """Write an alignment report with a PDF summary.


#     **Parameters**

#     **target**
#     > Path for PDF file.

#     **comparator**
#     > Comparator instance.

#     **assembly_path**
#     > Assembly to compare with. See `Comparator.compare_with_assembly()`.

#     **csv_path**
#     > Path to CSV output of results (`str`). Default no CSV.
#     """
#     comparator.perform_comparison(assembly_path)
#     comparator.figure_data = pdf_tools.figure_data(comparator.fig, fmt="svg")

#     if hasattr(comparator, "is_comparison_successful"):
#         if comparator.is_comparison_successful:
#             comparator.has_comparison = True
#             comparator.comparison_figure_data = pdf_tools.figure_data(
#                 comparator.comparison_figure, fmt="svg"
#             )
#         else:
#             comparator.comparison_figure_data = None

#     html = end_pug_to_html(REPORT_TEMPLATE, comparator=comparator,)

#     # if csv_path is not None:
#     #     analysis_summary.to_csv(csv_path, index=False)

#     write_report(html, target, extra_stylesheets=(STYLESHEET,))
