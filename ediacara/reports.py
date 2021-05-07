from datetime import datetime
import os

import pandas

from pdf_reports import (
    dataframe_to_html,
    pug_to_html,
    write_report,
)

from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "simulation_report.pug")
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


def write_pdf_report(target, comparator, csv_path=None):
    """Write an alignment report with a PDF summary.


    **Parameters**

    **target**
    > Path for PDF file.

    **comparator**
    > Comparator instance.

    **csv_path**
    > Path to CSV output of results (`str`). Default no CSV.
    """
    html = end_pug_to_html(
        REPORT_TEMPLATE,
        id=comparator.record.id,
        reference_length=len(comparator.record),
    )

    # if csv_path is not None:
    #     analysis_summary.to_csv(csv_path, index=False)

    write_report(html, target, extra_stylesheets=(STYLESHEET,))
