from datetime import datetime
import os

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
SEQCOLLECTION_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "seq_report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, "report_style.css")


def end_pug_to_html(template, **context):
    now = datetime.now().strftime("%Y-%m-%d")
    defaults = {
        "sidebar_text": "Generated on %s by Seq Report (version %s)"
        % (now, __version__),
    }
    for k in defaults:
        if k not in context:
            context[k] = defaults[k]
    return pug_to_html(template, **context)


def write_pdf_report(target, seqcollection):
    """Write a sequence collection report with a PDF summary.


    **Parameters**

    **target**
    > Path for PDF file.

    **seqcollection**
    > `SeqCollection` instance.
    """

    html = end_pug_to_html(SEQCOLLECTION_REPORT_TEMPLATE, seqcollection=seqcollection)
    write_report(html, target, extra_stylesheets=(STYLESHEET,))
