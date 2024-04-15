import os

import seqreport

data_dir = os.path.join("tests", "data")
seq_fasta = os.path.join(data_dir, "test.fa")


def test_write_pdf_report(tmpdir):
    seq_records = seqreport.read_fasta(seq_fasta)
    seq_coll = seqreport.SeqCollection(records=seq_records)
    pdf_path = os.path.join(str(tmpdir), "test_report.pdf")
    seqreport.write_pdf_report(target=pdf_path, seqcollection=seq_coll)

    with open(pdf_path, "rb") as f:
        filesize = len(f.read())
        assert filesize > 40000
