import os

import seqreport

data_dir = os.path.join("tests", "data")
seq_fasta = os.path.join(data_dir, "test.fa")


def test_SeqCollection(tmpdir):
    seq_coll = seqreport.SeqCollection(
        fasta=seq_fasta, projectname="EGF24", comments="This is a test sequence set."
    )
    assert seq_coll.n_seq == 3
    assert seq_coll.n_bp == 99
    assert seq_coll.projectname == "EGF24"
    assert seq_coll != ""
