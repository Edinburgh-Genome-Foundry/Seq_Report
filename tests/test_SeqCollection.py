import os

import seqreport

data_dir = os.path.join("tests", "data")
seq_fasta = os.path.join(data_dir, "test.fa")


def test_SeqCollection():
    seq_records = seqreport.read_fasta(seq_fasta)
    seq_coll = seqreport.SeqCollection(
        records=seq_records,
        projectname="EGF24",
        comments="This is a test sequence set.",
        min_length=20,
        max_length=40,
        name_length=10,
    )
    assert seq_coll.n_seq == 6
    assert seq_coll.n_bp == 175
    assert seq_coll.projectname == "EGF24"
    assert seq_coll.comments != ""
    assert len(seq_coll.too_short) == 1
    assert len(seq_coll.too_long) == 1
    assert len(seq_coll.long_names) == 1
    assert len(seq_coll.repeat_names) == 1
    assert len(seq_coll.repeat_seq) == 1
    assert len(seq_coll.reverse_complement_seq) == 2


def test_read_fasta():
    seq_records = seqreport.read_fasta(seq_fasta)
    assert len(seq_records) == 6


def test_seqcollection_from_csv():
    csv_path = os.path.join(data_dir, "values.csv")
    seq_coll = seqreport.seqcollection_from_csv(csv_file=csv_path)
    assert seq_coll.fasta_name == "test.fa"  # see in CSV file
    # (not tested passing records or param_dict above)
