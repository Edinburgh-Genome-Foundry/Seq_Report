from Bio import SeqIO


class SeqCollection:
    """Class to store basic properties of a collection of sequences.


    **Parameters**

    **records**
    > A list of Biopython SeqRecords.

    **cost_per_base**
    > Cost per nucleotide base.

    **cost_per_seq**
    > Fix overhead cost for each sequence part (cloning, delivery costs etc).

    **currency_symbol**
    > The currency symbol to display in the report.

    **projectname**
    > The name of the project (`str`).

    **comments**
    > Any comments to be included in the report (`str`).
    """

    def __init__(
        self,
        records,
        cost_per_base=0.25,
        cost_per_seq=0,
        currency_symbol="£",
        projectname="",
        comments="",
    ):
        self.sequences = records
        self.cost_per_base = cost_per_base
        self.cost_per_seq = cost_per_seq
        self.currency_symbol = currency_symbol
        self.n_seq = len(self.sequences)
        n_bp = 0
        for part in self.sequences:
            n_bp += len(part.seq)
        self.n_bp = n_bp
        self.cost = self.n_seq * self.cost_per_seq + self.n_bp * self.cost_per_base
        self.projectname = projectname
        self.comments = comments


def read_fasta(fasta):
    """Read a FASTA sequence file into a list of records.


    **Parameters**

    **fasta**
    > The FASTA filepath (`str`).
    """
    return list(SeqIO.parse(fasta, "fasta"))
