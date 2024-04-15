import dnacauldron


class SeqCollection:
    """Class to store basic properties of a collection of sequences.


    **Parameters**

    **fasta**
    > The FASTA file of the sequences.

    **cost_per_base**
    > Cost per nucleotide base.

    **cost_per_seq**
    > Fix overhead cost for each sequence part (cloning, delivery costs etc).

    **currency_symbol**
    > The currency symbol to display in the report.

    **projectname**
    > The name of the project (`str`).
    """

    def __init__(
        self,
        fasta,
        cost_per_base=0.25,
        cost_per_seq=0,
        currency_symbol="£",
        projectname="",
    ):
        self.fasta = fasta
        self.cost_per_base = cost_per_base
        self.cost_per_seq = cost_per_seq
        self.currency_symbol = currency_symbol
        self.sequences = dnacauldron.biotools.load_records_from_files(
            files=[self.fasta], use_file_names_as_ids=False
        )
        self.n_seq = len(self.sequences)
        n_bp = 0
        for part in self.sequences:
            n_bp += len(part.seq)
        self.n_bp = n_bp
        self.cost = self.n_seq * self.cost_per_seq + self.n_bp * self.cost_per_base
        self.projectname = projectname
