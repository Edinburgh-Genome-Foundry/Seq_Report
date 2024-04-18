import csv

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
    > The currency symbol (or string) to display in the report (`str`).

    **projectname**
    > The name of the project (`str`).

    **comments**
    > Any comments to be included in the report (`str`).

    **min_length**
    > Check that all sequences are at least this long (`int`).

    **max_length**
    > Check that all sequences are at most this long (`int`).

    **name_length**
    > Check which sequence IDs are longer than this cutoff. Genbank has a character
    limit.
    """

    def __init__(
        self,
        records,
        cost_per_base=0.25,
        cost_per_seq=0,
        currency_symbol="£",
        projectname="",
        comments="",
        min_length=0,
        max_length=0,
        name_length=15,  # max seq record name length. Genbank character limit.
    ):
        self.sequences = records
        self.cost_per_base = cost_per_base
        self.cost_per_seq = cost_per_seq
        self.currency_symbol = currency_symbol
        self.projectname = projectname
        self.comments = comments
        self.min_length = min_length
        self.max_length = max_length
        self.name_length = name_length
        self.calculate_values()

    def calculate_values(self):
        self.cost_per_base = float(self.cost_per_base)  # could be str after update
        self.cost_per_seq = float(self.cost_per_seq)
        self.n_seq = len(self.sequences)

        n_bp = 0
        for part in self.sequences:
            n_bp += len(part.seq)
        self.n_bp = n_bp
        self.cost = self.n_seq * self.cost_per_seq + self.n_bp * self.cost_per_base

        self.too_short = []
        self.too_long = []
        self.long_names = []
        for record in self.sequences:
            if len(record) < self.min_length:
                self.too_short += [record.id]
            if self.max_length > 0:  # otherwise skip if default
                if len(record) > self.max_length:
                    self.too_long += [record.id]
            if len(record.id) > self.name_length:
                self.long_names += [record.id]


def read_fasta(fasta):
    """Read a FASTA sequence file into a list of records.


    **Parameters**

    **fasta**
    > The FASTA filepath (`str`).
    """
    return list(SeqIO.parse(fasta, "fasta"))


def seqcollection_from_csv(csv_file, records=None, param_dict={}):
    """Create a SeqCollection, using parameters in a CSV file.

    The CSV file parameters override the default class parameters, and this function's
    parameters override the CSV file parameters. Either a FASTA file (in the CSV or in
    `param_dict`) or a list of `SeqRecords` must be specified.
    For the parameter descriptions, see the docstring of `SeqCollection`.
    """
    with open(csv_file, "r") as f:
        reader = csv.reader(f, skipinitialspace=True)
        next(reader)  # ignore header
        parameters = {}
        for row in reader:
            if row[1] == "":  # in case the full param list is in the csv
                continue
            else:
                parameters[row[0]] = row[1]

    for key, value in param_dict.items():  # override with function values
        parameters[key] = value
    if records:
        pass
    else:
        records = read_fasta(parameters["fasta"])

    seq_coll = SeqCollection(records=records)
    for key, value in parameters.items():
        setattr(seq_coll, key, value)
    seq_coll.calculate_values()
    return seq_coll
