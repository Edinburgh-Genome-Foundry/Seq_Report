import csv
import os

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
    limit (`int`).

    **assembly_plan**
    > Optional assembly plan CSV path for calculating savings by re-using DNA parts.

    **fasta_name**
    > Optional name of the fasta file for easy identification (`str`).
    """

    def __init__(
        self,
        records,
        cost_per_base=0.25,
        cost_per_seq=0,
        currency_symbol="£",
        projectname="",
        comments="",
        min_length=100,  # a good cutoff for min DNA synthesis length
        max_length=3000,  # a good cutoff for max DNA synthesis length
        name_length=15,  # max seq record name length. Genbank character limit.
        assembly_plan="",
        fasta_name="",
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
        self.assembly_plan = assembly_plan
        self.fasta_name = fasta_name
        self.calculate_values()

    def calculate_values(self):
        self.cost_per_base = float(self.cost_per_base)  # could be str after update
        self.cost_per_seq = float(self.cost_per_seq)
        self.n_seq = len(self.sequences)

        # Numbers section
        n_bp = 0
        for part in self.sequences:
            n_bp += len(part.seq)
        self.n_bp = n_bp
        self.cost = self.n_seq * self.cost_per_seq + self.n_bp * self.cost_per_base
        self.cost = round(self.cost)  # (in)accuracy is fine for our purposes

        # Lengths section
        self.too_short = []
        self.too_long = []
        self.long_names = []
        for record in self.sequences:
            if len(record) < self.min_length:
                self.too_short += [record.id]
            if len(record) > self.max_length:
                self.too_long += [record.id]
            if len(record.id) > self.name_length:
                self.long_names += [record.id]

        self.too_short = list(set(self.too_short))
        self.too_long = list(set(self.too_long))
        self.long_names = list(set(self.long_names))
        # Repeats will be checked further below.

        # For the PDF report:
        self.too_short_text = " ; ".join(self.too_short)
        self.too_long_text = " ; ".join(self.too_long)
        self.long_names_text = " ; ".join(self.long_names)

        # Repeats section
        self.repeat_names = []
        self.repeat_seq = []
        self.reverse_complement_seq = []
        for index, record in enumerate(self.sequences[:-1]):
            for other_record in self.sequences[index+1:]:
                if record.id == other_record.id:
                    self.repeat_names += [record.id]
                if record.seq == other_record.seq:
                    self.repeat_seq += [(record.id, other_record.id)]
                if record.seq == other_record.seq.reverse_complement():
                    self.reverse_complement_seq += [(record.id, other_record.id)]
                    
        self.repeat_names = list(set(self.repeat_names))

        # For the PDF report:
        self.repeat_names_text = " ; ".join(self.repeat_names)

        repeat_seq_text_list = []
        for identifier in self.repeat_seq:
            combined = " = ".join(identifier)
            repeat_seq_text_list += [combined]
        self.repeat_seq_text = " ; ".join(repeat_seq_text_list)

        reverse_complement_seq_text_list = []
        for identifier in self.reverse_complement_seq:
            combined = " = ".join(identifier)
            reverse_complement_seq_text_list += [combined]
        self.reverse_complement_seq_text = " ; ".join(reverse_complement_seq_text_list)

        # Savings section
        if self.assembly_plan:
            all_rows = []
            with open(self.assembly_plan, "r") as f:
                reader = csv.reader(f, skipinitialspace=True)
                next(reader)  # ignore header
                for row in reader:
                    all_rows += row[1:]  # first column is construct name
            self.savings_list = []
            self.total_savings = 0
            for record in self.sequences:
                count_in_plan = all_rows.count(record.id)
                if count_in_plan > 1:  # we count re-use savings only
                    self.total_savings += len(record.seq) * (count_in_plan - 1)  # ignore first synthesis
                    self.savings_list += [record.id]
            self.total_cost_savings = self.total_savings * self.cost_per_base  # ignore cost / seq
            self.total_cost_savings = round(self.total_cost_savings)  # (in)accuracy is fine for our purposes
            # For the PDF report:
            self.savings_list_text = " ; ".join(self.savings_list)
        else:
            self.savings_list_text = ""


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
    For the parameter descriptions, see parameters in the docstring of `SeqCollection`,
    except records and fasta_name.
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

    # To ensure that the correct filename is used:
    setattr(seq_coll, "fasta_name", os.path.basename(parameters["fasta"]))
    seq_coll.calculate_values()

    return seq_coll
