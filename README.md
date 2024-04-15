<p align="center">
<img alt="EGF logo" title="EGF" src="images/egf.png" width="120">
</p>

# Seq Report

<!-- [![build](https://github.com/Edinburgh-Genome-Foundry/Seq_Report/actions/workflows/build.yml/badge.svg)](https://github.com/Edinburgh-Genome-Foundry/Seq_Report/actions/workflows/build.yml) -->

<!-- [![coverage](https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/Seq_Report/badge.svg?branch=main)](https://coveralls.io/github/Edinburgh-Genome-Foundry/Seq_Report?branch=main) -->

Simple reporting on a set of sequences for documentation purposes.

## Install

```bash
pip install git+https://github.com/Edinburgh-Genome-Foundry/Seq_Report.git
```

## Usage

```python
import seqreport

seq_fasta = "seq.fa"
seq_coll = seqreport.SeqCollection(records=seqreport.read_fasta(seq_fasta), projectname="EGF24")
seqreport.write_pdf_report("seq_report.pdf", seq_coll)
```

<p align="center">
<img alt="Seq Report" title="Seq Report" src="images/seqreport_screenshot.png" width="360">
</p>

Alternatively, use a CSV file to specify parameters:

```python
seq_coll = seqreport.seqcollection_from_csv(csv_file="tests/data/values.csv")
```

Header is ignored, but there must be a header. Entries must match the SeqCollection parameters ([example](tests/data/values.csv)).

## Versioning

Seq Report uses the [semantic versioning](https://semver.org) scheme.

## License = MIT

Seq Report is free/libre and open-source software, which means the users have the freedom to run, study, change and distribute the software.

Seq Report was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp).

Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh
