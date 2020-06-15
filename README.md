# Installation

The best way to manage the dependency of this repo is through a conda virtual
environment. Use the following commands in your terminal

```
conda create intactness python=3
conda activate intactness
conda install --file requirements.txt
```

Please check this [site](https://www.anaconda.com/) for information.


# Running instruction

The whole repo is intended to be used as a template folder for reproducible
research purpose. Please create a copy of this folder each time your run the
pipeline.

The following steps should be taken in order, in the root folder of the project.

1. Store the DNA sequences in a FASTA file named "seqs.fasta" in `data` folder
1. Specify your email address on the 5th line in `intactness/default.cfg`
1. `conda activate intactness`
1. `python3 -m intactness`
1. The result is in `data/seqs/summary.csv`

For advanced users, the settings are stored in `intactness/default.cfg`.
