# AlphaFold2 structure dataset

Dataset composed of only AlphaFold2 structures

Created on: 2022-03-12

This dataset is licensed under CC BY 4.0. (https://creativecommons.org/licenses/by/4.0/)


## Directory structure

```txt
.
├── fasta
│   ├── target1.fasta
│   └── target2.fasta
├── native_pdb
│   ├── target1.pdb
│   └── target2.pdb
├── pdb
│   ├── target1
│   └── target2
└── score
    ├── all_score.csv
    └── target_list.csv
```

## File details

* fasta

    Target protein sequence in fasta format

* native_pdb

    Native structure in pdb format

* pdb

    Predicted structures for each target in pdb format

* score

    * all_score.csv

        csv file containing all scores such as labels for each predicted structure

    * target_list.csv

        csv file containing target information such as resolution, sequence length, and domain number.
