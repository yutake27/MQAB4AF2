# AlphaFold2 structure dataset

Dataset composed of only AlphaFold2 structures

Created on: 2022-03-12

This dataset is licensed under CC BY 4.0. (https://creativecommons.org/licenses/by/4.0/)

[AlphaFold](https://github.com/deepmind/alphafold) v2.0.1, [ColabFold](https://github.com/sokrypton/ColabFold), and [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold) v1.0.0 were used to predict structures.

- Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, Tunyasuvunakool K, Bates R, Žídek A, Potapenko A and Bridgland A. "Highly accurate protein structure prediction with AlphaFold." Nature. 2021 Aug;596(7873):583-9. doi: [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)

- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. "ColabFold - Making protein folding accessible to all." bioRxiv (2021) doi: [10.1101/2021.08.15.456425](https://www.biorxiv.org/content/10.1101/2021.08.15.456425v2)

- Moriwaki Y. "LocalColabFold" Available from: https://github.com/YoshitakaMo/localcolabfold.

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

- fasta

    Target protein sequence in fasta format

- native_pdb

    Native structure in pdb format

- pdb

    Predicted structures for each target in pdb format

- score

    - all_score.csv

        csv file containing all scores such as labels for each predicted structure

    - target_list.csv

        csv file containing target information such as resolution, sequence length, and domain number.
