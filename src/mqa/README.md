# MQA

## Directory structure

```
data/out/dataset
├── alphafold_output
├── fasta
├── native_pdb
├── pdb
└── score
    ├── label
    ├── mqa
    │   └── ProQ3D
    └── subsets
        └── subset_A
            ├── ProQ3D.csv
            └── label.csv
```

## Excecution

```bash
usage: throw_mqa_jobs.py [-h] -m {DOPE,ProQ3D,SBROD,P3CMQA,DeepAccNet,VoroCNN}
                         -t TARGET_LIST_CSV [-q] [-e]

Throw MQA job

optional arguments:
  -h, --help            show this help message and exit
  -m {DOPE,ProQ3D,SBROD,P3CMQA,DeepAccNet,VoroCNN}, --method {DOPE,ProQ3D,SBROD,P3CMQA,DeepAccNet,VoroCNN}
                        MQA method name
  -t TARGET_LIST_CSV, --target_list_csv TARGET_LIST_CSV
                        target list csv
  -q, --qsub            throw job using qsub
  -e, --execute         execute the job
```

```bash
python throw_mqa.py -m method -t target_subset_list.csv --qsub --execute
```

If execution for all targets are completed, all scores are concatenated.


## Calculate MQA scores only for resolved residues

For MQA methods (except pLDDT and pTM)
```bash
$ python calc_score_resolved_region.py --help
usage: calc_score_resolved_region.py [-h] target_csv method

positional arguments:
  target_csv  Path to target csv file
  method      Method name

optional arguments:
  -h, --help  show this help message and exit
```

For pLDDT and pTM,
```bash
$ cd ../alphafold/util
$ ../alphafold/colabfold-conda/bin/python3.7 calc_confidence_resolved_region.py target_csv_path
```
