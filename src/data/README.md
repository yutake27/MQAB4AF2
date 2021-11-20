# Target selection
## Download necessary files
```
bash download.sh
```
PDB sequence file and a file with clustered PDB entries will be downloaded.

## Select targets
```
python select_target.py
```
Target list will be created at `../../data/interim/target_list.csv`.

## Select subsets
```
> python select_subsets.py
usage: select_subset.py [-h] [-i INPUT_CSV_PATH] [--how {random,eq_random,head,tail}] [--target_num TARGET_NUM] [--random_state RANDOM_STATE] [-o OUTPUT_CSV_PATH]

Select subsets of the target from target list

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_CSV_PATH, --input_csv_path INPUT_CSV_PATH
                        path to the target list.
                        Default is ../../data/interim/target_list.csv
  --how {random,eq_random,head,tail}
                        how to select the subset. Default is eq_random
                        random: select a random subset of the target list.
                        eq_random: select a random subset of the target list
                        with the same size of the is_similar_AF2.
                        head:select the first target_num targets of the target list.
                        tail: select the last target_num targets of the target list.
  --target_num TARGET_NUM
                        number of target to select. Default is 100.
  --random_state RANDOM_STATE. Default is 0.
  -o OUTPUT_CSV_PATH, --output_csv_path OUTPUT_CSV_PATH
                        output path of subset of target list.
```
Subset list will be created at `../../data/interim/target_subset_how_{how}_num_{target_num}_seed_{seed}.csv`
