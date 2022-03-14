# MQAB4AF2

Model Quality Assessment benchmarking for AlphaFold2 structures

## Dataset generation
### Requirement

```txt
python=3.9.2
requests=2.26.0
numpy=1.21.3
pandas=1.3.4
tqdm=4.62.3
prody=2.0
joblib=1.1.0
dssp=3.0.0
pytest=6.2.5
```

### Procedure
1. [Target selection](./src/data/README.md)
2. [Run AlphaFold2](./src/alphafold/README.md)
3. [Optional] [Run MQA methods](./src/mqa/README.md)

## Download dataset
You can download AlphaFold2 structure data for 500 protein sequences from [http://www.cb.cs.titech.ac.jp/af2/af2_500_targets.tar.gz](http://www.cb.cs.titech.ac.jp/af2/af2_500_targets.tar.gz) (5.8GB).
http://www.cb.cs.titech.ac.jp/af2/af2_500_targets.tar.gz

The dataset includes protein target sequences, native and predicted structures of the targets, and labels.
For more information, see [here](./data/out/dataset/compress/README.md).

## Reference
Yuma Takei and Takashi Ishida, in preparation, 2022.


## Acknowledgement
[AlphaFold](https://github.com/deepmind/alphafold) v2.0.1, [ColabFold](https://github.com/sokrypton/ColabFold), and [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold) v1.0.0 were used to predict structures.

- Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, Tunyasuvunakool K, Bates R, Žídek A, Potapenko A and Bridgland A. "Highly accurate protein structure prediction with AlphaFold."<br />
  Nature. 2021 Aug;596(7873):583-9. doi: [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)

- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. "ColabFold - Making protein folding accessible to all." <br />
  bioRxiv (2021) doi: [10.1101/2021.08.15.456425](https://www.biorxiv.org/content/10.1101/2021.08.15.456425v2)

- Moriwaki Y. "LocalColabFold" Available from: https://github.com/YoshitakaMo/localcolabfold.
