# VoroCNN

## Install

```bash
# clone
git clone https://gitlab.inria.fr/grudinin/vorocnn.git
# Download precompliled executable file
# from https://drive.google.com/file/d/1OMR9w5z_8MdMJgHLBGXUwImtKOQ5G8d-/view
cd vorocnn
FILE_ID=1OMR9w5z_8MdMJgHLBGXUwImtKOQ5G8d-
FILE_NAME=vorocnn
curl -sc cookie "https://drive.google.com/uc?export=download&id=${FILE_ID}" > /dev/null
CODE="$(awk '/_warning_/ {print $NF}' cookie)"
curl -Lb cookie "https://drive.google.com/uc?export=download&confirm=${CODE}&id=${FILE_ID}" -o ${FILE_NAME}
rm cookie
```

## Run

```bash
vorocnn -i input_dir -o output_dir -v /path/to/voronota
```
