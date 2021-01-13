### Installation
1. Download conda
+ https://conda.io/projects/conda/en/latest/user-guide/install/linux.html
2. Install conda
```
bash Miniconda2-4.5.11-Linux-x86_64.sh
```
3. Activate .bashrc
```
source ~/.bashrc
```
### Create, activate and quit an environment
1. Create an environment
```
conda env create -n rna --file RNA_environment.yaml
# rna is the name for this specific environment
```
2. Activiate the enviroment rna
```
conda activate rna
```
3. Quit the present enviroment
```
conda deactivate
```
