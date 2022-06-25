# Medicinal Chemistry Interaction Finder

## environment
see environment.yml
```
openbabel >=3.1.0
smina
rdkit >=2020
pandas
pdbfixer
pymol-open-source >= 2.4
``` 

If you don't have a protein-ligand binding structure, you should use a docking program.
```
Qvina: A fork that improves the speed of Autodock vina.
https://github.com/QVina/qvina
(Conda installation is also supported, but this version is not very fast.
It is recommended to download and install directly from github.)
or
smina: A fork with improved convenience of Autodock vina.
conda install -c conda-forge smina 

autodocktools-prepare : It is used when converting a pdb file to a pdbqt file.
You can use obabel, but this sometimes misidentifies rotatable bonds

There are the following python3 conversion.
https://github.com/jaimergp/autodocktools-prepare-py3k 
or
https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3

```

## Install
``` 
git clone https://github.com/gicsaw/MCIFinder.git
```
Add the following to your .bashrc file
```
MCIFiner_PATH="$YOUR PATH/"
export PYTHONPATH="$MCIFiner_PATH:$PYTHONPATH"
export PATH="$MCIFiner_PATH/bin:$PATH"
source .bashrc
```


## example
### pdb preparation
There is a script that prepares the pdb file in the MCIFinder/script path.
Examples are prepared in MCIFinder/example/pdb_preparation.

```
mkdir test
cd test
cp $MCIFiner_PATH/example/pdb_preparation/vgfr2/list.txt .
cp $MCIFiner_PATH/script/* .
./auto_run.sh
```
It is recommended not to run auto_run.sh immediately and check the contents.

### MCIFinder
MCIscore and MCIdraw: single model

```
mkdir test
cd test
cp $MCIFiner_PATH//example/pifinder/* .

# docking preperation
smi="N(c1ccnc(Nc2ccc(cc2)[N+](=O)[O-])n1)c1c(cccc1)C(=O)O"
smi_0=`obabel -:"$smi" --neutralize -osmi`
smi_p=`obabel -:"$smi_0" -p7.4 -osmi`
echo $smi_p
# gen 3d conformation
obabel -:"$smi_p" --gen3d -O CHEMBL242865_ini.pdb
# pdb to pdbqt
prepare_ligand4.py -l CHEMBL242865_ini.pdb -o CHEMBL242865_ini.pdbqt
# docking
smina --config config.txt --ligand CHEMBL242865_ini.pdbqt --out CHEMBL242865_docking.pdbqt

# pdbqt to pdb
qt2pdb.py -i CHEMBL242865_docking.pdbqt -o CHEMBL242865_docking.pdb
# cal mciscore
mciscore.py -r 2OH4A_receptor.pdb -l CHEMBL242865_docking.pdb -p feature_receptor.txt -q feature_ligand.txt -i interaction.txt --print_option 1
# draw interaction
pi_draw.py -r 2OH4A_receptor.pdb -l CHEMBL242865_docking.pdb -p feature_receptor.txt -q feature_ligand.txt -i interaction.txt -o fig.pse -f fig.png
```

### MCIscore with PVSdock

PVSdock: 
https://github.com/gicsaw/PVSdock
install PVSdock 

```
# feature preparation
# extract pharmacophoric feature from "pdb/receptor.pdb" and write it to "pdb/pf_receptor.txt".
# Select a feature that interacts with the ligand "pdb/crystal_ligand.pdb" and write it to "pdb/pf_receptor_template.txt".
# Alse you can manually modify the "pdb/pf_receptor_template.txt" file. 
find_pf_receptor.py -v config.txt -r pdb/receptor.pdb -p pdb/pf_receptor.txt -t pdb/crystal_ligand.pdb -u pdb/pf_receptor_template.txt

# There are two options for using mciscore and PVSdock.
# One is to first find the docking pose with PVSdock and then calculate the MCIscore after PVSdock is over.
# The other is to compute the MCIscore together while the PVSdock is running. 

# If the docking structure using PVSdock is prepared in a directory such as "docking/MOL_ID[0:7]/dock_MOL_ID.pdb" 
# and the pandas dataframe as a result of docking is saved in "docking.pkl", MCIscore can be obtained using the following code. 
python mciscore_vhts.py -v config.txt -p pdb/pf_receptor.txt -u pdb/pf_receptor_template.txt -i workflow/master/docking.pkl -d docking --num_procs 2
# Example directory structure: 
# docking/CHEMBL1/dock_CHEMBL242865.pdb


# If you want to run PVSdock and MCIscore in conjunction, refer to the https://github.com/gicsaw/PVSdock. 

```
