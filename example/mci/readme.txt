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
# cal piscore
mciscore.py -r 2OH4A_receptor.pdb -l CHEMBL242865_docking.pdb -p feature_receptor.txt -q feature_ligand.txt -i interaction.txt --print_option 1
# draw interaction
mci_draw.py -r 2OH4A_receptor.pdb -l CHEMBL242865_docking.pdb -p feature_receptor.txt -q feature_ligand.txt -i interaction.txt -o fig.pse -f fig.png
