#!/usr/bin/env python
import sys
import numpy as np
import argparse


def read_pdb_protein(file_name):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    protein = dict()
    for line in lines:
        if line[0:6] != "ATOM  ":
            continue
        res_name = line[17:20].strip()
        res_num = int(line[22:26].strip())
        res_num2 = line[22:27].strip()
        chain_id = line[21]
        atom_name = line[12:16]
        atom_num = int(line[6:11])
        coor = np.array((float(line[30:38]), float(
            line[38:46]), float(line[46:54])))
        if chain_id not in protein:
            protein[chain_id] = dict()
        if res_num2 not in protein[chain_id]:
            protein[chain_id][res_num2] = dict()
            protein[chain_id][res_num2]["res_name"] = res_name
            protein[chain_id][res_num2]["atom"] = dict()

        protein[chain_id][res_num2]["atom"][atom_num] = [atom_name, coor]

    return protein


def read_pdb_ligand(file_name, ligand_name, hydrogen=False):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    ligand = dict()
    for line in lines:
        if line[0:6] != "HETATM":
            continue
        lig_name = line[17:20].strip()
        atom_name = line[12:16]
        atom_num = int(line[6:11])
        if lig_name != ligand_name and ligand_name != "":
            continue
        if line[14] == "H" and not hydrogen:
            continue
        coor = np.array((float(line[30:38]), float(
            line[38:46]), float(line[46:54])))
        ligand[atom_num] = [atom_name, coor]

    return ligand


def find_contact_residue(protein, ligand, dist_cut):

    contact_residue = list()
    for chain_name in protein:
        chain = protein[chain_name]
        for res_num in chain:
            residue = chain[res_num]
            res_name = residue['res_name']
            res_atoms = residue['atom']
            check_res = False
            for atom_num in res_atoms:
                protein_atom = res_atoms[atom_num]
                p_atom_name = protein_atom[0]
                p_atom_coor = protein_atom[1]
#                print(chain_name, res_num, res_name, atom_num, atom)
                for l_atom_num in ligand:
                    ligand_atom = ligand[l_atom_num]
                    l_atom_name = ligand_atom[0]
                    l_atom_coor = ligand_atom[1]

                    dist = np.linalg.norm(p_atom_coor-l_atom_coor)
                    if dist < dist_cut:
                        check_res = True
                        contact_residue += [[chain_name, res_num, res_name]]
                        break
                if check_res:
                    break
    return contact_residue


def main():

    parser = argparse.ArgumentParser(description='find binding residue')
    parser.add_argument('-r', '--receptor_file', type=str, required=True,
                        help='protein PDB file, ex: -r 2rgpA_receptor.pdb')
    parser.add_argument('-l', '--ligand_file', type=str, required=True,
                        help='ligand PDB file, ex: -l 2rgpA_HYZ.pdb')
    parser.add_argument('-n', '--ligand_name', type=str, default="",
                        required=False, help='ligand PDB file, ex: -n HYZ')
    parser.add_argument('-d', '--dist_cut', type=float, default=4.5,
                        required=False, help='dist_cutoff, default: 4.5')
    parser.add_argument('-o', '--output_option', type=int, default=3,
                        required=False, help='1:table, 2:csv for smina' +
                        ', 3: csv for prepare_flexreceptor4.py, default: 2')

    args = parser.parse_args()

    receptor_file = args.receptor_file
    ligand_file = args.ligand_file
    ligand_name = args.ligand_name.strip()
    dist_cut = args.dist_cut
#    hydrogen = args.hydrogen
    output_option = args.output_option
    if not isinstance(output_option, int):
        print("bad output option", output_option)
        sys.exit()
    elif output_option > 3:
        print("bad output option", output_option)
        sys.exit()
    receptor = read_pdb_protein(receptor_file)
    ligand = read_pdb_ligand(ligand_file, ligand_name)

#    print(protein)
#    print(ligand)

    contact_residue = find_contact_residue(receptor, ligand, dist_cut)

    receptor_name = receptor_file[:-4]

    receptor_name = receptor_file.split('/')[-1].split('.')[-2]
    if output_option == 1:
        for cr in contact_residue:
            chain, res_num, res_name = cr
            line_out = '%s:%s:%s' % (chain.strip(), res_num, res_name)
            print(line_out)
    elif output_option == 2:
        line_out = ""
        for cr in contact_residue:
            chain, res_num, res_name = cr
            line_out += '%s:%s,' % (chain.strip(), res_num)
        line_out = line_out.strip(',')
        print(line_out)
    elif output_option == 3:
        line_out = ""
        for cr in contact_residue:
            chain, res_num, res_name = cr
            line_out += '%s:%s:%s,' % (receptor_name,
                                       chain.strip(), res_name+res_num)
        line_out = line_out.strip(',')
        print(line_out)


if __name__ == "__main__":
    main()
