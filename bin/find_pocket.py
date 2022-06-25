#!/usr/bin/env python
import sys
import numpy as np
import argparse

Res31 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
         'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
         'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
         'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
         'ASX': 'N', 'GLX': 'Q', 'UNK': 'X', 'INI': 'K', 'AAR': 'R',
         'ACE': 'X', 'ACY': 'G', 'AEI': 'T', 'AGM': 'R', 'ASQ': 'D',
         'AYA': 'A', 'BHD': 'D', 'CAS': 'C', 'CAY': 'C', 'CEA': 'C',
         'CGU': 'E', 'CME': 'C', 'CMT': 'C', 'CSB': 'C', 'CSD': 'C',
         'CSE': 'C', 'CSO': 'C', 'CSP': 'C', 'CSS': 'C', 'CSW': 'C',
         'CSX': 'C', 'CXM': 'M', 'CYG': 'C', 'CYM': 'C', 'DOH': 'D',
         'EHP': 'F', 'FME': 'M', 'FTR': 'W', 'GL3': 'G', 'H2P': 'H',
         'HIC': 'H', 'HIP': 'H', 'HTR': 'W', 'HYP': 'P', 'KCX': 'K',
         'LLP': 'K', 'LLY': 'K', 'LYZ': 'K', 'M3L': 'K', 'MEN': 'N',
         'MGN': 'Q', 'MHO': 'M', 'MHS': 'H', 'MIS': 'S', 'MLY': 'K',
         'MLZ': 'K', 'MSE': 'M', 'NEP': 'H', 'NPH': 'C', 'OCS': 'C',
         'OCY': 'C', 'OMT': 'M', 'OPR': 'R', 'PAQ': 'Y', 'PCA': 'Q',
         'PHD': 'D', 'PRS': 'P', 'PTH': 'Y', 'PYX': 'C', 'SEP': 'S',
         'SMC': 'C', 'SME': 'M', 'SNC': 'C', 'SNN': 'D', 'SVA': 'S',
         'TPO': 'T', 'TPQ': 'Y', 'TRF': 'W', 'TRN': 'W', 'TRO': 'W',
         'TYI': 'Y', 'TYN': 'Y', 'TYQ': 'Y', 'TYS': 'Y', 'TYY': 'Y',
         'YOF': 'Y', 'FOR': 'X', '---': '-', 'PTR': 'Y', 'LCX': 'K',
         'SEC': 'D', 'MCL': 'K', 'LDH': 'K'}

typecode = "ARNDCQEGHILKMFPSTWYVX"
Ntype = 21
code_dict = dict()


def read_pdb(pdb_file, hydrogen=False):

    main_chain_atoms = ['N', 'CA', 'C', 'O']
    fp_file = open(pdb_file)
    lines = fp_file.readlines()
    fp_file.close()

    chain_dict = dict()
    for line in lines:
        if line[:6] == 'ATOM  ':
            chain_id = line[21]
            res_idx = line[22:27]  # not integer
            atom_idx = int(line[6:11].strip())
            atomname = line[12:16]
            atomname2 = atomname.strip()
            res_name = line[17:20]
            altloc = line[16]
            coor = np.array([float(line[30:38]), float(line[38:46]),
                             float(line[46:54])], dtype=np.float32)
            if altloc != ' ' and altloc != 'A':
                continue
            if atomname2[0] == 'H' and not hydrogen:
                continue
            if chain_id not in chain_dict:
                chain_dict[chain_id] = dict()
            if res_idx not in chain_dict[chain_id]:
                chain_dict[chain_id][res_idx] = {
                    'res_name': res_name, 'atom': dict()}
            chain_dict[chain_id][res_idx]['atom'][atom_idx] = {
                'atomname': atomname, 'coor': coor,
                'line': line}

    chain_keys = sorted(chain_dict.keys())
    for chain_id in chain_keys:
        chain = chain_dict[chain_id]
        res_idx_keys = sorted(chain.keys())
        for res_idx in res_idx_keys:
            residue = chain[res_idx]
            res_name = residue['res_name']
            atom_dict = residue['atom']
            atom_keys = sorted(atom_dict.keys())
            sg_coor = []
            for atom_idx in atom_keys:
                atom = atom_dict[atom_idx]
                atomname = atom['atomname']
                coor = atom['coor']
                if atomname == ' CA ':
                    chain_dict[chain_id][res_idx]['CA'] = coor
                    if res_name != 'GLY':
                        chain_dict[chain_id][res_idx]['CB'] = coor
                if atomname == ' CB ':
                    chain_dict[chain_id][res_idx]['CB'] = coor
                atomname2 = atomname.strip()
                if res_name != 'GLY':
                    if (atomname2 not in main_chain_atoms) and (atomname2[0] != 'H'):
                        sg_coor += [coor]
            if res_name != 'GLY':
                sg_coor = np.concatenate([sg_coor]).mean(axis=0)
            else:
                sg_coor = chain_dict[chain_id][res_idx]['CA']
            chain_dict[chain_id][res_idx]['SG'] = sg_coor
    return chain_dict


def read_pdb_ligand(file_name, ligand_name, hydrogen=False):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    coor_list = []
    for line in lines:
        if line[0:6] != "HETATM" and line[0:6] != "ATOM  ":
            continue
        lig_name = line[17:20].strip()
        if lig_name != ligand_name and ligand_name != '':
            continue
        if line[13] == 'H' and not hydrogen:
            continue
        coor = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        coor_list += [coor]
    coor = np.array(coor_list)

    return coor


def read_mol2_ligand(file_name, ligand_name, hydrogen=False):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    coor_list = []
    line = lines[3]

    state = ''
    line_count = 0
    for line in lines:
        if line.startswith('@<TRIPOS>'):
            state = line[9:].strip()
            line_count = 0
            continue
        if state == "MOLECULE":
            line_count += 1
            if line_count == 2:
                num_atoms = int(line[0:5])
                num_bonds = int(line[5:11])
            continue
        if state == "ATOM":
            lig_name = line[59:62].strip()
            if lig_name != ligand_name and ligand_name != '':
                continue
            if line[47] == 'H' and not hydrogen:
                continue
            coor = [float(line[16:26]), float(line[26:36]), float(line[36:46])]
            coor_list += [coor]
    coor = np.array(coor_list)

    return coor


def read_mol_ligand(file_name, ligand_name, hydrogen=False):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    coor_list = []
    line = lines[3]
    num_atoms = int(line[0:3])
    num_bonds = int(line[3:6])
    mol_version = line[34:39]

    for line in lines[4:4+num_atoms]:
        if line[31] == 'H' and not hydrogen:
            continue
        coor = [float(line[0:10]), float(line[10:20]), float(line[20:30])]
        coor_list += [coor]
    coor = np.array(coor_list)

    return coor


def find_pocket(chain_dict, ligand_coor, dist_cutoff=5):

    cmin = ligand_coor.min(axis=0) - dist_cutoff
    cmax = ligand_coor.max(axis=0) + dist_cutoff
    num_ligand_atoms = ligand_coor.shape[0]
    chain_pocket = dict()
    pocket_res_list = list()
    chain_keys = sorted(chain_dict.keys())
    for chain_id in chain_keys:
        chain = chain_dict[chain_id]
        res_idx_keys = sorted(chain.keys())
        for res_idx in res_idx_keys:
            pocket_residue = False
            residue = chain[res_idx]
            res_name = residue['res_name']
            atom_dict = residue['atom']
            atom_keys = sorted(atom_dict.keys())
            for atom_idx in atom_keys:
                atom = atom_dict[atom_idx]
                atomname = atom['atomname']
                atomname2 = atomname.strip()
                if atomname2[0] == 'H':
                    continue
                coor = atom['coor']
                if (coor < cmin).any() or (coor > cmax).any():
                    continue
                dist_min = 100
                for j in range(0, num_ligand_atoms):
                    dist = np.linalg.norm(coor-ligand_coor[j])
                    if dist < dist_min:
                        dist_min = dist
                        lc = ligand_coor[j]
                    if dist < dist_cutoff:
                        pocket_residue = True
                        break
                if pocket_residue:
                    break
            if pocket_residue:
                if chain_id not in chain_pocket:
                    chain_pocket[chain_id] = dict()
                chain_pocket[chain_id][res_idx] = residue
                pocket_res = '%s:%s' % (chain_id,
                                        Res31[res_name] + res_idx.strip())
                pocket_res_list += [pocket_res]
    return chain_pocket, pocket_res_list


def write_pocket(chain_dict, pocket_file):

    fp_out = open(pocket_file, 'w')
    chain_keys = sorted(chain_dict.keys())
    for chain_id in chain_keys:
        chain = chain_dict[chain_id]
        res_idx_keys = sorted(chain.keys())
        pocket_residue = False
        for res_idx in res_idx_keys:
            residue = chain[res_idx]
            res_name = residue['res_name']
            atom_dict = residue['atom']
            atom_keys = sorted(atom_dict.keys())
            for atom_idx in atom_keys:
                atom = atom_dict[atom_idx]
                line_out = atom['line']
                fp_out.write(line_out)
    fp_out.close()
    return


def write_binding_residue_file(pocket_res_list, binding_file):
    fp_out = open(binding_file, 'w')
    binding_chain = dict()
    for pocket_res in pocket_res_list:
        line_out = pocket_res + '\n'
        fp_out.write(line_out)
    fp_out.close()
    return


def main():
    parser = argparse.ArgumentParser(description='find_pocket')
    parser.add_argument('-p', '--protein_file', type=str,
                        required=True, help='protein file, ' +
                        'ex: -p protein.pdb')
    parser.add_argument('-o', '--output_file', type=str,
                        required=True, help='output protein file, ' +
                        'ex: -o pocket.pdb')
    parser.add_argument('-b', '--residue_file', type=str,
                        required=False, default='binding_residue.txt',
                        help='output file for binding_residue' +
                        ', ex: -o binding_residue.txt')

    parser.add_argument('-l', '--ligand_file', type=str,
                        required=True, help='ligand file, ' +
                        'ex: -l ligand.pdb, only pdb, mol2, mol, sdf')

    parser.add_argument('-n', '--ligand_name', type=str, default='',
                        required=False, help='ligand name in PDB, ' +
                        'ex: -n HYZ, default: -n ""')

    parser.add_argument('-c', '--dist_cutoff', type=float, default=8.0,
                        required=False, help='distant cutoff, '
                        'unit:A, default: 8.0 ')

    parser.add_argument('--hydrogen', action='store_true',
                        required=False,
                        help='--hydrogen: include hydrogen atoms or' +
                        ' exclude hydrogen atoms ')

    args = parser.parse_args()
    protein_file = args.protein_file
    ligand_file = args.ligand_file
    ligand_name = args.ligand_name.strip()
    pocket_file = args.output_file
    binding_residue_file = args.residue_file

    dist_cutoff = args.dist_cutoff
    hydrogen = args.hydrogen

    chain_dict = read_pdb(protein_file, hydrogen=hydrogen)

    file_format = ligand_file.strip().split('.')[-1].lower()
    if file_format == 'pdb':
        ligand_coor = read_pdb_ligand(ligand_file, ligand_name, hydrogen)
    elif file_format == 'mol' or file_format == 'sdf':
        ligand_coor = read_mol_ligand(ligand_file, ligand_name, hydrogen)
    elif file_format == 'mol2':
        ligand_coor = read_mol2_ligand(ligand_file, ligand_name, hydrogen)
    else:
        print('file format ', file_format, 'is not supported')

    chain_pocket, pocket_res_list = find_pocket(
        chain_dict, ligand_coor, dist_cutoff)
#    print(chain_pocket)
#    print(pocket_res_list)
    write_pocket(chain_pocket, pocket_file)
    write_binding_residue_file(pocket_res_list, binding_residue_file)

    binding_chain = dict()
    for pocket_res in pocket_res_list:
        chain_id = pocket_res[0]
        if chain_id not in binding_chain:
            binding_chain[chain_id] = 0
        binding_chain[chain_id] += 1
#    print(binding_chain)
    line_out = ''
    for chain_id in binding_chain:
        cnum = binding_chain[chain_id]
        line_out += '%s:%d,' % (chain_id, cnum)
    line_out = line_out.strip(',')
    print(line_out)


if __name__ == "__main__":
    main()
