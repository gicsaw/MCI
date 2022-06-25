#!/usr/bin/env python
import pymol
from pymol import cmd
# from pymol import stored
import argparse


def read_feature_file(feature_file, model_id=None):
    fp = open(feature_file)
    lines = fp.readlines()
    fp.close()
    feature_dict = dict()
    model_id0 = None
    for line in lines:
        if line.startswith('feature_type'):
            #            title = line.strip().split(':')
            continue
        if line.startswith('MODEL'):
            model_id0 = int(line.strip().split()[1])
            continue
        if model_id is None:
            if model_id0 is not None:
                if model_id0 > 1:
                    continue
        else:
            if model_id != model_id0:
                continue
        lis = line.strip().split(':')
        feature_type = lis[0].strip()
        atom_idx_list = tuple([int(x) for x in
                               lis[1].lstrip('(').rstrip(')').split(',')])
        center_coor = tuple([float(x) for x in
                             lis[3].lstrip('(').rstrip(')').split(',')])
        feature_dict[atom_idx_list] = (feature_type, center_coor)


    return feature_dict


def read_pdb_file(pdb_file, model_id):
    fp = open(pdb_file)
    lines = fp.readlines()
    fp.close()
    model_id0 = None
    atom_dict = dict()
    k = 0
    for line in lines:
        if line[0:6] == 'MODEL ':
            model_id0 = int(line.strip().split()[1])
        if model_id is None:
            if model_id0 is not None:
                if model_id0 > 1:
                    continue
        else:
            if model_id != model_id0:
                continue
        if line[0:6] == 'ENDMDL':
            if model_id is not None and model_id0 is not None:
                if model_id == model_id0:
                    break
        if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
            residue_name = line[17:20].strip()
            residue_num = int(line[22:26])
            chain_id = line[21].strip()
            atom_number = int(line[6:11])
            atom_name = line[12:16].strip()
            k += 1
            atom_dict[k] = (atom_number, atom_name,
                            residue_num, residue_name, chain_id)
    return atom_dict


def draw_mci_receptor(receptor_atom_dict, receptor_feature_dict,
             receptor_obj):

    pymol.one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E',
                        'GLN': 'Q', 'ASP': 'D', 'ASN': 'N', 'HIS': 'H',
                        'TRP': 'W', 'PHE': 'F', 'TYR': 'Y', 'ARG': 'R',
                        'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M',
                        'ALA': 'A', 'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}

    feature_name_dict = {'HBD': 'HBD', 'HBA': 'HBA',
                         'Cation': '+', 'Anion': '-',
                         'Aromatic': 'Aro', 'Metal': 'M', 'MBA': 'MBA'}

    feature_color_dict = {'HBD': 'cyan', 'HBA': 'magenta', '+': 'yellow',
                          '-': 'orange', 'Aro': 'wheat',
                          'M': 'palegreen', 'MBA': 'lightblue'}

    receptor_pseudo_idx = 0
    receptor_feature_idx = dict()

    cmd.set('sphere_scale', 0.5)
    cmd.set('label_size', 15.0)
    cmd.set('label_font_id', 7)

    receptor_obj2 = 'receptor_feature'
    for atom_idx_list in receptor_feature_dict:

        (feature_type, pseudo_coor) = receptor_feature_dict[atom_idx_list]

        if feature_type == 'Hydrophobic':
            continue
        fea_rec = feature_name_dict[feature_type]

        receptor_idx_tuple = atom_idx_list
        atom_receptor = receptor_atom_dict[receptor_idx_tuple[0]]
        (r_atom_number, r_atom_name, r_residue_num,
         r_residue_name, r_chain_id) = atom_receptor

        fc2_receptor = '/%s//%s/%s`%d' % (receptor_obj, r_chain_id,
                                          r_residue_name, r_residue_num)
#        fc2_receptor2 = '/%s//%s/%s`%d' % (receptor_obj2, r_chain_id,
#                                          r_residue_name, r_residue_num)

        cmd.show('sticks', fc2_receptor)
#        cmd.select(fc2_receptor)
#        cmd.copy_to(receptor_obj2, 'sele')

#        cmd.label(fc2_receptor+'/CA', 'chain+":"+one_letter[resn]+resi')
#        cmd.label(fc2_receptor+'/CA', 'one_letter[resn]+resi')
        cmd.label(fc2_receptor+'/CA', 'resn+resi')

        fc1_receptor = '/%s//%s/%s`%d/%s' % (receptor_obj2, r_chain_id,
                                             r_residue_name, r_residue_num,
                                             r_atom_name)

        num_receptor_fea_atom = len(receptor_idx_tuple)
        if num_receptor_fea_atom > 0:
            if fea_rec not in receptor_feature_idx:
                receptor_feature_idx[fea_rec] = 0
            receptor_feature_idx[fea_rec] += 1
            receptor_pseudo_idx += 1
#            f_name_r = 'aromatic%d' %(receptor_pseudo_idx)
#            f_name_r = receptor_obj
            cmd.pseudoatom(receptor_obj2, resi=1,
                           chain='P', color=feature_color_dict[fea_rec],
                           label=fea_rec,
                           pos=pseudo_coor)
            f_name_r = 'PS%d' % (receptor_pseudo_idx)
            fc1_receptor = '/%s/PSDO/P/PSD`1/%s' % (
                receptor_obj2, f_name_r)
            cmd.show('spheres', fc1_receptor)
#            fea_rec2 = '"' + fea_rec + '"'
#           cmd.label(fc1_receptor, fea_rec2)
            cmd.set('sphere_scale', 0.5, fc1_receptor)


#    fc1_ligand = '/%s/PSDO/%s/PSD`%d/PS1' % (
#        ligand_obj2, l_chain_id, l_residue_num)
#    cmd.remove(fc1_ligand)


def main():

    title_line = 'mci_draw'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-r', '--receptor', type=str, required=True,
                        #                        default='receptor.pdb',
                        help='receptor pdb file')

    parser.add_argument('-p', '--feature_receptor', type=str, required=False,
                        default='feature_receptor.txt',
                        help='receptor feature file')

    parser.add_argument('--show_cartoon', action='store_true', required=False,
                        help='show cartoon for receptor.')

    parser.add_argument('--show_all_hydrogen', action='store_true',
                        required=False, help='show all hydrogen.')

    parser.add_argument('-o', '--out_pymol', type=str, required=False,
                        default='fig.pse', help='output pymol')
    parser.add_argument('-f', '--out_fig', type=str, required=False,
                        default='fig.png', help='output figure')
    parser.add_argument('-m', '--model_id', type=int, required=False,
                        default=None, help='docking model id')

    args = parser.parse_args()

    receptor_file = args.receptor
    receptor_feature_file = args.feature_receptor
    pymol_file = args.out_pymol
    fig_file = args.out_fig
    show_cartoon = args.show_cartoon
    show_all_hydrogen = args.show_all_hydrogen

    receptor_atom_dict = read_pdb_file(receptor_file, model_id=None)
    receptor_feature_dict = read_feature_file(
        receptor_feature_file, model_id=None)

    receptor_obj = 'receptor'
    cmd.load(receptor_file, receptor_obj)

    object_list = cmd.get_object_list(selection='(all)')
    if receptor_obj not in object_list:
        print('Error:', receptor_obj, 'is not in', object_list)

    draw_mci_receptor(receptor_atom_dict, receptor_feature_dict,
             receptor_obj)

    cmd.hide('cartoon', receptor_obj)

    if not show_all_hydrogen:
        cmd.hide('(h. and (e. c extend 1))')
    cmd.orient('vis')
    if show_cartoon:
        cmd.show('cartoon', receptor_obj)

    cmd.bg_color('white')
    # cmd.space('cmyk')
    # cmd.space('pymol')
    cmd.space('rgb')
    cmd.util.performance(100)
    cmd.rebuild()
    cmd.save(pymol_file)
    cmd.png(fig_file, 0, 0, -1, ray=0)


if __name__ == '__main__':
    main()
