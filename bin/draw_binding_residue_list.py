#!/usr/bin/env python
# import pymol
from pymol import cmd
# from pymol import stored
import argparse


def draw(data_dir, pdb_id, ligand_list):

    receptor_obj = '%s_receptor' % (pdb_id)
    receptor_file = '%s/%s_receptor.pdb' % (data_dir, pdb_id)
    cmd.load(receptor_file, receptor_obj)

    for ligand_id in ligand_list:
        ligand_obj = '%s_%s' % (pdb_id, ligand_id)
        ligand_file = '%s/%s_%s.pdb' % (data_dir, pdb_id, ligand_id)
        cmd.load(ligand_file, ligand_obj)
        cmd.util.cbaw(ligand_obj)

    binding_residue_file = '%s/%s_binding_residue.txt' % (data_dir, pdb_id)

    object_list = cmd.get_object_list(selection='(all)')
    if receptor_obj not in object_list:
        print('Error:', receptor_obj, 'is not in', object_list)

    fp = open(binding_residue_file)
    lines = fp.readlines()
    fp.close()
    for line in lines:
        lis = line.strip().split(':')
        chain, res_num, res_name = lis
        ss = '/%s//%s/%s`%s' % (receptor_obj, chain, res_name, res_num)
        cmd.show('sticks', ss)
        cmd.label(ss+'/CA', 'resn+resi')

    show_all_hydrogen = True
    show_cartoon = True
    cmd.hide('cartoon', receptor_obj)
    if not show_all_hydrogen:
        cmd.hide('(h. and (e. c extend 1))')
    cmd.orient('vis')
    if show_cartoon:
        cmd.show('cartoon', receptor_obj)


def main():

    parser = argparse.ArgumentParser(description='draw_binding_residue_list')
    parser.add_argument('-l', '--list_file', type=str, required=True,
                        help='PDB LIGAND list file, ex: -l list.txt')
    parser.add_argument('-d', '--data_dir', type=str, required=True,
                        help='data directory, ex: pdb')
    parser.add_argument('-o', '--out_pymol', type=str, required=False,
                        default=None, help='output pymol file, ex: fig.pse')
    parser.add_argument('-p', '--out_png', type=str, required=False,
                        default=None, help='output png file, ex: fig.png')

    args = parser.parse_args()
    list_file = args.list_file
    data_dir = args.data_dir
    out_pymol = args.out_pymol
    out_png = args.out_png

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    for line in lines:
        if line.startswith('#'):
            continue
        lis = line.strip().split()
        pdb_id = lis[0]
        ligand_list = lis[1:]
        draw(data_dir, pdb_id, ligand_list)

    cmd.bg_color('white')
    # cmd.space('cmyk')
    # cmd.space('pymol')
    cmd.space('rgb')
    cmd.util.performance(100)
    cmd.rebuild()
    cmd.save(out_pymol)
    cmd.png(out_png, 0, 0, -1, ray=0)


if __name__ == '__main__':
    main()
