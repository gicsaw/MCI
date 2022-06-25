#!/usr/bin/env python
# import pymol
from pymol import cmd
# from pymol import stored
import argparse


def main():

    title_line = 'show_binding_residue'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-r', '--receptor', type=str, required=True,
                        help='receptor pdb file')
    parser.add_argument('-l', '--ligand', type=str, required=True,
                        help='ligand pdb file')
    parser.add_argument('-b', '--binding_residue', type=str, required=True,
                        help='binding residue list file')
    parser.add_argument('--show_cartoon', action='store_true', required=False,
                        help='show cartoon for receptor.')
    parser.add_argument('--show_all_hydrogen', action='store_true',
                        required=False, help='show all hydrogen.')

    parser.add_argument('-o', '--out_pymol', type=str, required=False,
                        default='fig.pse', help='output pymol')
    parser.add_argument('-f', '--out_fig', type=str, required=False,
                        default='fig.png', help='output figure')

    args = parser.parse_args()

    receptor_file = args.receptor
    ligand_file = args.ligand
    binding_residue_file = args.binding_residue
    pymol_file = args.out_pymol
    fig_file = args.out_fig
    show_cartoon = args.show_cartoon
    show_all_hydrogen = args.show_all_hydrogen

    receptor_obj = receptor_file.split('/')[-1].split('.')[-2]
    ligand_obj = ligand_file.split('/')[-1].split('.')[-2]

#    receptor_obj = 'receptor'
#    ligand_obj = 'ligand'

    cmd.load(receptor_file, receptor_obj)
    cmd.load(ligand_file, ligand_obj)
    cmd.util.cbaw(ligand_obj)

    object_list = cmd.get_object_list(selection='(all)')
    if receptor_obj not in object_list:
        print('Error:', receptor_obj, 'is not in', object_list)
    if ligand_obj not in object_list:
        print('Error:', ligand_obj, 'is not in', object_list)

    fp = open(binding_residue_file)
    lines = fp.readlines()
    fp.close()
    for line in lines:
        lis = line.strip().split(':')
        chain, res_num, res_name = lis
        ss = '/%s//%s/%s`%s' % (receptor_obj, chain, res_name, res_num)
        cmd.show('sticks', ss)
        cmd.label(ss+'/CA', 'resn+resi')

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
