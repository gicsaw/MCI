#!/usr/bin/env python
import sys
from mci import mcifinder


def main():

    import argparse
    title_line = 'Fixer for ligand pdb which is converted from pdbqt'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-r', '--mciscore_receptor', required=False,
                        default=None,
                        help='input receptor pdb file')
    parser.add_argument('-p', '--pf_receptor', required=False,
                        default=None,
                        help='output for pharmacophoric feature of receptor')
    parser.add_argument('-t', '--mcinfo_ligand', required=False,
                        default=None, help='mcinfo_ligand')
    parser.add_argument('-u', '--pf_receptor_info', required=False,
                        default=None,
                        help='output for template feature of receptor')
    parser.add_argument('-v', '--dock_config', type=str, required=False,
                        default=None, help='dock_config_file')
    parser.add_argument('-m', '--mci_cutoff', type=float, required=False,
                        default=6.5,
                        help='pharmacophoric interaction cutoff distance')
    parser.add_argument('-l', '--ligand_file', required=True,
                        help='input ligand pdb file')
    parser.add_argument('-q', '--pf_ligand_file', required=False,
                        default=None,
                        help='write pharmacophoric feature of ligand')
    parser.add_argument('--include_hydrophobic', action='store_true',
                        required=False,
                        help='include hydrophobic feature for template')
    parser.add_argument('-i', '--interaction_file', required=False,
                        default=None, help='write interaction')
    parser.add_argument('--print_option', type=int, required=False, default=0,
                        help='print_option 0:None, 1: one_line, 2: multi_line')

#    draw_ligand = False
    args = parser.parse_args()
    if args.mciscore_receptor is None and args.pf_receptor is None:
        parser.print_usage()
        print('error mciscore_receptor and pf_receptor are None')
        sys.exit()

    pharma = mcifinder.set_mcifinder(args)
    use_mcinfo = pharma.use_mcinfo

    ligand_file = args.ligand_file
    pf_ligand_file = args.pf_ligand_file
    interaction_file = args.interaction_file
    print_option = args.print_option

    result = pharma.cal_mciscore(ligand_file, pf_ligand_file, interaction_file)
    total_count_dict, total_count_info_dict = result

#    if draw_ligand:
#        size = (200, 200)
#        fig_dir = 'fig'
#        output_name = fig_dir + '/%s.png' % (ligand_file[:-3])
#        m_rdkit_ligand = Chem.MolFromPDBFile(ligand_file, removeHs=True)
#        pharma.draw_ligand(m_rdkit_ligand, PF_dict_ligand, output_name, size)

    if print_option > 0:
        lines_count = pharma.print_interaction(total_count_dict,
                                               total_count_info_dict,
                                               use_mcinfo,
                                               print_option=print_option)
        print(lines_count.rstrip())


if __name__ == '__main__':
    main()
