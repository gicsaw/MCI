#!/usr/bin/env python
from mci import mcifinder
import argparse
import sys


def main():

    title_line = 'picount for vhts'
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
    parser.add_argument('--include_hydrophobic', action='store_true',
                        required=False,
                        help='include hydrophobic feature for template')

    args = parser.parse_args()
    if args.mciscore_receptor is None and args.pf_receptor is None:
        parser.print_usage()
        print('error mciscore_receptor and pf_receptor are None')
        sys.exit()

    pharma = mcifinder.set_mcifinder(args)


if __name__ == "__main__":
    main()
