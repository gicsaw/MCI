#!/usr/bin/env python
import os
import sys
from multiprocessing import Manager
from multiprocessing import Process
from multiprocessing import Queue
from mci import mcifinder

import numpy as np
import pandas as pd
import argparse


def creator(q, data, num_procs):
    """
        put data to queue
        input: queue
            data = [(idx1, molid1, smi1), (idx2, molid2, smi2), ...]
            num_procs (for end signal)
    """
    for d in data:
        idx = d[0]
        q.put((idx, d[1]))

    for i in range(0, num_procs):
        q.put('DONE')


def sub_worker(q, return_dict, i_dir, pharma):
    # pid = os.getpid()

    while True:
        qqq = q.get()
        if qqq == 'DONE':
            # print('proc =', os.getpid())
            break

        (idx, mol_id) = qqq
        # print screening processing in every pout step
#        pout = 1000
#        if pout != 0:
#            if idx % pout == pout-1:
#                print("processing: ", idx+1, flush=True)
        p1 = mol_id[0:7]
        i_dir1 = i_dir + '/' + p1
        input_dir = i_dir1
        ligand_file = input_dir + '/dock_' + mol_id + '.pdb'
        pf_ligand_file = input_dir + '/' + mol_id + '_pf_ligand.txt'
        interaction_file = input_dir + '/' + mol_id + '_interaction.txt'
        if not os.path.exists(ligand_file):
            result_dict = None
            return_dict[idx] = result_dict
            continue

        result_mciscore = pharma.cal_mciscore(
            ligand_file, pf_ligand_file, interaction_file)
        total_count_dict, total_count_info_dict = result_mciscore

        keys = list(total_count_dict.keys())
        mciscore = list()
        if pharma.use_mcinfo:
            mcinfo = list()
        for key in keys:
            count_dict = total_count_dict[key]
            pis = count_dict['Score']
            mciscore += [pis]
            if pharma.use_mcinfo:
                count_info_dict = total_count_info_dict[key]
                pin = count_info_dict['Score']
                mcinfo += [pin]
        mciscore = np.array(mciscore)
        if pharma.use_mcinfo:
            mcinfo = np.array(mcinfo)

        return_dict[idx] = (mciscore, mcinfo)


def main():

    title_line = 'mciscore for vhts'
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
    parser.add_argument('-v', '--dock_config', type=str, required=True,
                        help='dock_config_file')
    parser.add_argument('-m', '--mci_cutoff', type=float, required=False,
                        default=6.5,
                        help='pharmacophoric interaction cutoff distance')
    parser.add_argument('--include_hydrophobic', action='store_true',
                        required=False,
                        help='include hydrophobic feature for template')
    parser.add_argument('-i', '--docking_file', type=str, required=True,
                        help='docking result file, pickle')
    parser.add_argument('-d', '--docking_dir', type=str, required=True,
                        help='docking result dir')
    parser.add_argument('--num_procs', type=int, required=False, default=1,
                        help='num_procs')
#    parser.add_argument('-t', '--template_list_file', type=str, required=True,
#                        help='template list file')

    args = parser.parse_args()
    if args.mciscore_receptor is None and args.pf_receptor is None:
        parser.print_usage()
        print('error mciscore_receptor and pf_receptor are None')
        sys.exit()

    pharma = mcifinder.set_mcifinder(args)

    use_mcinfo = pharma.use_mcinfo
    num_procs = args.num_procs
    i_dir = args.docking_dir
    data_file = args.docking_file

    df = pd.read_pickle(data_file)
#    df = df.loc[0:2]
    num_data = df.shape[0]
    num_procs = min(num_procs, num_data)
    print('num of cpu', num_procs, flush=True)
    print('number of protein', num_data, flush=True)

    data = list()
    for idx in range(num_data):
        df_0 = df.iloc[idx]
        mol_id = df_0['MOL_ID']
        data += [[idx, mol_id]]

    q1 = Queue()
    manager = Manager()
    return_dict = manager.dict()
    proc_master = Process(target=creator,
                          args=(q1, data, num_procs))
    proc_master.start()

    # create slave process
    procs = []
    for sub_id in range(0, num_procs):
        proc = Process(target=sub_worker,
                       args=(q1, return_dict, i_dir, pharma))
        procs.append(proc)
        proc.start()

    q1.close()
    q1.join_thread()
    proc_master.join()
    for proc in procs:
        proc.join()
#    keys = sorted(return_dict.keys())

    pis_list = list()
    if use_mcinfo:
        mcinfo_list = list()
    for idx in range(num_data):
        if idx not in return_dict:
            pis_list += [[-999]]
            if use_mcinfo:
                mcinfo_list += [[0]]
            continue
        return0 = return_dict[idx]
        if return0 is None:
            pis_list += [[-999]]
            if use_mcinfo:
                mcinfo_list += [[0]]
            continue
        mciscore, mcinfo = return0
        df_0 = df.iloc[idx]
        pis_list += [mciscore]
        if use_mcinfo:
            mcinfo_list += [mcinfo]
    df['PIscore'] = pis_list
    if use_mcinfo:
        df['Pinfo'] = mcinfo_list

    out_file = 'mciscore.pkl'
    df.to_pickle(out_file)

    out_file = 'mciscore.csv'
    df.to_csv(out_file, sep=',', float_format='%.3f', index=False)


if __name__ == '__main__':
    main()
