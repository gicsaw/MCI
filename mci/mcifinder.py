#!/usr/bin/env python
import numpy as np
from openbabel import pybel
# import time
from rdkit import Chem
from rdkit.Chem import Draw
from pbi.pdbtools import PDBtools


def cal_angle_from_points(dist_lp, dist_dh, dist_ah):
    cos_theta = (np.power(dist_dh, 2) + np.power(dist_ah, 2) -
                 np.power(dist_lp, 2)) / (2 * dist_dh * dist_ah)
    theta = np.arccos(cos_theta)

    return theta


def cal_angle_from_vectors(A, B):
    cos_theta = np.dot(A, B)
    theta = np.arccos(cos_theta)
    return theta


class Pharmacophore(object):
    """
    pharmacophore
    """
    feature_type_dict = {
        'HBD': 0,
        'HBA': 1,
        'Anion': 2,
        'Cation': 3,
        'Hydrophobic': 4,
        'Aromatic': 5
    }
    feature_type_list = ['HBD', 'HBA', 'Anion',
                         'Cation', 'Hydrophobic', 'Aromatic']

    complementary_feature = {
        'HBD': [['HBA', 3.4]],
        'HBA': [['HBD', 3.4]],
        'Anion': [['Cation', 5.0]],
        'Cation': [['Anion', 5.0], ['Aromatic', 5.0]],
        'Hydrophobic': [['Hydrophobic', 6.5]],
        'Aromatic': [['Cation', 5.0], ['Aromatic', 6.0]],
        'MBA': [['Metal', 3.1]],
        'Metal': [['MBA', 3.1]],
    }

    num_feature_type = len(feature_type_list)
    SMARTS_pattern_dict = {
        'HBD': {
            0: '[N&H2&v3]',
            1: '[N&H1&v3]',
            2: '[N&!H0&+1&v4]',
            3: '[n&H1&+0]',
            4: '[n&H1&+1]',
            5: '[O;H1;+0]',
            6: '[S;H1;+0]',
            7: '[O;H2;+0]',
        },
        'HBA': {
            0: '[$([O,S;-])]',
            1: '[$([O,S;H0;v2])]',
            2: '[$([O,S;H1;v2]-[!$(*=[O,N,P,S])])]',
            3: '[$([N;v3;!$(N-*=!@[O,N,P,S])])]',
            4: '[$([nH0,o,s;+0])]',
            5: '[$([N;H0]#[C&v4])]',
            6: '[N&v3;H0;$(Nc)]',
            7: '[F;$(F-[#6])]',
            8: '[O;H2;+0]',

        },
        'Anion': {
            0: '[-]',
            1: '[SX4](=O)(=O)(-[O;H1,H0&-1])',
            2: '[PX4](=O)(-[O;H1,H0&-1])([!O])([!O])',
            3: '[PX4](=O)(-[O;H1,H0&-1])(-[O;H1,H0&-1])',
            4: '[CX3,SX3](=[O,S,P])-[O;H1,H0&-1]',
        },
        'Cation': {
            0: '[+]',
            1: '[NX3]=[CX3]([NX3])[!N]',
            2: 'NC(=N)N',
            3: 'c1ncnc1'
        },
        'Hydrophobic': {
            0: '[D3,D4;#6;+0;!$([#6][#7,#8,#9])]',
            1: '[R0;D2;#6;+0;!$([#6][#7,#8,#9])]',
            2: '[CX4](F)(F)(F)',
            3: '[#6;R]',
            4: '[#17,#35,#53]',
        },
        'Aromatic': {
            0: 'a1:a:a:a:1',
            1: 'a1:a:a:a:a:1',
            2: 'a1:a:a:a:a:a:1',
            3: 'a1:a:a:a:a:a:a:1',
            4: 'a1:a:a:a:a:a:a:a:1',
        }
    }

    metal_name_dict = {
        'MG': 0,
        'K': 1,
        'MN': 2,
        'FE': 3,
        'ZN': 4,
    }
    metal_type_list = ['MG', 'K', 'MN', 'FE', 'ZN']
    metal_atomic_name_dict = {12: 'MG', 19: 'K', 25: 'MN', 26: 'FE', 30: 'ZN'}
    metal_bind_name_dict = {
        7: 0,
        8: 1,
        16: 2,
    }
    metal_interaction_cutoff = {
        0: {2: 2.6, 3: 2.5, 4: 2.56},
        1: {0: 3.1, 1: 2.4, 2: 2.7, 3: 2.5, 4: 2.8},
        2: {4: 2.6},
    }

    weight_mcis = {
        'HBD:HBA': 0.129,
        'HBA:HBD': 0.129,
        'Anion:Cation': 0.31,
        'Cation:Anion': 0.31,
        'Cation:Aromatic': 0.039,
        'Hydrophobic:Hydrophobic': 0.00585,
        'Aromatic:Cation': 0.039,
        'Aromatic:Aromatic': 0.064,
        'MBA:Metal': 1.0
    }
    bias_mcis = 4.0
    cutoff = 6.5
    use_mcinfo = False

    def __init__(self, params=None):
        """
            initialize Pharmacophore
            set parameters
        """

        self.interaction_feature = list(self.weight_mcis.keys())
        if 'weight_mcis' in params:
            self.weight_mcis = params['weight_mcis']
        if 'bias_mcis' in params:
            self.bias_mcis = params['bias_mcis']
        if 'cutoff' in params:
            self.cutoff = params['cutoff']
        cmin = -99999.999
        cmax = 99999.999

        use_box = False
        if 'dock_config' in params:
            box_center, box_size = self.read_dock_config(params['dock_config'])
            box_min = box_center - box_size/2
            box_max = box_center + box_size/2
            cmin = box_min - self.cutoff
            cmax = box_max + self.cutoff
            use_box = True

        if 'ligand_file' in params:
            ligand_file = params['ligand_file']
            ligand_model_dict = PDBtools.read_coor_pdb(
                ligand_file, exclude_Hs=True)
            ligand_dict = ligand_model_dict[1]
            ligand_coor = list(ligand_dict.values())[0]
            cmin, cmax = PDBtools.cal_ligand_size(ligand_coor)
            cmin = cmin - self.cutoff
            cmax = cmax + self.cutoff
            use_box = True

        if 'include_hydrophobic' in params:
            include_hydrophobic = params['include_hydrophobic']

        mciscore_receptor = params['mciscore_receptor']
        pf_receptor = params['pf_receptor']

        if mciscore_receptor is not None:
            ms = pybel.readfile('pdb', mciscore_receptor)
            m_protein = list(ms)[0]
            receptor_bond_dict = self.get_bond_info(m_protein)
            PF_dict_protein = self.find_PF(m_protein, receptor_bond_dict,
                                           is_protein=True)
            PF_coor_protein = self.find_PF_coor(m_protein, PF_dict_protein,
                                receptor_bond_dict, is_protein=True)
            metal_coor = self.find_metal_ion(m_protein, is_protein=True)
            PF_coor_protein['Metal'] = metal_coor
            if use_box:
                self.PF_coor_protein = self.box_protein(
                    PF_coor_protein, cmin, cmax)
            else:
                self.PF_coor_protein = PF_coor_protein
            if pf_receptor is not None:
                PF_coor_protein_dict = {0: self.PF_coor_protein}
                self.write_PF(PF_coor_protein_dict,
                              pf_receptor, is_protein=True)
        elif pf_receptor is not None:
            self.PF_coor_protein = self.read_PF(pf_receptor,
                                                is_protein=True)[0]

        mcinfo_ligand = params['mcinfo_ligand']
        pf_receptor_info = params['pf_receptor_info']
        if mcinfo_ligand is not None:
            self.use_mcinfo = True
            template_dict = self.find_template(self.PF_coor_protein, mcinfo_ligand,
                                               include_hydrophobic=include_hydrophobic)
            self.PF_coor_info = self.select_template_protein(self.PF_coor_protein,
                                                             template_dict, cut_num=1)
            if pf_receptor_info is not None:
                PF_coor_info_dict = {0: self.PF_coor_info}
                self.write_PF(PF_coor_info_dict, pf_receptor_info,
                              is_protein=True)
        elif pf_receptor_info is not None:
            self.use_mcinfo = True
            self.PF_coor_info = self.read_PF(
                pf_receptor_info, is_protein=True)[0]

    def get_bond_info(self, m):
        """
            find neighbor atoms in molecule
            input:
                m: pybel molecule object
            output:
                bond_dict: neighbor atom dictionary
        """
        bond_dict = dict()
        mol = m.OBMol

        for i in range(mol.NumBonds()):
            bb = mol.GetBondById(i)
            if bb is None:
                continue
            begin = bb.GetBeginAtomIdx()
            end = bb.GetEndAtomIdx()
            if begin not in bond_dict:
                bond_dict[begin] = [end]
            else:
                bond_dict[begin] += [end]
            if end not in bond_dict:
                bond_dict[end] = [begin]
            else:
                bond_dict[end] += [begin]
        return bond_dict

    def find_PF(self, m, bond_dict, is_protein=False):
        """
            Find pharmacophoric features
            input:
                m: pybel molecule object
                is_protein: protein or not
            output: PF dictionary
        """
        PF_dict = dict()
        atoms = m.atoms
#        bond_dict = self.get_bond_info(m)

        # search pharmacophoric feature using SMARTS pattern
        for feature_type in self.feature_type_list:
            patterns = self.SMARTS_pattern_dict[feature_type]
            PF_dict[feature_type] = list()
            pattern_keys = sorted(patterns.keys())
            for pattern_idx in pattern_keys:
                smarts_pattern = patterns[pattern_idx]
                smarts = pybel.Smarts(smarts_pattern)
                atom_idx_group_list = smarts.findall(m)
                Nfeatures = len(atom_idx_group_list)
                for k in range(Nfeatures):
                    atom_idx_group = atom_idx_group_list[k]

                    fea_dict = dict()
                    fea_dict['atom_idx_group'] = atom_idx_group
                    fea_dict['pattern_idx'] = pattern_idx
                    if feature_type == 'HBD':
                        h_atoms_idx = list()
                        atom_idx = atom_idx_group[0]
                        neighbor_atoms_idx = bond_dict[atom_idx]
                        for neighbor_atom_idx in neighbor_atoms_idx:
                            neighbor_atom = atoms[neighbor_atom_idx - 1]
                            if neighbor_atom.atomicnum != 1:
                                continue
                            h_atoms_idx += [neighbor_atom_idx]
                        if len(h_atoms_idx) == 0:
                            continue
                        fea_dict['h_atoms_idx'] = h_atoms_idx

                    atom_idx = atom_idx_group[0]
                    atom = atoms[atom_idx-1]
                    fea_dict['atom_type'] = atom.type
                    if is_protein:
                        residue = atom.residue
                        chain_id = residue.OBResidue.GetChain()
                        fea_dict['chain_id'] = chain_id
                        residue_name = residue.name
                        fea_dict['residue_name'] = residue_name
                        residue_num = residue.OBResidue.GetNum()
                        fea_dict['residue_num'] = residue_num
                        fea_dict['weight'] = 1.0

                    PF_dict[feature_type] += [fea_dict]

        return PF_dict

    def find_PF_coor(self, m, PF_dict, bond_dict, is_protein=False):
        """
            Find pharmacophoric features
            input:
                m: pybel molecule object
                PF_dict: pharmacophoric feature dictionary
                is_protein: True or False
            output:
                PF_coor: dictionary
        """

        atoms = m.atoms
        PF_coor = dict()
        for feature_type in self.feature_type_list:
            if feature_type not in PF_dict:
                continue
            fea_dict_list = PF_dict[feature_type]
            PF_coor[feature_type] = list()
            Nfeatures = len(fea_dict_list)
            p00_atoms_old = list()
            num_member = list()
            for k in range(Nfeatures):
                fea_dict = fea_dict_list[k]
                num_member += [len(fea_dict['atom_idx_group'])]
            num_member = -np.array(num_member)
            index = num_member.argsort()
            for k in index:
                fea_dict = fea_dict_list[k]
                atom_idx_group = fea_dict['atom_idx_group']
                pattern_idx = fea_dict['pattern_idx']

                intersection = set(atom_idx_group).intersection(p00_atoms_old)
                p00_atoms_old += atom_idx_group
                if len(intersection) > 0 and feature_type != 'Aromatic':
                    continue

                pseudo_atom = []
                for atom_index in atom_idx_group:
                    coor = np.array(atoms[atom_index - 1].coords)
                    pseudo_atom += [coor]
                pseudo_atom = np.array(pseudo_atom)
                pseudo_atom_coor = pseudo_atom.mean(axis=0)

                fea_dict_new = dict()
                fea_dict_new['atom_idx_group'] = atom_idx_group
                fea_dict_new['pattern_idx'] = pattern_idx
                fea_dict_new['pseudo_atom_coor'] = pseudo_atom_coor

                if feature_type == 'HBD':
                    h_atoms = list()
                    atom_idx = atom_idx_group[0]
                    neighbor_heavy_atom_list = list()
                    neighbor_atoms_idx = bond_dict[atom_idx]
                    neighbor_atom_list = list()
                    for neighbor_atom_idx in neighbor_atoms_idx:
                        neighbor_atom = atoms[neighbor_atom_idx - 1]
                        if neighbor_atom.atomicnum == 1:
                            continue
                        nei_coor = np.array(atoms[neighbor_atom_idx - 1].coords)
                        neighbor_atom_list += [[neighbor_atom_idx, nei_coor]]

                    h_atoms_idx = fea_dict['h_atoms_idx']
                    num_neighbor_atom = len(neighbor_atom_list)
                    for h_atom_idx in h_atoms_idx:
                        h_coor = np.array(atoms[h_atom_idx - 1].coords)
                        h_atoms += [[h_atom_idx, h_coor]]
                        #if len(neighbor_atom_list)==1:
                            #nei_coor = neighbor_atom_list[0][1]
                            #vec_a = pseudo_atom_coor - nei_coor
                            #vec_b = h_coor - pseudo_atom_coor
                            #vec_c = h_coor - nei_coor
                            #a = np.linalg.norm(vec_a)
                            #b = np.linalg.norm(vec_b)
                            #c = np.linalg.norm(vec_c)
                            #d = (c**2 - a**2 - b**2)/(2*a**2)
                            #m_coor = pseudo_atom_coor + vec_a * d/a
                            #h_coor2 = h_coor + 2 * (m_coor - h_coor)
                            #h_atoms += [[h_atom_idx, h_coor2]]

                    fea_dict_new['h_atoms'] = h_atoms
                    fea_dict_new['num_neighbor_atom'] = num_neighbor_atom



                elif feature_type == 'Aromatic':
                    v_a = pseudo_atom[1] - pseudo_atom[0]
                    v_b = pseudo_atom[2] - pseudo_atom[1]
                    v_c = np.cross(v_a, v_b)
                    v_n = v_c / np.linalg.norm(v_c)
                    fea_dict_new['v_n'] = v_n
                fea_dict_new['atom_type'] = fea_dict['atom_type']
                if is_protein:
                    fea_dict_new['chain_id'] = fea_dict['chain_id']
                    fea_dict_new['residue_name'] = fea_dict['residue_name']
                    fea_dict_new['residue_num'] = fea_dict['residue_num']
                    fea_dict_new['weight'] = fea_dict['weight']

                PF_coor[feature_type] += [fea_dict_new]

        return PF_coor

    def find_metal_ion(self, m, is_protein=True):
        """
            input:
                m: pybel molecule object
                is_protein: True or False
            output:
                metal_coor: list
        """
        metal_coor = list()
        atoms = m.atoms
        for atom in atoms:
            atomic_num = atom.atomicnum

            if atomic_num not in self.metal_atomic_name_dict:
                continue
            atomic_name = self.metal_atomic_name_dict[atomic_num]
            atom_idx = atom.idx
            atom_idx_group = (atom_idx,)
            coor = np.array(atom.coords)
            if atomic_name not in self.metal_name_dict:
                continue
            pattern_idx = self.metal_name_dict[atomic_name]

            fea_dict = dict()
            fea_dict['atom_type'] = atom.type
            fea_dict['atom_idx_group'] = atom_idx_group
            fea_dict['pattern_idx'] = pattern_idx
            fea_dict['pseudo_atom_coor'] = coor
            if is_protein:
                residue = atom.residue
                chain_id = residue.OBResidue.GetChain()
                fea_dict['chain_id'] = chain_id
                residue_name = residue.name
                fea_dict['residue_name'] = residue_name
                residue_num = residue.OBResidue.GetNum()
                fea_dict['residue_num'] = residue_num
                fea_dict['weight'] = 1.0

            metal_coor += [fea_dict]

        return metal_coor

    def find_mba(self, m, is_protein=False):
        """
            metal binding atom for ligand
            input:
                m: pybel molecule object
                is_protein: True or False
            output:
                mba_coor: list
        """

        mba_coor = list()
        atoms = m.atoms
        for atom in atoms:
            atom_idx = atom.idx
            atomic_num = atom.atomicnum
            if atomic_num != 7 and atomic_num != 8 and atomic_num != 16:
                continue

            # if atom.formalcharge>0:
            #    continue
            # Excluded because protonation state calculation is incomplete

            atom_idx_group = (atom_idx,)
            coor = np.array(atom.coords)
            pattern_idx = self.metal_bind_name_dict[atomic_num]

            fea_dict = dict()
            fea_dict['atom_type'] = atom.type
            fea_dict['atom_idx_group'] = atom_idx_group
            fea_dict['pattern_idx'] = pattern_idx
            fea_dict['pseudo_atom_coor'] = coor
            if is_protein:
                residue = atom.residue
                chain_id = residue.OBResidue.GetChain()
                fea_dict['chain_id'] = chain_id
                residue_name = residue.name
                fea_dict['residue_name'] = residue_name
                residue_num = residue.OBResidue.GetNum()
                fea_dict['residue_num'] = residue_num
                fea_dict['weight'] = 1.0
            mba_coor += [fea_dict]

        return mba_coor

    def draw_ligand(self, m_ligand, PF_dict, output_name, size):
        type_list = self.feature_type_list
        patom_list = list()
        for fea in type_list:
            patoms = list()
            fff = PF_dict[fea]
            for ff in fff:
                atom_idx_group = ff[0]
                for idx in atom_idx_group:
                    patoms += [idx-1]
            patom_list += [patoms]
        m_ligand_noH = Chem.RemoveHs(m_ligand)
        m_ligand_noH.RemoveConformer(0)
        mol_list = [m_ligand_noH]*6
        img = Draw.MolsToGridImage(mol_list, legends=type_list,
                                   highlightAtomLists=patom_list,
                                   subImgSize=size, molsPerRow=2)
        img.save(output_name)

    def box_protein(self, PF_coor_protein, cmin, cmax):
        """
            select PF of PF_coor_protein in the box
            input:
                PF_coor_protein: dict
                cmin : np.array
                cmax : np.array
            output:
                PF_coor_box: dict
        """

        type_list = self.feature_type_list + ['Metal']
        PF_coor_box = dict()
        for fea in type_list:
            if fea not in PF_coor_protein:
                continue
            fea_dict_list = PF_coor_protein[fea]
            PF_coor_box[fea] = list()
            for fea_dict in fea_dict_list:
                pseudo_atom_coor = fea_dict['pseudo_atom_coor']
                if (pseudo_atom_coor > cmax).any():
                    continue
                elif (pseudo_atom_coor < cmin).any():
                    continue
                PF_coor_box[fea] += [fea_dict]
        return PF_coor_box

    def find_template(self, PF_coor_protein, template_ligand_file, include_hydrophobic=False):
        template_dict = dict()
        feature_type_list_receptor = self.feature_type_list + ['Metal']
        for feature_type_receptor in feature_type_list_receptor:
            if not include_hydrophobic and feature_type_receptor == 'Hydrophobic':
                continue
            template_dict[feature_type_receptor] = dict()

        file_format = template_ligand_file.split('.')[-1]
        ms = list(pybel.readfile(file_format, template_ligand_file))
        m_ligand = ms[0]
        ligand_bond_dict = self.get_bond_info(m_ligand)
        PF_dict_ligand = self.find_PF(m_ligand, ligand_bond_dict)

        num_model = len(ms)
        if num_model > 1:
            num_model = num_model - 1
        PF_coor_ligand_dict = dict()
        interaction_dict = dict()
        model_idx = 0
        m_ligand = ms[model_idx]
        PF_coor_ligand = self.find_PF_coor(m_ligand, PF_dict_ligand, ligand_bond_dict)
        mba_coor = self.find_mba(m_ligand, is_protein=False)
        PF_coor_ligand['MBA'] = mba_coor

        PF_coor_ligand_dict[model_idx] = PF_coor_ligand
        interaction = self.find_interaction(PF_coor_protein, PF_coor_ligand,
                                            ligand_bond_dict)
        interaction_dict[model_idx] = interaction
        total_rec_pf_dict = self.count_interaction_receptor(interaction_dict)
        rec_pf_dict = total_rec_pf_dict[0]
        for feature_type_receptor in rec_pf_dict.keys():
            if not include_hydrophobic and feature_type_receptor == 'Hydrophobic':
                continue
            rec_pf = rec_pf_dict[feature_type_receptor]
            atom_idx_group_list = rec_pf.keys()
            for atom_idx_group in atom_idx_group_list:
                if atom_idx_group not in template_dict[feature_type_receptor]:
                    template_dict[feature_type_receptor][atom_idx_group] = 0
                template_dict[feature_type_receptor][atom_idx_group] += 1
        return template_dict

    def select_template_protein(self, PF_coor_protein, template_dict,
                                cut_num=0):
        """
            select PF of PF_coor_protein in template
            input:
                PF_coor_protein: dict
                template_dict: dict
            output:
                PF_coor_info: dict
        """
        type_list = self.feature_type_list + ['Metal']
        PF_coor_info = dict()
        for fea in type_list:
            if fea not in PF_coor_protein:
                continue
            if fea not in template_dict:
                continue
            fea_dict_list = PF_coor_protein[fea]
            PF_coor_info[fea] = list()
            for fea_dict in fea_dict_list:
                atom_idx_group = fea_dict['atom_idx_group']
                if atom_idx_group not in template_dict[fea]:
                    continue
                num_in_template = template_dict[fea][atom_idx_group]
                if num_in_template < cut_num:
                    continue
                PF_coor_info[fea] += [fea_dict]
        return PF_coor_info

    def find_interaction(self, PF_coor_protein, PF_coor_ligand,
                         ligand_bond_dict):
        """
            find pharmacophoric interaction
            input:
                PF_coor_protein: dict
                PF_coor_ligand: dict
                ligand_bond_dict: dict
            output:
                interaction: dict
        """
        type_list = self.feature_type_list + ['MBA']
#        type_list = self.feature_type_list

        coor_ligand_list = list()
        for fea_ligand in type_list:
            if fea_ligand not in PF_coor_ligand:
                continue
            fea_dict_list_ligand = PF_coor_ligand[fea_ligand]
            for fea_dict_ligand in fea_dict_list_ligand:
                coor_ligand = fea_dict_ligand['pseudo_atom_coor']
                coor_ligand_list.append(coor_ligand)
        coor_ligand_list = np.array(coor_ligand_list)
        cmin = coor_ligand_list.min(axis=0) - self.cutoff
        cmax = coor_ligand_list.max(axis=0) + self.cutoff

        PF_coor_protein_box = self.box_protein(PF_coor_protein, cmin, cmax)

        interaction = dict()
        for fea_ligand in type_list:
            if fea_ligand not in interaction:
                interaction[fea_ligand] = dict()
            if fea_ligand not in PF_coor_ligand:
                continue
            fea_dict_list_ligand = PF_coor_ligand[fea_ligand]
            cfea_list = self.complementary_feature[fea_ligand]
            for fea_dict_ligand in fea_dict_list_ligand:
                idx_ligand = fea_dict_ligand['atom_idx_group']
                pattern_ligand = fea_dict_ligand['pattern_idx']
                coor_ligand = fea_dict_ligand['pseudo_atom_coor']
                for cfea_t in cfea_list:
                    fea_protein = cfea_t[0]
                    if fea_protein not in interaction[fea_ligand]:
                        interaction[fea_ligand][fea_protein] = list()
                    pair_cutoff = cfea_t[1]
                    if fea_protein not in PF_coor_protein_box:
                        continue
                    fea_dict_list_protein = PF_coor_protein_box[fea_protein]
                    for fea_dict_protein in fea_dict_list_protein:
                        idx_protein = fea_dict_protein['atom_idx_group']
                        pattern_protein = fea_dict_protein['pattern_idx']
                        coor_protein = fea_dict_protein['pseudo_atom_coor']

                        r_lp = coor_protein - coor_ligand
                        dist_lp = np.linalg.norm(r_lp)
#                        if dist_lp > self.cutoff:
                        if dist_lp > pair_cutoff:
                            continue

                        interact_ij = dict()
                        interact_ij['idx_ligand'] = idx_ligand
                        interact_ij['idx_protein'] = idx_protein
                        interact_ij['pattern_ligand'] = pattern_ligand
                        interact_ij['pattern_protein'] = pattern_protein

                        interact_ij['dist_lp'] = dist_lp
                        interact_ij['chain_id'] = fea_dict_protein['chain_id']
                        interact_ij['residue_name'] = fea_dict_protein['residue_name']
                        interact_ij['residue_num'] = fea_dict_protein['residue_num']
                        interact_ij['weight'] = fea_dict_protein['weight']

                        if fea_ligand == 'HBD' and fea_protein == 'HBA':
                            num_neighbor_atom = fea_dict_ligand['num_neighbor_atom']
                            donor_hydrogen_list = fea_dict_ligand['h_atoms']
                            if num_neighbor_atom <= 1:
                                theta = 0
                                interact_ij['theta'] = theta
                            else:
                                dist_ah = 100.0
                                for donor_hydrogen in donor_hydrogen_list:
                                    hydrogen_coor0 = donor_hydrogen[1]
                                    dist_ah0 = np.linalg.norm(
                                        hydrogen_coor0 - coor_protein)
                                    if dist_ah0 < dist_ah:
                                        dist_ah = dist_ah0
                                        hydrogen_coor = hydrogen_coor0
                                dist_dh = np.linalg.norm(
                                    hydrogen_coor - coor_ligand)
                                theta = cal_angle_from_points(
                                    dist_lp, dist_dh, dist_ah)
                                if theta < 0.6*np.pi:
                                    continue
                                interact_ij['theta'] = theta

                        elif fea_ligand == 'HBA' and fea_protein == 'HBD':
                            num_neighbor_atom = fea_dict_protein['num_neighbor_atom']
                            donor_hydrogen_list = fea_dict_protein['h_atoms']
                            if num_neighbor_atom <= 1:
                                theta = 0
                                interact_ij['theta'] = theta
                            else:
                                dist_ah = 100.0
                                for donor_hydrogen in donor_hydrogen_list:
                                    hydrogen_coor0 = donor_hydrogen[1]
                                    dist_ah0 = np.linalg.norm(
                                        hydrogen_coor0 - coor_ligand)
                                    if dist_ah0 < dist_ah:
                                        dist_ah = dist_ah0
                                        hydrogen_coor = hydrogen_coor0
                                dist_dh = np.linalg.norm(
                                    hydrogen_coor - coor_protein)
                                theta = cal_angle_from_points(
                                    dist_lp, dist_dh, dist_ah)
                                if theta < 0.6*np.pi:
                                    continue
                                interact_ij['theta'] = theta

                        elif fea_ligand == 'Aromatic' and fea_protein == 'Cation':
                            n_vector_ligand = fea_dict_ligand['v_n']
                            n_vector_lp = r_lp/dist_lp
                            theta = np.abs(cal_angle_from_vectors(
                                           n_vector_ligand, n_vector_lp))
                            if theta > np.pi/2:
                                theta = np.pi-theta
                            if theta > 0.5*np.pi/2:
                                continue
                            interact_ij['theta'] = theta

                        elif fea_ligand == 'Cation' and fea_protein == 'Aromatic':
                            n_vector_protein = fea_dict_protein['v_n']
                            n_vector_lp = r_lp/dist_lp
                            theta = np.abs(cal_angle_from_vectors(
                                           n_vector_protein, n_vector_lp))
                            if theta > np.pi/2:
                                theta = np.pi-theta
                            if theta > 0.5*np.pi/2:
                                continue
                            interact_ij['theta'] = theta

                        elif fea_ligand == 'Aromatic' and fea_protein == 'Aromatic':
                            n_vector_ligand = fea_dict_ligand['v_n']
                            n_vector_protein = fea_dict_protein['v_n']
                            n_vector_lp = r_lp/dist_lp
                            theta = cal_angle_from_vectors(
                                n_vector_ligand, n_vector_protein)
                            if theta > np.pi/2:
                                theta = np.pi-theta
                            interact_ij['theta'] = theta

                            alpha = np.abs(cal_angle_from_vectors(
                                           n_vector_ligand, n_vector_lp))
#                            dist_v = dist_lp * np.cos(alpha)
#                            dist_h = dist_lp * np.sin(alpha)
#                            print(np.sqrt(dist_v*dist_v+dist_h*dist_h), dist_lp)
#                            interact_ij['dist_v'] = dist_v
#                            interact_ij['dist_h'] = dist_h
                            if theta > 0.4*np.pi/2 and theta < 0.8*np.pi/2:
                                continue
                            interact_ij['alpha'] = alpha

                        elif fea_ligand == 'MBA' and fea_protein == 'Metal':
                            if pattern_ligand not in self.metal_interaction_cutoff:
                                continue
                            if pattern_protein not in self.metal_interaction_cutoff[pattern_ligand]:
                                continue
                            metal_cutoff = self.metal_interaction_cutoff[
                                pattern_ligand][pattern_protein]
                            if dist_lp > metal_cutoff:
                                continue
                        interaction[fea_ligand][fea_protein] += [interact_ij]

        if 'MBA' in interaction:
            if 'Metal' in interaction['MBA']:
                metal_interaction = interaction['MBA']['Metal']
                d_list = list()
                d_list = np.array([x['dist_lp'] for x in metal_interaction])
                idx_sorted = np.argsort(d_list)

                metal_interaction_new = list()
                batom_idx = list()
                for idx in idx_sorted:
                    interact_ij = metal_interaction[idx]
                    idx_ligand = interact_ij['idx_ligand'][0]
#                    dist_lp = interact_ij['dist_lp']
                    neighbor_atoms_idx = ligand_bond_dict[idx_ligand]
                    neighbor_check = False
                    for neighbor_atom_idx in neighbor_atoms_idx:

                        if neighbor_atom_idx in batom_idx:
                            neighbor_check = True
                    if not neighbor_check:
                        metal_interaction_new += [interact_ij]
                        batom_idx += [idx_ligand]

                interaction['MBA']['Metal'] = metal_interaction_new

        return interaction

    def write_PF(self, PF_coor_dict, out_file, is_protein=False):
        """
            write pharmacophoric feature
            input:
                PF_coor_dict: dict
                out_file: filename, str
                is_protein: bool
        """

        if is_protein:
            feature_type_list = self.feature_type_list + ['Metal']
        else:
            feature_type_list = self.feature_type_list + ['MBA']
        fp = open(out_file, 'w')
        line_out = 'feature_type:atom_idx_group:pattern_idx:pseudo_atom_coor'
        line_out += ':atom_type:etc'
        if is_protein:
            line_out += ':chain_id:residue_name:residue_num:weight'
        line_out += '\n'
        fp.write(line_out)
        model_idx_list = PF_coor_dict.keys()
        num_model = len(model_idx_list)
        for model_idx in model_idx_list:
            PF_coor = PF_coor_dict[model_idx]
            if num_model > 1:
                line_out = 'MODEL %d\n' % (model_idx + 1)
                fp.write(line_out)
            for feature_type in feature_type_list:
                if feature_type not in PF_coor:
                    continue
                fea_dict_list = PF_coor[feature_type]
                for fea_dict in fea_dict_list:
                    line_out = '%s' % feature_type
#                    line_out += ':%s' % (str(fea_dict['atom_idx_group']))
                    idx_line = ''
                    for atom_idx in fea_dict['atom_idx_group']:
                        idx_line += '%d,' % atom_idx
                    line_out += ':(%s)' % idx_line.strip(',')

                    line_out += ':%s' % (fea_dict['pattern_idx'])
                    coor = fea_dict['pseudo_atom_coor']
                    line_out += ':(%.3f,%.3f,%.3f)' % (
                        coor[0], coor[1], coor[2])
                    line_out += ':%s' % (fea_dict['atom_type'])

                    if 'h_atoms' in fea_dict:
                        hatoms = fea_dict['h_atoms']
                        hatom_line = ''
                        if 'num_neighbor_atom' in fea_dict:
                            num_neighbor_atom = fea_dict['num_neighbor_atom']
                            hatom_line += '%d;' %num_neighbor_atom
                        for hatom in hatoms:
                            line_f = '%d=(%.3f,%.3f,%.3f);'
                            hatom_line += line_f % (hatom[0], hatom[1][0],
                                                    hatom[1][1], hatom[1][2])
                        line_out += ':%s' % hatom_line.strip(';')

                    elif 'v_n' in fea_dict:
                        coor_v_n = fea_dict['v_n']
                        line_out += ':(%.3f,%.3f,%.3f)' % (
                            coor_v_n[0], coor_v_n[1], coor_v_n[2])

                    else:
                        line_out += ':'

                    if is_protein:
                        line_out += ':%s' % (fea_dict['chain_id'])
                        line_out += ':%s' % (fea_dict['residue_name'])
                        line_out += ':%d' % (fea_dict['residue_num'])
                        line_out += ':%.3f' % (fea_dict['weight'])

                    line_out += '\n'
                    fp.write(line_out)
        fp.close()

    def read_PF(self, feature_file, is_protein=False):
        PF_coor_dict = dict()
        fp = open(feature_file)
        lines = fp.readlines()
        fp.close()

        model_idx = 0
        PF_coor_dict[model_idx] = dict()

        for line in lines:
            if line.startswith('feature_type'):
                # title = line.strip().split(':')
                continue

            if line.startswith('MODEL'):
                model_idx = int(line.strip().split()[1]) - 1
                fea_dict = dict()
                if model_idx not in PF_coor_dict:
                    PF_coor_dict[model_idx] = dict()
                continue
            lis = line.strip().split(':')
            if len(lis) < 1:
                continue
            fea_dict = dict()
            feature_type = lis[0]
            atom_idx_group = lis[1]
            pattern_idx = lis[2]
            pseudo_atom_coor = lis[3]
            atom_type = lis[4]
            etc = lis[5]
            if feature_type not in PF_coor_dict[model_idx]:
                PF_coor_dict[model_idx][feature_type] = list()

            fea_dict['feature_type'] = feature_type
            aline = atom_idx_group.lstrip('(').rstrip(')')
            fea_dict['atom_idx_group'] = np.array(aline.split(','), dtype=int)
            fea_dict['pattern_idx'] = int(pattern_idx)
            coor = pseudo_atom_coor.lstrip('(').rstrip(')').split(',')
            fea_dict['pseudo_atom_coor'] = np.array(coor, dtype=np.float32)
            fea_dict['atom_type'] = atom_type
            if feature_type == 'HBD':
                hline_list = etc.split(';')
                h_list = list()
                num_neighbor_atom = float(hline_list[0])
                fea_dict['num_neighbor_atom'] = num_neighbor_atom
                for hline in hline_list[1:]:
                    aa = hline.split('=')
                    h_idx = int(aa[0])
                    h_coor = aa[1].lstrip('(').rstrip(')').split(',')
                    h_list += [[h_idx, np.array(h_coor, dtype=np.float32)]]
                    fea_dict['h_atoms'] = h_list
            if feature_type == 'Aromatic':
                v_n = etc.lstrip('(').rstrip(')').split(',')
                fea_dict['v_n'] = np.array(v_n, dtype=np.float32)
            if is_protein:
                chain_id = lis[6]
                residue_name = lis[7]
                residue_num = lis[8]
                fea_dict['chain_id'] = chain_id
                fea_dict['residue_name'] = residue_name
                fea_dict['residue_num'] = int(residue_num)
                if len(lis) >= 10:
                    weight = lis[9]
                    fea_dict['weight'] = float(weight)
                else:
                    fea_dict['weight'] = 1.0

            PF_coor_dict[model_idx][feature_type] += [fea_dict]

        return PF_coor_dict

    def write_interaction(self, interaction_dict, out_file):
        """
            write pharmacophoric interaction
            input:
                interaction_dict: dict()
                out_file: str
        """
        fp = open(out_file, 'w')
        line_out = 'feature_type_ligand:feature_type_protein'
        line_out += ':atom_idx_group_ligand:atom_idx_group_protein'
        line_out += ':pattern_ligand:pattern_protein'
        line_out += ':chain_id:residue_name:residue_num'
        line_out += ':dist_lp:theta:alpha:weight\n'
        fp.write(line_out)

        model_idx_list = interaction_dict.keys()
        num_model = len(model_idx_list)
        for model_idx in model_idx_list:
            interaction = interaction_dict[model_idx]
            if num_model > 1:
                line_out = 'MODEL %d\n' % (model_idx + 1)
                fp.write(line_out)

            for feature_type_ligand in interaction.keys():
                inter_i = interaction[feature_type_ligand]
                for feature_type_protein in inter_i.keys():
                    inter_ij_list = inter_i[feature_type_protein]
                    for inter_ij in inter_ij_list:
                        line_out = '%s' % feature_type_ligand
                        line_out += ':%s' % feature_type_protein

                        idx_line = ''
                        for idx in inter_ij['idx_ligand']:
                            idx_line += '%d,' % (idx)
                        line_out += ':(%s)' % (idx_line.strip(','))

                        idx_line = ''
                        for idx in inter_ij['idx_protein']:
                            idx_line += '%d,' % (idx)
                        line_out += ':(%s)' % (idx_line.strip(','))
                        line_out += ':%s' % (inter_ij['pattern_ligand'])
                        line_out += ':%s' % (inter_ij['pattern_protein'])

                        line_out += ':%s' % (inter_ij['chain_id'])
                        line_out += ':%s' % (inter_ij['residue_name'])
                        line_out += ':%d' % (inter_ij['residue_num'])
                        line_out += ':%.3f' % (inter_ij['dist_lp'])
                        if 'theta' in inter_ij:
                            line_out += ':%.3f' % (inter_ij['theta'])
                        else:
                            line_out += ':'
                        if 'alpha' in inter_ij:
                            line_out += ':%.3f' % (inter_ij['alpha'])
                        else:
                            line_out += ':'
                        line_out += ':%.3f' % (inter_ij['weight'])

                        line_out += '\n'
                        fp.write(line_out)

        fp.close()

    def count_interaction_receptor(self, interaction_dict):

        total_rec_pf_dict = dict()
        model_idx_list = interaction_dict.keys()
        feature_type_list_receptor = self.feature_type_list + ['Metal']
        feature_type_list_ligand = self.feature_type_list + ['MBA']

        for model_idx in model_idx_list:
            line_out_title = 'Model_idx'
            interaction = interaction_dict[model_idx]

            rec_pf_dict = dict()
            for feature_type_ligand in feature_type_list_ligand:
                line_out_title += ' %s' % (feature_type_ligand)
                if feature_type_ligand not in interaction:
                    continue
                inter_i = interaction[feature_type_ligand]
                for feature_type_protein in feature_type_list_receptor:
                    if feature_type_protein not in rec_pf_dict:
                        rec_pf_dict[feature_type_protein] = dict()
                    if feature_type_protein not in inter_i.keys():
                        continue
                    inter_ij_list = inter_i[feature_type_protein]
                    for inter_ij in inter_ij_list:
                        idx_protein = inter_ij['idx_protein']
                        if idx_protein not in rec_pf_dict[feature_type_protein]:
                            rec_pf_dict[feature_type_protein][idx_protein] = 0
                        rec_pf_dict[feature_type_protein][idx_protein] += 1
            total_rec_pf_dict[model_idx] = rec_pf_dict

        return total_rec_pf_dict

    def print_interaction_receptor(self, total_rec_pf_dict):
        line_out_total = ''
        model_idx_list = total_rec_pf_dict.keys()
        feature_type_list_receptor = self.feature_type_list + ['Metal']

        line_out_title = 'Model_idx'
        for feature_type_receptor in feature_type_list_receptor:
            line_out_title += ' %s' % (feature_type_receptor)
        line_out_total += line_out_title + '\n'

        for model_idx in model_idx_list:
            line_out_value = '%d' % (model_idx + 1)
            rec_pf_dict = total_rec_pf_dict[model_idx]
            for feature_type_protein in feature_type_list_receptor:
                if feature_type_protein in rec_pf_dict:
                    count = len(rec_pf_dict[feature_type_protein].keys())
                else:
                    count = 0
                line_out_value += ' %d' % (count)
            line_out_total += line_out_value + '\n'

        return line_out_total

    def count_interaction_ligand(self, interaction_dict):
        total_lig_pf_dict = dict()
        model_idx_list = interaction_dict.keys()
        feature_type_list = self.feature_type_list + ['MBA']

        for model_idx in model_idx_list:
            interaction = interaction_dict[model_idx]
            lig_pf_dict = dict()
            for feature_type_ligand in feature_type_list:
                if feature_type_ligand not in interaction:
                    continue
                inter_i = interaction[feature_type_ligand]
                lig_pf_dict[feature_type_ligand] = dict()
                for feature_type_protein in inter_i.keys():
                    inter_ij_list = inter_i[feature_type_protein]
                    for inter_ij in inter_ij_list:
                        idx_ligand = inter_ij['idx_ligand']
                        if idx_ligand not in lig_pf_dict[feature_type_ligand]:
                            lig_pf_dict[feature_type_ligand][idx_ligand] = 0
                        lig_pf_dict[feature_type_ligand][idx_ligand] += 1
            total_lig_pf_dict[model_idx] = lig_pf_dict

        return total_lig_pf_dict

    def print_interaction_ligand(self, total_lig_pf_dict):
        line_out_total = ''
        model_idx_list = total_lig_pf_dict.keys()
        feature_type_list = self.feature_type_list + ['MBA']

        line_out_title = 'Model_idx'
        for feature_type_ligand in feature_type_list:
            line_out_title += ' %s' % (feature_type_ligand)
        line_out_total += line_out_title + '\n'

        for model_idx in model_idx_list:
            line_out_value = '%d' % (model_idx + 1)
            lig_pf_dict = total_lig_pf_dict[model_idx]
            for feature_type_ligand in feature_type_list:
                if feature_type_ligand in lig_pf_dict:
                    count = len(lig_pf_dict[feature_type_ligand].keys())
                else:
                    count = 0
                line_out_value += ' %d' % (count)
            line_out_total += line_out_value + '\n'

        return line_out_total

    def count_interaction(self, interaction_dict, use_weight_mcis=True,
                          include_hydrophobic=True):
        """
            input:
                interaction_dict: dict()
            output:
                total_count_dict: dict()
        """

        total_count_dict = dict()

        feature_type_list = self.feature_type_list + ['MBA']
        model_idx_list = interaction_dict.keys()
        for model_idx in model_idx_list:
            interaction = interaction_dict[model_idx]
            count_dict = dict()
            w_count_dict = dict()
            for feature_type_ligand in feature_type_list:
                cfea_list = self.complementary_feature[feature_type_ligand]
                for cfea in cfea_list:
                    feature_type_protein = cfea[0]

                    key = '%s:%s' % (feature_type_ligand, feature_type_protein)
                    if key not in count_dict:
                        count_dict[key] = 0
                        w_count_dict[key] = 0
                    if feature_type_ligand not in interaction:
                        continue
                    inter_i = interaction[feature_type_ligand]
                    if feature_type_protein not in inter_i:
                        continue
                    if ((not include_hydrophobic)
                            and feature_type_protein == 'Hydrophobic'):
                        continue

                    inter_ij_list = inter_i[feature_type_protein]
                    for inter_ij in inter_ij_list:
                        count_dict[key] += 1
                        weight = inter_ij['weight']
                        w_count_dict[key] += weight

            score = 0
            for key in self.interaction_feature:
                if use_weight_mcis:
                    score += w_count_dict[key]*self.weight_mcis[key]
                else:
                    score += w_count_dict[key]
            if use_weight_mcis:
                score += self.bias_mcis

            count_dict['Score'] = score

            total_count_dict[model_idx] = count_dict

        return total_count_dict

    def print_interaction(self, total_count_dict, total_count_info_dict,
                          use_mcinfo, print_option=1):
        """
            input:
                total_count_dict: from count_interaction
                print_option:
                    if 1, print interction in one line for one model.
                    if 2, print interaction in multi lines.
            output:
                line_out_total: str
        """

        line_out_total = ''

        model_idx_list = total_count_dict.keys()
        num_model = len(model_idx_list)
        if print_option == 1:
            line_out_title = 'Model_idx PIscore'
            if use_mcinfo:
                line_out_title += ' Pinfo'
            for key in self.interaction_feature:
                line_out_title += ' %s' % (key)
            line_out_total += line_out_title + '\n'

        for model_idx in model_idx_list:
            count_dict = total_count_dict[model_idx]
            if num_model > 1 and print_option == 2:
                line_out = 'MODEL %d' % (model_idx + 1)
                line_out_total += line_out + '\n'
            if print_option == 1:
                line_out_value = '%d' % (model_idx + 1)
                line_out_value += ' %.4f' % (count_dict['Score'])
                if use_mcinfo:
                    count_info_dict = total_count_info_dict[model_idx]
                    line_out_value += ' %.4f' % (count_info_dict['Score'])

            for key in self.interaction_feature:
                count = count_dict[key]
                if print_option == 1:
                    line_out_value += ' %d' % (count)

                if count > 0 and print_option == 2:
                    line_out = '%s %d\n' % (key, count)
                    line_out_total += line_out

            if print_option == 2:
                line_out = 'PIscore: %.4f\n' % (count_dict['Score'])
                if use_mcinfo:
                    count_info_dict = total_count_info_dict[model_idx]
                    line_out += 'Pinfo %.4f\n' % (count_info_dict['Score'])

                line_out_total += line_out

            if print_option == 1:
                line_out_total += line_out_value + '\n'

        return line_out_total

    def read_dock_config(self, dock_config):
        fp = open(dock_config)
        lines = fp.readlines()
        fp.close()

        for line in lines:
            lis = line.strip().split('=')
            if lis[0] == 'center_x':
                center_x = float(lis[1])
            elif lis[0] == 'center_y':
                center_y = float(lis[1])
            elif lis[0] == 'center_z':
                center_z = float(lis[1])
            elif lis[0] == 'size_x':
                size_x = float(lis[1])
            elif lis[0] == 'size_y':
                size_y = float(lis[1])
            elif lis[0] == 'size_z':
                size_z = float(lis[1])

        center = np.array((center_x, center_y, center_z), dtype=float)
        size = np.array((size_x, size_y, size_z), dtype=float)
        return center, size

    def cal_mciscore(self, ligand_file, pf_ligand_file, interaction_file):

        ms = list(pybel.readfile('pdb', ligand_file))
        num_model = len(ms)
        if num_model > 1:
            num_model = num_model - 1
        if num_model < 1:
            return None

        m_ligand = ms[0]

        PF_coor_ligand_dict = dict()
        interaction_dict = dict()
        if self.use_mcinfo:
            interaction_info_dict = dict()
        ligand_bond_dict = self.get_bond_info(m_ligand)
        PF_dict_ligand = self.find_PF(m_ligand, ligand_bond_dict)
        for model_idx in range(num_model):
            m_ligand = ms[model_idx]
            PF_coor_ligand = self.find_PF_coor(m_ligand, PF_dict_ligand, ligand_bond_dict)
            mba_coor = self.find_mba(m_ligand, is_protein=False)
            PF_coor_ligand['MBA'] = mba_coor
            PF_coor_ligand_dict[model_idx] = PF_coor_ligand
            interaction = self.find_interaction(self.PF_coor_protein,
                                                PF_coor_ligand,
                                                ligand_bond_dict)

            interaction_dict[model_idx] = interaction
            if self.use_mcinfo:
                interaction_info = self.find_interaction(self.PF_coor_info,
                                                         PF_coor_ligand,
                                                         ligand_bond_dict)
                interaction_info_dict[model_idx] = interaction_info

        if pf_ligand_file is not None:
            self.write_PF(PF_coor_ligand_dict, pf_ligand_file)

        if interaction_file is not None:
            self.write_interaction(interaction_dict, interaction_file)

        total_count_dict = self.count_interaction(interaction_dict,
                                                  use_weight_mcis=True,
                                                  include_hydrophobic=True)
        if self.use_mcinfo:
            total_count_info_dict = self.count_interaction(
                interaction_info_dict,
                use_weight_mcis=False,
                include_hydrophobic=True)
        else:
            total_count_info_dict = None
        return total_count_dict, total_count_info_dict


def set_mcifinder(args):

    mciscore_receptor = args.mciscore_receptor
    pf_receptor = args.pf_receptor
    mcinfo_ligand = args.mcinfo_ligand
    pf_receptor_info = args.pf_receptor_info
    mci_cutoff = args.mci_cutoff
    include_hydrophobic = args.include_hydrophobic
    dock_config = args.dock_config

    params = dict()
    weight_mcis = {
        'HBD:HBA': 0.129,
        'HBA:HBD': 0.129,
        'Anion:Cation': 0.31,
        'Cation:Anion': 0.31,
        'Cation:Aromatic': 0.039,
        'Hydrophobic:Hydrophobic': 0.00585,
        'Aromatic:Cation': 0.039,
        'Aromatic:Aromatic': 0.064,
        'MBA:Metal': 1.0
    }
    params['weight_mcis'] = weight_mcis
    params['bias_mcis'] = 4.0
    params['cutoff'] = mci_cutoff
    if dock_config is not None:
        params['dock_config'] = dock_config
    params['include_hydrophobic'] = include_hydrophobic
    params['mciscore_receptor'] = mciscore_receptor
    params['pf_receptor'] = pf_receptor
    params['mcinfo_ligand'] = mcinfo_ligand
    params['pf_receptor_info'] = pf_receptor_info

    pharma = Pharmacophore(params)

    return pharma


def main():

    import sys
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

    pharma = set_pifinder(args)
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
