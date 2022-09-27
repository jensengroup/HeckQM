# MIT License
#
# Copyright (c) 2022 Nicolai Ree
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import copy
import numpy as np
from operator import itemgetter
from concurrent.futures import ThreadPoolExecutor

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

import heckqm.molecule_formats as molfmt
import heckqm.run_xTB as run_xTB
import heckqm.run_orca as run_orca

# CPU usage
# -- Note, that ORCA is set to use 4 cpu cores and 2 conformers are running in parallel
#    resulting in a total of 8 cpu cores.
num_cpu_parallel = 2 # number of parallel jobs.
num_cpu_single = 4 # number of cpus per job


def confsearch_xTB(conf_complex_mols, template, conf_names, chrg=0, spin=0, method='ff', solvent='', conf_cutoff=50, contrained=True, precalc_path=None):
    
    global num_cpu_single
    
    # Run a constrained xTB optimizations  
    confsearch_args = []
    for i in range(len(conf_names)):
        if contrained:
            if precalc_path:
                confsearch_args.append((conf_complex_mols[i], template, conf_names[i]+'_contrained_gfn'+method.replace(' ', ''), chrg, spin, method, solvent, precalc_path[i]))
            else:
                confsearch_args.append((conf_complex_mols[i], template, conf_names[i]+'_contrained_gfn'+method.replace(' ', ''), chrg, spin, method, solvent, None))
        else:
            if precalc_path:
                confsearch_args.append((conf_names[i]+'_full_opt.xyz', conf_complex_mols[i], chrg, spin, method, solvent, True, precalc_path[i]))
            else:
                confsearch_args.append((conf_names[i]+'_full_opt.xyz', conf_complex_mols[i], chrg, spin, method, solvent, True, None))

    with ThreadPoolExecutor(max_workers=num_cpu_single) as executor:
        if contrained:
            results = executor.map(run_xTB.run_contrained_xTB, confsearch_args)
        else:
            results = executor.map(run_xTB.run_xTB, confsearch_args)

    conf_energies = []
    conf_paths = []
    for result in results:
        conf_energy, path_contrained_opt = result
        conf_energies.append(conf_energy)
        conf_paths.append(path_contrained_opt)

    # Find the conformers below cutoff #kJ/mol (12.6 kJ/mol = 3 kcal/mol)
    rel_conf_energies = np.array(conf_energies) - np.min(conf_energies) #covert to relative energies
    below_cutoff = (rel_conf_energies <= conf_cutoff).sum() #get number of conf below cutoff

    conf_tuble = list(zip(conf_names, conf_complex_mols, conf_paths, conf_energies, rel_conf_energies)) #make a tuble
    conf_tuble = sorted(conf_tuble, key=itemgetter(4))[:below_cutoff] #get only the best conf below cutoff

    conf_names, conf_complex_mols, conf_paths, conf_energies, rel_conf_energies = zip(*conf_tuble) #unzip tuble
    conf_names, conf_complex_mols, conf_paths, conf_energies = list(conf_names), list(conf_complex_mols), list(conf_paths), list(conf_energies) #tubles to lists
    mol_files = [os.path.join(item, item.split('/')[-1] + '_opt.sdf') for item in conf_paths] #list of paths to optimized structures in .sdf format

    # Find only unique conformers
    conf_names, conf_complex_mols, conf_paths, conf_energies = zip(*molfmt.find_unique_confs(list(zip(conf_names, conf_complex_mols, conf_paths, conf_energies)), mol_files, threshold=0.5)) #find unique conformers
    conf_names, conf_complex_mols, conf_paths, conf_energies = list(conf_names), list(conf_complex_mols), list(conf_paths), list(conf_energies) #tubles to lists

    return conf_names, conf_complex_mols, conf_paths, conf_energies


def calculateEnergy(args):
    """ Embed the post-insertion complex and calculate the ground-state free energy 
    return: energy [kJ/mol]
    """
    
    global num_cpu_single

    complex_mol, template, name, chrg = args
    method=' 2'  # <--- change the method for accurate calculations ('ff', ' 0', ' 1', ' 2')
    solvent = '--alpb Phenol' # <--- change the solvent ('--gbsa solvent_name', '--alpb solvent_name', or '')

    # Obtain product by deattaching it from the catalyst
    product = molfmt.run_rxn((complex_mol,), '[13#6:4][#6:1]([#1:3])[#6:2][#46:10]>>[13#6:4][#6:1]=[#6:2].[#46:10]')[0]
    chrg = chrg + Chem.GetFormalCharge(product) # checking and adding formal charge to "chrg"
    # OBS! Spin is hardcoded to zero!

    # Get number of conformers
    rot_bond = Chem.rdMolDescriptors.CalcNumRotatableBonds(product)
    n_conformers = min(1 + 3 * rot_bond, 80)

    # Embed the post-insertion complex
    selenium_indicies = complex_mol.GetSubstructMatches(Chem.MolFromSmarts('[Se:1]([*:2])[*:3]'))
    if selenium_indicies:
        complex_mol = molfmt.replaceSelenium(complex_mol, selenium_indicies, replace_Se=True) # replace Selenium with Sulfur
        Chem.SanitizeMol(complex_mol)
        complex_mol_conformers = molfmt.ConstrainedEmbedMultipleConfs(complex_mol, template, numConfs=n_conformers, useTethers=True, numThreads=num_cpu_single)
        complex_mol_conformers = molfmt.replaceSelenium(complex_mol_conformers, selenium_indicies, replace_Se=False) # reintroduce Selenium
    else:
        complex_mol_conformers = molfmt.ConstrainedEmbedMultipleConfs(complex_mol, template, numConfs=n_conformers, useTethers=True, numThreads=num_cpu_single)

    # Unpack confomers and assign conformer names
    conf_complex_mols = [Chem.Mol(complex_mol_conformers, False, i) for i in range(complex_mol_conformers.GetNumConformers())]
    conf_names = [name + f'_conf{str(i+1).zfill(2)}' for i in range(complex_mol_conformers.GetNumConformers())] #change zfill(2) if more than 99 conformers
    conf_names_copy = copy.deepcopy(conf_names)

    # Run a constrained GFN-FF optimization
    conf_names, conf_complex_mols, conf_paths, conf_energies = confsearch_xTB(conf_complex_mols, template, conf_names, chrg=chrg, spin=0, method='ff', solvent=solvent, conf_cutoff=10, contrained=True, precalc_path=None)
    
    # Run a constrained GFN?-xTB optimization
    conf_names, conf_complex_mols, conf_paths, conf_energies = confsearch_xTB(conf_complex_mols, template, conf_names, chrg=chrg, spin=0, method=method, solvent=solvent, conf_cutoff=10, contrained=True, precalc_path=conf_paths)
    
    # Run a full GFN?-xTB optimization
    conf_names, conf_complex_mols, conf_paths, conf_energies = confsearch_xTB(conf_complex_mols, template, conf_names, chrg=chrg, spin=0, method=method, solvent=solvent, conf_cutoff=1e100, contrained=False, precalc_path=conf_paths)

    # Run Orca single point calculations
    final_conf_energies = []
    final_conf_complex_mols = []
    for conf_name, conf_complex_mol, conf_path, conf_energy in zip(conf_names, conf_complex_mols, conf_paths, conf_energies):
        if conf_energy != 60000.0: # do single point calculations on all unique conformers
            conf_energy = run_orca.run_orca('xtbopt.xyz', chrg, os.path.join("/".join(conf_path.split("/")[:-2]), 'full_opt', conf_name+'_full_opt'), ncores=num_cpu_single, mem=10000)
        final_conf_energies.append(conf_energy)
        final_conf_complex_mols.append(conf_complex_mol)
    
    # Get only the lowest energy conformer
    minE_index = np.argmin(final_conf_energies)
    best_conf_complex_mol = final_conf_complex_mols[minE_index]
    best_conf_energy = final_conf_energies[minE_index] # uncomment when doing single point calculations on all unique conformers
    # best_conf_energy = run_orca.run_orca('xtbopt.xyz', chrg, os.path.join(os.getcwd(), 'calc', conf_names[minE_index]+'_full_opt'), ncores=num_cpu_single, mem=20000) # comment when doing single point calculations on all unique conformers, otherwise this runs a Orca single point calculation on the lowest xTB energy conformer

    # Remove isotope infomation
    best_conf_complex_mol = MolStandardize.rdMolStandardize.IsotopeParent(best_conf_complex_mol, skipStandardize=True)
    product = MolStandardize.rdMolStandardize.IsotopeParent(product, skipStandardize=True)

    ### START - CLEAN UP ###
    for conf_name in conf_names_copy:

        conf_path = os.path.join(os.getcwd(), 'calc', conf_name.split('_')[0], conf_name.split('_')[1])
        
        if os.path.isfile(os.path.join(conf_path, conf_name+'_full_opt.xyz')):
            os.remove(os.path.join(conf_path, conf_name+'_full_opt.xyz'))
        
        # Remove GFNFF-xTB folder
        folder_path = os.path.join(conf_path, 'contrained_gfnff', conf_name + '_contrained_gfnff')
        if os.path.exists(folder_path):
            for file_remove in os.listdir(folder_path):
                if os.path.isfile(f'{folder_path}/{file_remove}'):
                    os.remove(f'{folder_path}/{file_remove}')
            # checking whether the folder is empty or not
            if len(os.listdir(folder_path)) == 0:
                os.rmdir(folder_path)
            else:
                print("Folder is not empty")
    
        # Remove GFN?-xTB folder
        folder_path = os.path.join(conf_path, 'contrained_gfn' + method.replace(' ', ''), conf_name + '_contrained_gfn' + method.replace(' ', ''))
        if os.path.exists(folder_path):
            for file_remove in os.listdir(folder_path):
                if os.path.isfile(f'{folder_path}/{file_remove}'):
                    os.remove(f'{folder_path}/{file_remove}')
            # checking whether the folder is empty or not
            if len(os.listdir(folder_path)) == 0:
                os.rmdir(folder_path)
            else:
                print("Folder is not empty")

        # Clean full opt folder
        folder_path = os.path.join(conf_path, 'full_opt', conf_name + '_full_opt')
        file_remove_list = ['charges', 'coordprot.0', 'lmocent.coord', 'orca_calc_atom46.densities',
                        'orca_calc_atom46.out', 'orca_calc_atom46_property.txt', 'orca_calc_atom53.densities',
                        'orca_calc_atom53.out', 'orca_calc_atom53_property.txt', 'orca_calc.cpcm',
                        'orca_calc.densities', 'orca_calc.gbw', 'wbo', 'xtblmoinfo', 'xtbopt.log',
                        '.xtboptok', 'xtbrestart', 'xtbscreen.xyz']
        if os.path.exists(folder_path):
            for file_remove in file_remove_list:
                if os.path.isfile(f'{folder_path}/{file_remove}'):
                    os.remove(f'{folder_path}/{file_remove}')
    
    folder_path = os.path.join(conf_path, 'contrained_gfnff')
    if os.path.exists(folder_path):
        os.rmdir(folder_path)
    
    folder_path = os.path.join(conf_path, 'contrained_gfn'+method.replace(' ', '')) 
    if os.path.exists(folder_path):
        os.rmdir(folder_path)
    ### END - CLEAN UP ###

    return best_conf_energy, best_conf_complex_mol, product


def run_heck_reaction(alkene_mol, halogen_mol, name='mol', chrg=0, neutral_path=True):
    """ Creates post-insertion complexes and calculates the energies of these complexes 
    in order to predict the most probable product of a Heck reaction
    return: energies [kJ/mol] (post-insertion complex energies), complexes (rdkit Mol objects), products (rdkit Mol objects), legends (relative energies)
    """
    
    global num_cpu_parallel

    ### Select the correct templates and reaction SMARTS depending on chosen the reaction path ###
    if neutral_path:
        # NEUTRAL REACTION TEMPLATES
        template_alpha = Chem.SDMolSupplier('heckqm/template_cat_iso/first_temp_cat_31.sdf', removeHs=False, sanitize=True)[0] # alpha / branched (alkene mol associate trans to the phosphine)
        template_beta = Chem.SDMolSupplier('heckqm/template_cat_iso/first_temp_cat_29.sdf', removeHs=False, sanitize=True)[0] # beta / linear (alkene mol associate trans to the phosphine)
        template_alpha_th = Chem.SDMolSupplier('heckqm/template_cat_iso/first_temp_cat_37.sdf', removeHs=False, sanitize=True)[0] # alpha / branched (alkene mol associate trans to the halide)
        template_beta_th = Chem.SDMolSupplier('heckqm/template_cat_iso/first_temp_cat_35.sdf', removeHs=False, sanitize=True)[0] # beta / linear (alkene mol associate trans to the halide)
        templates = [template_alpha, template_beta, template_alpha_th, template_beta_th]
        template_names = ['4', '3', '2', '1'] #['31', '29', '37', '35']

        # NEUTRAL REACTION SMARTS
        if halogen_mol:
            inter_one_alpha_complex = '[*;!#1:3][#6:1]([#1:4])=&!@[#6:2]([#1:5])[#1:6].[#6;$([#6]=[#6]),$(c:c):7][Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:8].[Cl:9][#46:10][#6:11][#6:12]([#1:13])([#6:14])[#6:15]>>[*:3][#6:1]([#1:4])([13#6:7])[#6:2]([#1:5])([#46:10][*:8])[#1:6].[Cl:9]'
            inter_one_beta_complex = '[*;!#1:3][#6:1]([#1:4])=&!@[#6:2]([#1:5])[#1:6].[#6;$([#6]=[#6]),$(c:c):7][Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:8].[Cl:9][#46:10][#6:11][#6:12]([#1:13])[#6:14]>>[*:3][#6:1]([#1:4])([#46:10][*:8])[#6:2]([#1:5])[13#6:7].[Cl:9]'
            smarts_list = [inter_one_alpha_complex, inter_one_beta_complex, inter_one_alpha_complex, inter_one_beta_complex]
        else:
            intra_one_alpha_complex = '([*;!#1:3][#6:1]([#1:4])=&!@[#6:2]([#1:5])[#1:6].[#6;$([#6]=[#6]),$(c:c):7][Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:8]).[Cl:9][#46:10][#6:11][#6:12]([#1:13])([#6:14])[#6:15]>>[*:3][#6:1]([#1:4])([13#6:7])[#6:2]([#1:5])([#46:10][*:8])[#1:6].[Cl:9]'
            intra_one_beta_complex = '([*;!#1:3][#6:1]([#1:4])=&!@[#6:2]([#1:5])[#1:6].[#6;$([#6]=[#6]),$(c:c):7][Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:8]).[Cl:9][#46:10][#6:11][#6:12]([#1:13])[#6:14]>>[*:3][#6:1]([#1:4])([#46:10][*:8])[#6:2]([#1:5])[13#6:7].[Cl:9]'
            smarts_list = [intra_one_alpha_complex, intra_one_beta_complex, intra_one_alpha_complex, intra_one_beta_complex]
    else:
        # CATIONIC REACTION TEMPLATES
        template_alpha = Chem.SDMolSupplier('heckqm/template_cat_iso/first_temp_cat_17.sdf', removeHs=False, sanitize=True)[0] # alpha / branched
        template_beta = Chem.SDMolSupplier('heckqm/template_cat_iso/first_temp_cat_09.sdf', removeHs=False, sanitize=True)[0] # beta / linear
        templates = [template_alpha, template_beta]
        template_names = ['5', '6'] #['17', '09']

        # CATIONIC REACTION SMARTS
        if halogen_mol:
            inter_one_alpha_complex = '[*;!#1:3][#6:1]([#1:4])=&!@[#6:2]([#1:5])[#1:6].[#6;$([#6]=[#6]),$(c:c):7][Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:8].[#46:9][#6:10][#6:11]([#1:12])([#6:13])[#6:14]>>[*:3][#6:1]([#1:4])([13#6:7])[#6:2]([#1:5])([#46:9])[#1:6]'
            inter_one_beta_complex = '[*;!#1:3][#6:1]([#1:4])=&!@[#6:2]([#1:5])[#1:6].[#6;$([#6]=[#6]),$(c:c):7][Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:8].[#46:9][#6:10][#6:11]([#1:12])[#6:13]>>[*:3][#6:1]([#1:4])([#46:9])[#6:2]([#1:5])[13#6:7]'
            smarts_list = [inter_one_alpha_complex, inter_one_beta_complex]
        else:
            intra_one_alpha_complex = '([*;!#1:3][#6:1]([#1:4])=&!@[#6:2]([#1:5])[#1:6].[#6;$([#6]=[#6]),$(c:c):7][Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:8]).[#46:9][#6:10][#6:11]([#1:12])([#6:13])[#6:14]>>[*:3][#6:1]([#1:4])([13#6:7])[#6:2]([#1:5])([#46:9])[#1:6]'
            intra_one_beta_complex = '([*;!#1:3][#6:1]([#1:4])=&!@[#6:2]([#1:5])[#1:6].[#6;$([#6]=[#6]),$(c:c):7][Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:8]).[#46:9][#6:10][#6:11]([#1:12])[#6:13]>>[*:3][#6:1]([#1:4])([#46:9])[#6:2]([#1:5])[13#6:7]'
            smarts_list = [intra_one_alpha_complex, intra_one_beta_complex]


    ### Generate post-insertion complexes ###
    complex_mols_list = []
    modified_templates = []
    names = []
    for smarts, template, template_name in zip(smarts_list, templates, template_names):
        
        # Obtain position of the dummy atom and its neighbor, and replace the dummy atom with a carbon atom 
        templateIndex1, templateIndex2 = molfmt.getAttachmentVector(template)
        modified_template = molfmt.replaceAtom(template, templateIndex1, templateIndex2, atom_type='C')

        # Run rxn to obtain a SMILES of the post-insertion complex
        if halogen_mol:
            reactants = (alkene_mol, halogen_mol, modified_template)
        else:
            reactants = (alkene_mol, modified_template)
        complex_mols = molfmt.run_rxn(reactants, smarts)
        complex_mols = [Chem.AddHs(complex_mol) for complex_mol in complex_mols]

        # Append to lists
        complex_mols_list.extend(complex_mols)
        names.extend([f'{name}_{template_name}_{i}' for i in range(len(complex_mols))])

        for complex_mol in complex_mols:

            modified_template_copy = copy.deepcopy(modified_template)

            # Change the connected alkene R-group atom in the template if not carbon
            alkene_atom_index = complex_mol.GetSubstructMatch(Chem.MolFromSmarts(smarts.split('>>')[-1].split('.')[0]))[0]
            alkene_atom_type = complex_mol.GetAtomWithIdx(alkene_atom_index).GetSymbol()
            if alkene_atom_type != 'C':
                modified_template_copy = molfmt.replaceAtom(modified_template_copy, templateIndex1, templateIndex2, atom_type=alkene_atom_type)

            # Change the halide atom in the template if not Chlorine (in case of neutral pathway)
            if neutral_path:
                halide_atom_index = complex_mol.GetSubstructMatch(Chem.MolFromSmarts('[Cl,Br,$(OS(=O)(=O)C(F)(F)F),I:1][#46]'))[0] # extract the reacting halide of halogen_mol (in case of neutral pathway)
                halide_atom_type = complex_mol.GetAtomWithIdx(halide_atom_index).GetSymbol()
                
                if halide_atom_type != 'Cl':
                    indexAtom, indexNeighbor = modified_template_copy.GetSubstructMatch(Chem.MolFromSmarts('[Cl][#46]'))
                    modified_template_copy = molfmt.replaceAtom(modified_template_copy, indexAtom, indexNeighbor, atom_type=halide_atom_type)
                
            #     modified_template_copy = AllChem.ReactionFromSmarts('[#6:2][#46:10]([P:11]([#6:12])([#6:13])[#6:14])[Cl,Br,O,I:15]>>[#6:2][#46:10]([P:11])[*:15].[*:12].[*:13].[*:14]').RunReactants((modified_template_copy,))[0][0] # limit constraining core
            # else:
            #     modified_template_copy = AllChem.ReactionFromSmarts('[#6:2][#46:10]([P:11]([#6:12])([#6:13])[#6:14])[P:15]([#6:16])([#6:17])[#6:18]>>[#6:2][#46:10]([P:11])[P:15].[*:12].[*:13].[*:14].[*:16].[*:17].[*:18]').RunReactants((modified_template_copy,))[0][0] # limit constraining core

            modified_templates.append(modified_template_copy)


    ### Calculate energies for all of the post-insertion complexes ###
    args = [(complex_mol, modified_template, name, chrg) for complex_mol, modified_template, name in zip(complex_mols_list, modified_templates, names)]
    with ThreadPoolExecutor(max_workers=num_cpu_parallel) as executor:
        results = executor.map(calculateEnergy, args)

    energies = []
    complexes = []
    products = []
    for result in results:
        energy, complex_mol, product = result
        energies.append(energy)
        complexes.append(complex_mol)
        products.append(product)

    # Calculate the relative energies of all products --> product legends
    minE = min(energies)
    legends = [f'{energy - minE:.2f} kJ/mol [{name}]' for energy, name in zip(energies, template_names)]

    # Sort lists (energies, complexes, products, and legends) according to the energies in ascending order
    all_lists = list(zip(energies, complexes, products, legends))
    all_lists = sorted(all_lists, key=itemgetter(0))
    energies = [item[0] for item in all_lists]
    complexes = [item[1] for item in all_lists]
    products = [item[2] for item in all_lists]
    legends = [item[3] for item in all_lists]

    return energies, complexes, products, legends


def control_calcs(df):

    df['energies_neu'] = pd.Series(dtype='object')
    df['complexes_neu'] = pd.Series(dtype='object')
    df['products_neu'] = pd.Series(dtype='object')
    df['legends_neu'] = pd.Series(dtype='object')
    df['energies_cat'] = pd.Series(dtype='object')
    df['complexes_cat'] = pd.Series(dtype='object')
    df['products_cat'] = pd.Series(dtype='object')
    df['legends_cat'] = pd.Series(dtype='object')
    for idx, row in df.iterrows():

        rxn_id = row['Reaction ID']
        reactants = row['Reactant(s)']
        alkene_smi = reactants[0]
        alkene_mol = Chem.AddHs(Chem.MolFromSmiles(alkene_smi))
        if len(reactants) > 1:
            halogen_smi = reactants[1]
            halogen_mol = Chem.AddHs(Chem.MolFromSmiles(halogen_smi))
        else:
            halogen_mol = None
        
        ### RUN CALCULATIONS ###
        # Neutral reaction path
        try:
            energies_neu, complexes_neu, products_neu, legends_neu = run_heck_reaction(alkene_mol, halogen_mol, name=str(rxn_id), chrg=0, neutral_path=True)
            df.at[idx, 'energies_neu'] = energies_neu
            df.at[idx, 'complexes_neu'] = complexes_neu
            df.at[idx, 'products_neu'] = products_neu
            df.at[idx, 'legends_neu'] = legends_neu
        except Exception as e:
            print(f'WARNING! HeckQM Neutral failed for {rxn_id}')
            print(e)

        # Cationic reaction path
        try:
            energies_cat, complexes_cat, products_cat, legends_cat = run_heck_reaction(alkene_mol, halogen_mol, name=str(rxn_id), chrg=1, neutral_path=False)
            df.at[idx, 'energies_cat'] = energies_cat
            df.at[idx, 'complexes_cat'] = complexes_cat
            df.at[idx, 'products_cat'] = products_cat
            df.at[idx, 'legends_cat'] = legends_cat
        except Exception as e:
            print(f'WARNING! HeckQM Cationic failed for {rxn_id}')
            print(e)

    return df


if __name__ == "__main__":
    
    import pandas as pd
    import submitit

    df = pd.read_pickle('data_cabri/data_cabri.pkl') #Cabri dataset
    # df = pd.read_pickle('data_curation/df_all_data_curated.pkl') #all curated data

    executor = submitit.AutoExecutor(folder="submitit_heckqm")
    executor.update_parameters(
        name="heckQM",
        cpus_per_task=int(num_cpu_parallel*num_cpu_single),
        mem_gb=20, # remember to change this for the ORCA calculations!
        timeout_min=6000,
        slurm_partition="kemi1",
        slurm_array_parallelism=50,
    )
    print(executor)

    jobs = []
    with executor.batch():
        chunk_size = 1 
        for start in range(0, df.shape[0], chunk_size):
            df_subset = df.iloc[start:start + chunk_size]
            job = executor.submit(control_calcs, df_subset)
            jobs.append(job)
