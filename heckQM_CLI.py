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
import numpy as np
import submitit
import argparse

from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from heckQM import run_heck_reaction

# CPU usage
# -- change this in heckQM.py. Note, that ORCA is set to use 4 cpu cores and 2 conformers are running in parallel
#    resulting in a total of 8 cpu cores.



def parse_args():
    """
    Argument parser so this can be run from the command line
    """
    parser = argparse.ArgumentParser(description='Run regioselectivity predictions from the command line')
    parser.add_argument('-sa', '--alkene_smi', default='[H]C([H])=C([H])c1c([H])c([H])c([H])c([H])c1[H]',
                        help='SMILES input for alkene containing molecule')
    parser.add_argument('-sh', '--halogen_smi', default='[H]c1c([H])c([H])c(Cl)c([H])c1[H]',
                        help='SMILES input for halogen containing molecule')
    parser.add_argument('-n', '--name', default='testRXN', help='The name of the reaction (one word without "_")')
    return parser.parse_args()


def html_output(svg_neu, svg_cat):

    html = """<!DOCTYPE html>
            <html lang="en">
            <head>
                <title>HeckQM</title>

                <meta charset="utf-8">
                <meta name="google" content="notranslate" />
                <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no" />

            </head>
            <body>
                <h1 style="color:black; font-family:verdana; text-align:center;">
                    HeckQM - Predicted Product(s)
                </h1>

                <h2 style="color:#e41a1c; font-family:verdana; text-align:center;">
                    NEUTRAL PATHWAY
                </h2>
                <div style="text-align:center;">
                    {0}
                </div>
                
                <hr size="2" noshade>
                
                <h2 style="color:#e41a1c; font-family:verdana; text-align:center;">
                    CATIONIC PATHWAY
                </h2>
                <div style="text-align:center;">
                    {1}
                </div>
            </body>
            </html>""".format(svg_neu, svg_cat)
            
    return html



if __name__ == "__main__":

    args = parse_args()
    
    ### SLURM SETTINGS ###
    executor = submitit.AutoExecutor(folder="submitit_heckqm_cli")
    executor.update_parameters(
        name="heckQM",
        cpus_per_task=24,
        mem_gb=40,
        timeout_min=6000,
        slurm_partition="kemi1",
        slurm_array_parallelism=50,
    )
    ### END ###


    ### LOAD REACTANTS ###
    # Load alkene containing molecule
    alkene_mol = Chem.AddHs(Chem.MolFromSmiles(args.alkene_smi))
    
    # Load halogen containing molecule
    if args.halogen_smi:
        halogen_mol = Chem.AddHs(Chem.MolFromSmiles(args.halogen_smi))
        cddd_embedding = cddd_server.smiles_to_cddd([Chem.MolToSmiles(alkene_mol), Chem.MolToSmiles(halogen_mol)])
    else:
        halogen_mol = None
        cddd_embedding = cddd_server.smiles_to_cddd([Chem.MolToSmiles(alkene_mol), Chem.MolToSmiles(alkene_mol)])
    ### END ###
    

    ### RUN CALCULATIONS ###
    jobs = []
    with executor.batch():
        # Neutral reaction path
        job = executor.submit(run_heck_reaction, alkene_mol, halogen_mol, name=args.name, chrg=0, neutral_path=True)
        jobs.append(job)

        # Cationic reaction path
        job = executor.submit(run_heck_reaction, alkene_mol, halogen_mol, name=args.name, chrg=1, neutral_path=False)
        jobs.append(job)
    ### END ###


    ### MAKE HECKQM OUTPUT ###
    # Read results
    energies_neu, complexes_neu, products_neu, legends_neu = jobs[0].result()
    energies_cat, complexes_cat, products_cat, legends_cat = jobs[1].result()
    
    # Generate graphical output 
    """ The generated html code can be displayed in a notebook with the following commands:
    >>> from IPython.core.display import HTML 
    >>> HTML(html) """
    svg_neu = MolsToGridImage(products_neu, legends=legends_neu, molsPerRow=4, subImgSize=(300,300), useSVG=True)#.data
    svg_cat = MolsToGridImage(products_cat, legends=legends_cat, molsPerRow=4, subImgSize=(300,300), useSVG=True)#.data

    html = html_output(svg_neu, svg_cat) # generate html code

    os.makedirs(os.path.join(os.getcwd(), 'results'), exist_ok=True) # save graphical output in the results folder
    with open(f'results/{args.name}.html', 'w') as f: 
        f.write(html)
    ### END ###
