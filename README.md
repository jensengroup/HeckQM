This folder contains supporting material for the HeckSQM paper such as code for dataset curation, regioselectivity predictions using automated QM, and machine learning training/testing.

Unfortunately, the dataset cannot be shared directly due to licensing restrictions, but it can easily  be extracted from Reaxys via a single Reaxys query and a list of Reaction IDs separated by semicolons.
In order to perform such a query, a user should:
1. Go to the “Query builder” tab in Reaxys.
2. Use the “Search fields” functionality located in the top-right part of the window.
3. Input “Reacton ID".
4. Select "Reaction ID” from the filtered list. At the center of the page, the input field will appear.
5. Paste the provided list of Reaction IDs into the input field.
6. Clicking on the “Reactions” button will redirect the user to the result page with all the reactions corresponding to Reaction IDs that have been provided.

Otherwise follow the procedure provided in "data_curation/info.txt".

-------------------------------------------------------------------------------------------

### INSTALLATION ###

We recommend using anaconda to install the Python 3 environment:
conda env create -f environment.yml && conda activate hecksqm

Then download the binaries of xtb version 6.4.1:
mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.4.1/xtb-6.4.1.tar.xz; tar -xvf ./xtb-6.4.1.tar.xz; cd ..

Furthermore, ORCA version 5.0.1 must be installed following the instructions found here: https://sites.google.com/site/orcainputlibrary/setting-up-orca

OBS! The path to ORCA must be modified in "hecksqm/run_orca.py".
