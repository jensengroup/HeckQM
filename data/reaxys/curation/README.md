# FILE INFORMATION:

1) data_readifier.ipynb:

Run all cells in the notebook to convert the raw Reaxys data into a Pandas dataframe with curated data, which will be saved as "df_reaxys_data_curated.pkl". The notebook must be modified, if the names of the downloaded batches are different than the ones given above. 

2) rxnIDs.txt:

A list of unique Reaction IDs of the downloaded Reaxys data.
This content can be used to extract the data from Reaxys, 
via a single Reaxys query by following these instructions:

<ul>
    <li> Go to the “Query builder” tab in Reaxys.</li>
    <li> Use the “Search fields” functionality located in the top-right part of the window.</li>
    <li> Input “Reacton ID".</li>
    <li> Select "Reaction ID” from the filtered list. At the center of the page, the input field will appear.</li>
    <li> Paste the provided list of Reaction IDs separated by semicolons into the input field.</li>
    <li> Clicking on the “Reactions” button will redirect the user to the result page with all the reactions corresponding to Reaction IDs that have been provided.</li>
</ul>


# Data from the Reaxys repository

We have provided an importable .json file of the search query along with screen-dumbs of the query in the "queries" folder. 

For import of the .json file: 
1. Go to the “Query builder” tab on the Reaxys website (www.reaxys.com)
2. Select "Import" in the upper left corner, and choose the file "queries/query_2022_06_14_12_55_00.json".

Reaxys has a limit to the number of reactions that can be directly exported, so the data was exported as "Tab-delimited text" in several batches and saved inside the "raw_reaxys_data" folder.
Data export info:
1.  batch: 1 - 5,000 	    (raw_reaxys_data/data_1_5000.xls)
2.  batch: 5,001 - 10,000 	(raw_reaxys_data/data_5001_10000.xls)
3.  batch: 10,001 - 15,000 	(raw_reaxys_data/data_10001_15000.xls)
4.  batch: 15,001 - 20,000 	(raw_reaxys_data/data_15001_20000.xls)
5.  batch: 20,001 - 25,000	(raw_reaxys_data/data_20001_25000.xls)
6.  batch: 25,001 - 30,000	(raw_reaxys_data/data_25001_30000.xls)
7.  batch: 30,001 - 35,000	(raw_reaxys_data/data_30001_35000.xls)
8.  batch: 35,001 - 40,000	(raw_reaxys_data/data_35001_40000.xls)
9.  batch: 40,001 - 45,000	(raw_reaxys_data/data_40001_45000.xls)
10. batch: 25,001 - 30,000	(raw_reaxys_data/data_45001_50000.xls)
11. batch: 30,001 - 35,000	(raw_reaxys_data/data_50001_55000.xls)
12. batch: 35,001 - 40,000	(raw_reaxys_data/data_55001_60000.xls)
13. batch: 40,001 - 45,000	(raw_reaxys_data/data_60001_65000.xls)
14. batch: 40,001 - 45,000	(raw_reaxys_data/data_65001_66783.xls)

---------------------------------------------------------------------------------------------------

Another search for additional heck reactions resulting was carried out to find missing rxns from doi: 10.1002/anie.202109801.

A .json file with the query is provided in the "queries" folder and it is named "Selective Approaches to α- and β-Arylated Vinyl Ethers_query_2022_06_16_09_24_01.json".

The data was exported as "Tab-delimited text" and saved as "raw_reaxys_data/data_10.1002ange.202109801.xls".

---------------------------------------------------------------------------------------------------
