# Chemical Abundances of Metal-Poor Stars (CAMPS) Catalogue

## Overview
The chemical abundances of stars are valuable data for many fields of astrophysics. In particular, the chemical abundances of metal-poor stars are especially unique because these stars tend to be quite old and are from very different environments. However, most largely available chemical abundance data are for more metal-rich, sun-like stars. A large fraction of metal-poor star data is found in research papers which can be conveniently compiled into literature catalogues. Unfortunately, current literature compilations are not able to be updated are gradually becoming outdated. The Chemical Abundances of Metal-Poor Stars (CAMPS) catalogue builds upon preexisting literature compilations by allowing for papers to be added overtime without compromising the running naming conventions. CAMPS is designed to retrieve data from the Vizier database, merge tables from the same catalog, and compile them into a single catalog file along with other preexisting literature compilations like SAGA and JINAbase. 
CAMPS contains abundance data for ~13,000 unique stars with a mean [Fe/H] ~ -1.6 dex.

## Installation

The tools to compile CAMPS are available in this repository. The data for SAGA and JINAbase are not available in this repository but are available for download from their respective sources, located below.

### Manual Compiling
 To use the tool, download the repository and open it in your Python environment (I opened it in a JupyterHub workspace), then open the `cat_builder.ipynb` file to build a catalogue. The tool requires the following packages:

- Numpy
- Scipy
- Matplotlib
- Astropy
- Pandas
- json
- astroquery
- collections

### Compiling External Sources
CAMPS has the option to add on preexisting literature compilations, SAGA, JINAbase, and Dr. Alexander Ji's dwarf galaxy catalogue. Dr. Ji's catalogue is located in the repository as `dwarf_lit_all_coo.org`, SAGA and JINAbase are not located in the repository, they can be downloaded here:

SAGA: http://sagadatabase.jp/

JINAbase: https://jinabase.pythonanywhere.com/ 

Dr. Ji's dwarf galaxy catalogue was downloaded from https://github.com/alexji/alexmods/tree/master/alexmods/data 

SAGA can be added with the `Add_SAGA` boolean set to True. Path names may need to be adjusted in the `addsurveys.py` file, under the `saga` function
JINAbase can be added with the `Add_JINA` boolean set to True. Path names may need to be adjusted in the `addsurveys.py` file, under the `saga` function
Dr. Ji's dwarf galaxy catalogue can be added with the `Add_Ji` boolean set to True


### Including Gaia Identification
For included stars to be given their Gaia ID, if it exists, the `Add_Gaia` boolean should be set to True. 
Gaia cross-matching was performed with the Whole Sky Database (WSDB), constructed by Dr. Sergey Koposov. This database and its functions are not public, so `Add_Gaia` should be set to false if you do not have access. 

## Directory
- CAMPS_Code: folder for all the filed required to compile CAMPS
    - config.json: configuration file. This file keeps track of the folder paths and the column headers of the data files. When importing a new catalog, the tool checks this file to see if it has already seen the column headers. If it has, it will automatically rename the column headers to be consistent with the cache file. If it hasn't, it will prompt you to rename the column headers. This file is essential for the tool to function properly.
    - dataloader.py: the main script for loading data. This script contains the functions to retrieve data.
    - utilfuncs.py: utility functions
    - addsurveys.py: the main script for combining data from various sources. This script contains the functions to combine, clean, and format the data.
    - CAMPS_Metadata_Dict.csv: Contains information on all of the papers in CAMPS; their shortened name, bibcode, where the data was compiled from, the solar scale used in the paper, and the instrument and its resolution.
    - dwarf_lit_all_coo.org: Dr. Ji's dwarf galaxy table
- Original_CACT: The original paper compiler made by Shuhan Zheng

## How to use
1. Fill in booleans and the `paper_list` dictionary:
- Set booleans to true or false depending on which catalogues you would like to add
- Fill in the paper list dictionary. For each paper specify which tables from Vizier you want, what the Vizier tag is, temporary storage path, shortened paper name, name of the column containing the stars, name of the column containing the galaxies (if no such column exists, write None and specify a galaxy name), and the paper type 'Vizier'
- If you would like to add extra catalogues stored as an Astropy table locally, include them in the custom dictionary
2. Run cells below:
- If an incoming paper has column names that are not recognized, you will be prompted to input a column name. Input _ if you do not wish for this column to be included
3. Run metadata table and/or collapse table:
- Metadata table generates one like CAMPS_Metadata_Dict.csv
- Collapsed table will remove abundance species differentiation by keeping the first instance of any measurement of an abundance.

## Known Problems
- CAMPS cannot handle importing papers that aren't in Vizier or added separately in an Astropy format.
- Incoming tables must have each row corresponding to a star. They cannot have stars along the columns.
- Data is not standardized to the same solar scale, however the scale is specified for each star.
- Information on a stars solar scale, facility, and resolution are not automatically added and are done based on what is in the metadata file.

## Included Catalogs & Referencing
All included papers' bibcodes are in the CAMPS_Metadata_Dict.csv. The bibtex file for all relevant papers is located in the file `CAMPS_citations.bibtex`

Any paper used should be cited directly.

If the paper comes from SAGA, follow the SAGA citation requirements specified here: http://sagadatabase.jp/

If the paper comes from JINAbase, cite the following paper: https://iopscience.iop.org/article/10.3847/1538-4365/aadfe9


This work was co-created with Shuhan Zheng.
