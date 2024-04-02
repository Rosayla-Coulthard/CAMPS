# Catalog Acquisition and Compilation Tool (CACT)

## How to use
1. Retreive data from Vizier:
Run `retrieve_catalog()` following the instructions in the docstring and the prompts in the terminal. Save the output as a separate .json file. If `config.json` has no record of a column header contained in the retrieved table, it will prompt you to rename the column as it will appear in the cache file. If the column header is already in `config.json`, it will automatically rename the column as it appears in the cache file. Once all new column headers have been renamed and added, the script will terminate and you will need to run it again. 
2. Merge tables:
If the Vizier prompt retrieves multiple tables, save them as separate `.json` files and run `merge_tables()` to merge them into one table.
3. Create catalog compilation:
Copy one of the existing `.json` catalog files and make it the cache file. Run `add_catalog()` to add one of the already retrieved catalogs to the cache file.

## Included catalogs
- Kirby 2009 (J/ApJS/191/352):
    - Title: Multi-element Abundance Measurements from Medium-resolution Spectra. II. Catalog of Stars in Milky Way Dwarf Satellite Galaxies
    - Members: 2947 stars in 8 galaxies.
    - Telescope: Keck
    - Instruments: DEIMOS (R ~ 6,000)
    - Fitting method: synthetic spectra, equivalent widths

- Kirby 2017 (J/ApJ/838/83):
    - Title: Triangulum II. Not Especially Dense After All
    - Members: 13 stars in 1 galaxies.
    - Telescope: Keck
    - Instruments: DEIMOS (R ~ 6,000)
    - Fitting method: equivalent widths

- Kirby 2015 (J/ApJ/801/125):
    - Title: Carbon in Red Giants in Globular Clusters and Dwarf Spheroidal Galaxies
    - Members: 745 stars in 7 galaxies.
    - Telescope: Keck
    - Instruments: DEIMOS (R ~ 6,000)
    - Fitting method: synthetic spectra, equivalent widths

- Kirby 2017-A (J/ApJ/834/9):
    - Title: Chemistry and Kinematics of the Late-forming Dwarf Irregular Galaxies Leo A, Aquarius, and Sagittarius DIG
    - Members: 197 stars in 3 galaxies.
    - Telescope: Keck
    - Instruments: DEIMOS (R ~ 6,000)
    - Fitting method: synthetic spectra, equivalent widths

- Duggan 2018 (J/ApJ/869/50):
    - Title: Neutron Star Mergers are the Dominant Source of the r-process in the Early Evolution of Dwarf Galaxies
    - Members: 243 stars in 5 galaxies.
    - Telescope: Keck
    - Instruments: DEIMOS (R ~ 6,000)
    - Fitting method: synthetic spectra

- de los Reyes 2020 (J/ApJ/891/85):
    - Title: Manganese Indicates a Transition from Sub- to Near-Chandrasekhar Type Ia Supernovae in Dwarf Galaxies
    - Members: 241 stars in 9 galaxies.
    - Telescope: Keck
    - Instruments: DEIMOS (R ~ 6,000)
    - Fitting method: synthetic spectra


- Theler 2020 (J/A+A/642/A176):
    - Title: The chemical evolution of the dwarf spheroidal galaxy Sextans
    - Members: 87 stars in 1 galaxies.
    - Telescope: VLT
    - Instruments: FLAMES/GIRAFFE (R ~ 20,000), UVES (R ~ 47,000)
    - Fitting method: equivalent widths, synthetic spectra

- Skúladóttir 2019 (J/A+A/631/A171):
    - Title: Neutron-capture elements in dwarf galaxies
    - Members: 98 stars in 1 galaxies.
    - Telescope: VLT
    - Instruments: FLAMES/GIRAFFE (R ~ 19,300 ~ 22,500), UVES (R ~ 47,000)
    - Fitting method: synthetic spectra

- Reichert 2020 (J/A+A/641/A127):
    - Title: Neutron-capture elements in dwarf galaxies: A homogenized analysis of 13 dwarf spheroidal and ultra-faint galaxies
    - Members: 380 stars in 13 galaxies.
    - Telescope: ESO archive and Keck archive
    - Instruments: FLAMES/GIRAFFE (295 observations), UVES (56 observations), X-shoooter (2 observations), HIRES (27 observations)
    - Fitting method: equivalent widths, synthetic spectra

- Kirby 2020 (J/AJ/159/46):
    - Title: Elemental Abundances in M31: The Kinematics and Chemical Evolution of Dwarf Spheroidal Satellite Galaxies
    - Members: 256 stars in 5 galaxies.
    - Telescope: Keck
    - Instruments: DEIMOS (R ~ 6,000)
    - Fitting method: equivalent widths, synthetic spectra

