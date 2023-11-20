"""Reads data from a catalog file and stores it locally in a dictionary."""
import json
import numpy as np
from astroquery.vizier import Vizier
from config import COLUMNS, TEMP_FOLDER
from utilfuncs import find_nearest_galaxy

def retrieve_catalog(prompt:str, savepath = None, UI = 'Name', gal_title = 'dSph'):
    """Retrieves data from Vizier via one prompt"""
    V = Vizier(
        columns = COLUMNS,
        row_limit = -1
    )
    catalog_list = V.get_catalogs(prompt)
    catalog_keys = list(catalog_list.keys())

    for count, i in enumerate(catalog_keys):
        print(f'[{count}]', i)

    catalog_choice = int(input("Which catalog above would you like to use? Type index:"))

    catalog = catalog_list[catalog_keys[catalog_choice]]
    
    print(catalog.colnames)
    output = {}
    # catalog is a table object containing the data
    print("Adding stars...")

    for star in catalog:
        star_gal_name = star[gal_title]
        try:
            output[star_gal_name][star[UI]] = {}
        except KeyError:
            output[star_gal_name] = {star[UI]: {}}
            print("Adding newfound galaxy...")

        # Add catalog under star
        output[star_gal_name][star[UI]][prompt] = {}

        for col in catalog.colnames:
            if col == "RAJ2000" or col == "DEJ2000":
                star_col = str(star[col]).split(" ")
                star_cord = float(star_col[0])
                star_cord += float(star_col[1])/60
                star_cord += float(star_col[2])/3600
                
                star_col = star_cord
            
            elif col == gal_title:
                star_col = star[col]
            elif col == UI:
                continue
            elif star[col] == "--":
                star_col = 'NaN'
            else:
                star_col = float(star[col])

            output[star_gal_name][star[UI]][prompt][col] = star_col

    print("Done adding stars.")
    if savepath is not None:
        print("Saving catalog to JSON...")
        with open(savepath, 'w', encoding="utf-8") as f:
            json.dump(output, f, indent = 4)
    
    print("Done saving catalog.")
    return output

# Write a function to sort new stars into galaxies
def add_catalog(new_prompt, cache_path, catalog_name):
    """Adds a new catalog from Vizier to an existing cache."""
    # Load JSON file
    V = Vizier(
        columns = COLUMNS,
        row_limit = -1
    )
    catalog_list = V.get_catalogs(new_prompt)
    catalog_keys = list(catalog_list.keys())

    for count, i in enumerate(catalog_keys):
        print(f'[{count}]', i)

    catalog_choice = int(input("Which catalog above would you like to use? Type index:"))

    catalog = catalog_list[catalog_keys[catalog_choice]]

    # Load the cache and find the galaxy radii
    

# Execution code
if __name__ == "__main__":
    Kirby_2009 = retrieve_catalog("J/ApJS/191/352/abun",
                                  f"{TEMP_FOLDER}/Kirby_2009.json")

    test = retrieve_catalog("J/A+A/641/A127")

    # print(Kirby_2009)
    # add_catalog(test, f"{TEMP_FOLDER}/Kirby_2009.json", "J/ApJS/191/352/abun")
