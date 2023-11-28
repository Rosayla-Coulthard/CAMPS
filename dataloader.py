"""Reads data from a catalog file and stores it locally in a dictionary."""
import json
import numpy as np
from astroquery.vizier import Vizier
from datetime import datetime
from config import COLUMNS, TEMP_FOLDER, META_DATA_KEY
from utilfuncs import find_nearest_galaxy

def retrieve_catalog(prompt:str, savepath = None, UI = 'Name', gal_title = 'dSph'):
    """Retrieves data from Vizier via one prompt.
    
    Args:
        prompt: The prompt to search for in Vizier.
        savepath: The path to save the catalog to.
        UI: The column name for the unique identifier.
        gal_title: The column name for the galaxy title.
    """
    V = Vizier(
        columns = COLUMNS,
        row_limit = -1
    )
    catalog_list = V.get_catalogs(prompt)
    catalog_keys = list(catalog_list.keys())

    for count, i in enumerate(catalog_keys):
        print(f'[{count}]', i)

    input_msg = "Which catalog above would you like to use? Type index: (Type anything else to quit)"
    catalog_choice = input(input_msg)

    try:
        catalog_choice = int(catalog_choice) 
        catalog = catalog_list[catalog_keys[catalog_choice]]
    except ValueError:
        print(f"[{datetime.now()}] WARNING: Invalid input. Operation aborted, no changes made.")
        quit()

    print(f"[{datetime.now()}] Found the following columns: {catalog.colnames}")
    output = {f"{META_DATA_KEY}": {"Included titles": {prompt: catalog.colnames}}}
    # catalog is a table object containing the data
    print(f"[{datetime.now()}] Adding stars...")

    for star in catalog:
        try:
            star_gal_name = star[gal_title]
        except KeyError:
            err_msg = "The catalog does not have a galaxy title column. "
            err_msg += "Please use a different catalog that has one."
            print(err_msg)
            print(f"[{datetime.now()}] WARNING: Operation aborted, no changes made.")
            quit()

        try:
            output[star_gal_name][star[UI]] = {}
        except KeyError:
            output[star_gal_name] = {star[UI]: {}}
            print(f"[{datetime.now()}] Adding newfound galaxy {star_gal_name}...")

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

    print(f"[{datetime.now()}] Done adding stars.")
    if savepath is not None:
        print(f"[{datetime.now()}] Saving catalog to JSON...")
        with open(savepath, 'w', encoding="utf-8") as f:
            json.dump(output, f, indent = 4)
    
    print(f"[{datetime.now()}] Done saving catalog.")
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
    # for gal in catalog.keys():


# Execution code
if __name__ == "__main__":
    Kirby_2009 = retrieve_catalog("J/ApJS/191/352/abun",
                                  f"{TEMP_FOLDER}/Kirby_2009.json")

    test = retrieve_catalog("J/ApJ/838/83")

    # print(Kirby_2009)
    # add_catalog(test, f"{TEMP_FOLDER}/Kirby_2009.json", "J/ApJS/191/352/abun")


# =============================================================================
# 1. Check if the coords from stars in that catalog match any of the existing database members

# Precondition for new incoming catalogs:
# they usually either look at several galaxies at once, or one galaxy at a time.
# If it's the latter, then we can just manually indicate which galaxy is that,
# and manually decide if that's a new galaxy.

#NOTE: Read the papers for each new catalog.

#NOTE: If the catalog doesn't mention galaxy names, it's a single-galaxy catalog.


