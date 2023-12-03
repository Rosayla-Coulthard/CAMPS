"""Reads data from a catalog file and stores it locally in a dictionary."""
import json
import numpy as np
from astroquery.vizier import Vizier
from datetime import datetime
from config import COLUMNS, TEMP_FOLDER, META_DATA_KEY
from utilfuncs import find_nearest_galaxy
from utilfuncs import dict_depth, galaxy_crossmatch

def retrieve_catalog(prompt:str, savepath = None, UI = 'Name', gal_title = 'dSph'):
    """Retrieves data from Vizier via one prompt.
    
    Args:
        prompt: The prompt to search for in Vizier.
        savepath: The path to save the catalog to.
        UI: The column name for the unique identifier.
        gal_title: The column name for the galaxy title.
    """
    V = Vizier(
        # columns = COLUMNS,
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
    
    # Ask if the catalog is a single-galaxy catalog
    msg = "Is this a single-galaxy catalog? (y/n): "
    single_galaxy = input(msg)
    if single_galaxy == 'y':
        single_galaxy = True
        msg = "Please specify a galaxy name: "
        gal_name = input(msg)
        catalog[gal_title] = gal_name
    elif single_galaxy == 'n':
        single_galaxy = False
    else:
        print(f"[{datetime.now()}] WARNING: Invalid input. Operation aborted, no changes made.")
        quit()
    
    print(f"[{datetime.now()}] Adding stars...")

    count = 0
    for star in catalog:
        # Fix this part
        try:
            star_gal_name = star[gal_title]
        except KeyError:
            # If the galaxy name is not specified, find the nearest galaxy
            print(f"[{datetime.now()}] Single galaxy catalog detected.")
            star_gal_name = input("Please specify a galaxy name: ")

        # Check if the star has a name, make one up if not
        if UI not in star.keys():
            star_name = f"{star_gal_name}_{count:04d}"
        else:
            star_name = star[UI]
        
        # Check if the star's galaxy is specified. If not, add one.
        try:
            output[star_gal_name][star_name] = {}
        except KeyError:
            if single_galaxy:
                output[star_gal_name] = {star_name: {}}
                print(f"[{datetime.now()}] Adding newfound galaxy {star_gal_name}...")

        # Add catalog under star
        output[star_gal_name][star_name][prompt] = {}

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
                try:
                    star_col = float(star[col])
                except ValueError:
                    star_col = star[col]

            output[star_gal_name][star_name][prompt][col] = star_col

        count += 1

    output[META_DATA_KEY]['Member count'] = count

    print(f"[{datetime.now()}] Done adding stars.")
    if savepath is not None:
        print(f"[{datetime.now()}] Saving catalog to JSON...")
        with open(savepath, 'w', encoding="utf-8") as f:
            json.dump(output, f, indent = 4)
    
    print(f"[{datetime.now()}] Done saving catalog.")
    return output

# Write a function to sort new stars into galaxies
def add_catalog(catalog_path, cache_path):
    """Takes a single-catalog JSON file and adds it to the cache."""
    # Load the catalog
    catalog = json.load(open(catalog_path, encoding="utf-8"))
    catalog_contents = list(catalog[META_DATA_KEY]['Included titles'].keys())

    print(catalog_contents)
    # Check if the incoming file is a single-catalog file
    if len(catalog_contents) > 1:
        raise ValueError("Incoming catalog must be a single-catalog file.")

    # Load the cache
    cache = json.load(open(cache_path, encoding="utf-8"))

    # Do galaxy crossmatching
    print(f"[{datetime.now()}] Crossmatching galaxies...")
    catalog_galaxies = list(catalog.keys())[1:]
    cache_galaxies = list(cache.keys())[1:]

    # Iterate through each galaxy in the catalog
    for gal_name in catalog_galaxies:
        # Check if the galaxy is already in the cache
        for gal_name_ in cache_galaxies:
            # If the galaxy is already in the cache, add the stars to the cache
            if gal_name == gal_name_:
                status_msg = f"[{datetime.now()}] Galaxy {gal_name} found in cache. " 
                status_msg += "Adding stars..."
                print(status_msg)
                for star_name in catalog[gal_name].keys():
                    if star_name == META_DATA_KEY:
                        continue
                    cache[gal_name][star_name] = catalog[gal_name][star_name]
                print(f"[{datetime.now()}] Done adding stars.")
                break


# Execution code
if __name__ == "__main__":
    add_catalog("Temp/Kirby_2009.json", "Data/cache.json")
    


# =============================================================================
# 1. Check if the coords from stars in that catalog match any of the existing database members

# Precondition for new incoming catalogs:
# they usually either look at several galaxies at once, or one galaxy at a time.
# If it's the latter, then we can just manually indicate which galaxy is that,
# and manually decide if that's a new galaxy.

#NOTE: Read the papers for each new catalog.

#NOTE: If the catalog doesn't mention galaxy names, it's a single-galaxy catalog.


