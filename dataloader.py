"""Reads data from a catalog file and stores it locally in a dictionary."""
import json
import numpy as np
from astroquery.vizier import Vizier
from datetime import datetime
from utilfuncs import find_nearest_galaxy, save_to_json
from utilfuncs import dict_depth, galaxy_crossmatch

# Load the configs
config = json.load(open("config.json", encoding="utf-8"))
ALT_COL_NAMES = config["ALT_COL_NAMES"]
META_DATA_KEY = config["META_DATA_KEY"]

def retrieve_catalog(prompt:str, savepath = None, UI = 'Name', gal_title = 'dSph'):
    """Retrieves data from Vizier via one prompt.
    
    Args:
        prompt: The prompt to search for in Vizier.
        savepath: The path to save the catalog to.
        UI: The column name for the unique identifier.
        gal_title: The column name for the galaxy title.
    """
    V = Vizier(
        row_limit = -1
    )
    catalog_list = V.get_catalogs(prompt)
    catalog_keys = list(catalog_list.keys())
    
    for count, i in enumerate(catalog_keys):
        print(f'[{count}]', i)

    input_msg = "Which table above would you like to use? Type index: (Type anything else to quit)"
    catalog_choice = input(input_msg)

    try:
        catalog_choice = int(catalog_choice) 
        catalog = catalog_list[catalog_keys[catalog_choice]]
    except ValueError:
        print(f"[{datetime.now()}] WARNING: Invalid input. Operation aborted, no changes made.")
        quit()

    print(f"[{datetime.now()}] Found the following columns: {catalog.colnames}")

    # Refer to config file for the column name conversion
    col_list = {}
    for col in catalog.colnames:
        try:
            col_list[col] = ALT_COL_NAMES[col]
        except KeyError:
            msg = f"Found column {col} not in config file. Please specify a name: "
            col_list[col] = input(msg)

    output = {f"{META_DATA_KEY}": {"Included titles": {prompt: col_list}}}
    
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

    # Configure naming scheme
    if isinstance(UI, list) is False:
        UI_new = [UI]
    else:
        UI_new = UI

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
        if UI not in star.colnames or UI is not str:
            star_name = f"{star_gal_name}_{count:04d}"
        else:
            # If the star has a name, use it
            star_name = star[UI]
        
        # Check if the star's galaxy is specified. If not, add one.
        try:
            output[star_gal_name][star_name] = {}
        except KeyError:
            if single_galaxy:
                print(f"[{datetime.now()}] Adding newfound galaxy {star_gal_name}...")
            output[star_gal_name] = {star_name: {}}

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
            else:
                try:
                    star_col = float(star[col])
                except ValueError:
                    star_col = star[col]

            if col == 'dSph':
                new_col_name = 'dSph'
                star_col = star_gal_name
            else:
                new_col_name = col_list[col]

            # If the star has been labled as not a member, remove the star
            star_abandoned = False
            if new_col_name == 'Member' and star_col != 'Y' and star_col != 'y':
                output[star_gal_name].pop(star_name)
                star_abandoned = True
                break
            output[star_gal_name][star_name][prompt][new_col_name] = star_col

        if star_abandoned:
            continue

        star = output[star_gal_name][star_name]
        # Refresh star name
        if isinstance(UI_new, list):
            star = output[star_gal_name][star_name]
            ui_temp = [star[prompt][i] for i in UI_new]
            name_list = [str(i) for i in ui_temp]
            star_name_new = '_'.join(name_list)

            try:
                output[star_gal_name][star_name_new] = output[star_gal_name].pop(star_name)
            except KeyError:
                pass

        count += 1

    output[META_DATA_KEY]['Member count'] = count
    print(f"[{datetime.now()}] Done adding stars.")
    if savepath is not None:
        print(f"[{datetime.now()}] Saving catalog to JSON...")
        with open(savepath, 'w', encoding="utf-8") as f:
            json.dump(output, f, indent = 4)

    # Update column name settings
    with open("config.json", 'r', encoding="utf-8") as f:
        config_file = json.load(f)
    config_file["ALT_COL_NAMES"].update(col_list)
    save_to_json(config_file, "config.json")
    
    print(f"[{datetime.now()}] Done saving catalog.")
    return output


# Write a function to sort new stars into galaxies
def add_catalog(catalog_path, catalog_name, cache_path, ref_catalog):
    """Takes a single-catalog JSON file and adds it to the cache."""
    # Load the catalog
    catalog = json.load(open(catalog_path, encoding="utf-8"))
    catalog_contents = list(catalog[META_DATA_KEY]['Included titles'].keys())

    # print(catalog_contents)
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
        galaxy_found = False
        for gal_name_ in cache_galaxies:
            # If the galaxy is already in the cache, add the stars to the cache
            if gal_name == gal_name_:
                galaxy_found = True
                status_msg = f"[{datetime.now()}] Galaxy {gal_name} found in cache. "
                status_msg += "Adding stars..."
                print(status_msg)
                # NOTE: star crossmatch happens here, and add the stars to the cache

                break
            # If the galaxy is not in the cache, find the nearest galaxy
            # and add the stars to that galaxy

        if galaxy_found is False:
            status_msg = f"[{datetime.now()}] Galaxy {gal_name} not found in cache. "
            status_msg += "Finding nearest galaxy..."
            print(status_msg)

            # Find the nearest galaxy
            nearest_galaxy = galaxy_crossmatch(catalog[gal_name],catalog_name,ref_catalog,cache)

            if nearest_galaxy is None:
                # This happens when the galaxy was not previously recorded
                msg = f"[{datetime.now()}] No galaxy in the cach matches galaxy {gal_name}. "
                msg += "Adding galaxy to cache..."
                print(msg)
                cache[gal_name] = catalog[gal_name]
                continue
            
            #NOTE: star crossmatch happens here, AGAIN

        time = datetime.now()
        print(f"[{time}] Done adding stars.")

    # Add catalog to the list of included titles
    catalog_metadata = catalog[META_DATA_KEY]['Included titles'][catalog_name]
    cache[META_DATA_KEY]['Included titles'][catalog_name] = catalog_metadata

    # Update member count
    cache[META_DATA_KEY]['Member count'] += catalog[META_DATA_KEY]['Member count']

    # Save the cache
    print(f"[{datetime.now()}] Saving cache...")
    save_to_json(cache, cache_path)

    return cache

def merge_catalogs(catalog1_path, catalog2_path, catalog1_name, catalog2_name, alter_savepath = None):
    """Naively merges two catalogs by combining their data.
    If both catalogs have stars and galaxies under the same names,
        their data will be combined.
    If both catalogs have stars and galaxies under different names,
        their data will be logged separately. 
    
    Args:
        catalog1: Save path to the receiving catalog
        catalog2: Save path to the catalog to be merged into the receiving catalog
        catalog1_name: The name of the receiving catalog
        catalog2_name: The name of the other catalog
    
    Returns:
        The merged catalog
    """
    # Load the catalogs
    catalog1 = json.load(open(catalog1_path, encoding="utf-8"))
    catalog2 = json.load(open(catalog2_path, encoding="utf-8"))

    # Check each of their depths
    catalog1_depth = dict_depth(catalog1)
    catalog2_depth = dict_depth(catalog2)

    if catalog1_depth != catalog2_depth:
        raise ValueError("Catalogs must have the same depth.")
    if catalog1_depth != 4:
        msg = "Need 4-level catalogs. "
        msg += f"Catalog 1 has {catalog1_depth} levels. "
        msg += f"Catalog 2 has {catalog2_depth} levels."
        raise ValueError(msg)
    
    for gal_name in catalog2:
        if gal_name == META_DATA_KEY:
            catalog1_columns = catalog1[gal_name]['Included titles'][catalog1_name]
            catalog2_columns = catalog2[gal_name]['Included titles'][catalog2_name]

            new_columns = catalog1_columns | catalog2_columns

            catalog1[gal_name]['Included titles'][catalog1_name] = new_columns

            continue
        
        gal = catalog2[gal_name]
        for star2_name in gal:
            star2 = gal[star2_name]
            try:
                catalog1_gal = catalog1[gal_name]
            except KeyError:
                continue
            if star2_name in catalog1_gal.keys():
                star1 = catalog1_gal[star2_name]
                new_star = star1 | star2
                catalog1_gal[star2_name] = new_star
            else:
                catalog1_gal[star2_name] = star2

    # Save the merged catalog
    savepath = alter_savepath if alter_savepath is not None else catalog1_path
    save_to_json(catalog1, savepath)

    return catalog1
                

# Execution code
if __name__ == "__main__":
    add_catalog("Temp/J_ApJ_838_83.json", "J/ApJ/838/83", 'Data/cache.json', "J/ApJS/191/352/abun")


# =============================================================================
# 1. Check if the coords from stars in that catalog match any of the existing database members

# Precondition for new incoming catalogs:
# they usually either look at several galaxies at once, or one galaxy at a time.
# If it's the latter, then we can just manually indicate which galaxy is that,
# and manually decide if that's a new galaxy.

#NOTE: Read the papers for each new catalog.

#NOTE: If the catalog doesn't mention galaxy names, it's a single-galaxy catalog.


