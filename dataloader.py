"""Reads data from a catalog file and stores it locally in a dictionary."""
import json
from json.decoder import JSONDecodeError
import numpy as np
from astroquery.vizier import Vizier
from datetime import datetime
from utilfuncs import save_to_json, member_count, star_crossmatch
from utilfuncs import dict_depth, galaxy_crossmatch

# Load the configs
config = json.load(open("config.json", encoding="utf-8"))
ALT_COL_NAMES = config["ALT_COL_NAMES"]
META_DATA_KEY = config["META_DATA_KEY"]

def retrieve_catalog(prompt:str, savepath = None, UI = 'Name', gal_title = 'Galaxy',
                     data_source = "Vizier",
                     table_num = None, paper_name = None, single_galaxy = None, gal_name = None):
    """Retrieves data from Vizier via one prompt.
    
    Args:
        prompt: The prompt to search for in Vizier.
        savepath: The path to save the catalog to.
        UI: The column name for the unique identifier. If there are multiple columns,
            use a list of strings. Use the same names as Vizier.
        gal_title: The column name for the galaxy title.

    Returns:
        The catalog as a dictionary.e
    """
    # Retrieve data from Vizier
    print(f"Retrieving catalog {prompt}...")

    if data_source == "Vizier":
        V = Vizier(
        row_limit = -1,
        columns=['*', '_RAJ2000', '_DEJ2000']
        )
    else:
        raise ValueError("Data source not recognized. Not implemented.")

    catalog_list = V.get_catalogs(prompt)
    catalog_keys = list(catalog_list.keys())

    # Ask user to choose a table. If there is only one table, use that table.
    for count, i in enumerate(catalog_keys):
        print(f'[{count}]', i)
        count_temp = count
        if count_temp > config['MAX_TABLES']:
            msg = "WARNING: Too many tables are found."
            msg += "Please refine your search, or increase MAX_TABLES in config.json."
            print(msg)
            print("Operation aborted, no changes made.")
            raise ValueError(msg)

    if count_temp == 0:
        msg = "Only one table found. Using that table."
        catalog_choice = 0
        print(msg)
    elif table_num is not None:
        catalog_choice = table_num
    else:
        input_msg = "Which table above would you like to use? Type index: (Type anything else to quit)"
        catalog_choice = input(input_msg)

    # Check if the input is valid. Either use that table, or quit.
    try:
        catalog_choice = int(catalog_choice)
        table_name = catalog_keys[catalog_choice]
        catalog = catalog_list[table_name]
    except ValueError:
        print("WARNING: Invalid input. Operation aborted, no changes made.")
        raise ValueError("Invalid input.") from None

    # Indicates the columns in the catalog
    print(f"Found the following columns: {catalog.colnames}")

    # Refer to config file for the column name conversion
    col_list = {}
    new_column_found = False
    for col in catalog.colnames:
        try:
            col_list[col] = ALT_COL_NAMES[col]
        except KeyError:
            # If the column is not in the config file, ask the user to specify a name
            new_column_found = True
            msg_ = f"(enter _ to not include column '{col}' in the catalog)"
            msg = f"Found column '{col}' not in config file. Please specify a name {msg_}: "
            user_input = input(msg)
            col_list[col] = user_input
            config["ALT_COL_NAMES"][col] = col_list[col]

    # Update the config file with the newfound column names
    if new_column_found:
        save_to_json(config, "config.json")
        print("Updated config file with new column names.")
        retrieve_catalog(prompt, savepath, UI, gal_title, table_num=catalog_choice)

    # Initiate the output dictionary with metadata
    # if paper_name is not None:
    if paper_name is None:
        cat_title = prompt
    else:
        cat_title = paper_name

    output = {f"{META_DATA_KEY}": {"Included titles": {cat_title: {}}}}
    output[META_DATA_KEY]["Included titles"][cat_title]["Columns"] = list(col_list.values())
    
    # Add the table info to the metadata
    table_name = table_name.split("/")
    catalog_name = "/".join(table_name[:-1])
    table_name = table_name[-1]

    table_info = {"Data source": data_source, "Paper tag": catalog_name, "Tables used": [table_name]}

    output[META_DATA_KEY]["Included titles"][cat_title]["Table info"] = table_info

    # Ask if the catalog is a single-galaxy catalog
    msg = "Is this a single-galaxy catalog? (y/n): "
    if single_galaxy is None:
        single_galaxy = input(msg)
    if single_galaxy == 'y' or single_galaxy == True:
        # If it is, ask for the galaxy name
        single_galaxy = True
        msg = "Please specify a galaxy name: "
        if gal_name is None:
            gal_name = input(msg)
        # catalog[gal_title] = gal_name
    elif single_galaxy == 'n' or single_galaxy == False:
        # If it is not, default to the galaxy name specified in the catalog
        single_galaxy = False
    else:
        # If the input is invalid, quit
        print("WARNING: Invalid input. Operation aborted, no changes made.")
        raise ValueError("Invalid input.") from None

    # Configure naming scheme
    if isinstance(UI, list) is False:
        UI_new = [UI]
    else:
        UI_new = UI

    # Adding stars after this point
    print("Adding stars...")
    count = 0
    for star in catalog:
        if single_galaxy:
            # If the star belongs to a single galaxy catalog, the galaxy name is
            # already specified above.
            star_gal_name = gal_name
        elif single_galaxy is False:
            # If the star belongs to a multi-galaxy catalog, the galaxy name is
            # specified in the catalog.
            try:
                star_gal_name = star[gal_title].replace(' ', '')
            except KeyError:
                print(gal_title)
                err_msg = "Either this is not a multi-galaxy catalog,"
                err_msg += "or the galaxy name is not specified in the input file."
                raise KeyError(err_msg) from None

        # Check if the star has a name, make one up if not
        if UI_new is None:
            star_name = f"{gal_name}{count:04d}"
        ui_temp = [star[i] for i in UI_new]
        name_list = [str(i) for i in ui_temp]
        star_name = '_'.join(name_list)

        # If the star belongs to a multi-galaxy catalog, and the galaxy name is not
        # previously recorded, make a new dictionary entry for the galaxy.
        if star_gal_name not in output.keys():
            output[star_gal_name] = {star_name: {cat_title: {"Name": star_name}}}
        elif star_gal_name in output.keys():
            output[star_gal_name][star_name] = {cat_title: {"Name": star_name}}

        # Add catalog under star
        for col in catalog.colnames:
            if col == "RAJ2000"or col == "DEJ2000":
                star_col = star[col]
            elif col == "Name":
                star_col = star_name
            elif col == 'dSph':
                new_col_name = 'dSph'
                star_col = star_gal_name
            else:
                try:
                    star_col = float(star[col])
                except ValueError:
                    star_col = star[col]

            new_col_name = col_list[col]

            # If the star has been labled as not a member, remove the star
            star_abandoned = False
            if new_col_name == 'Member' and star_col != 'Y' and star_col != 'y':
                output[star_gal_name].pop(star_name)
                star_abandoned = True
                break
            if new_col_name != "_":
                output[star_gal_name][star_name][cat_title][new_col_name] = star_col

        if star_abandoned:
            continue

        output[star_gal_name][star_name][cat_title]["Galaxy"] = star_gal_name

    # Update member count to metadata
    output, _, _ = member_count(output)
    galaxy_members = list(output.keys())[1:]
    output[META_DATA_KEY]["Included titles"][cat_title]["Members"] = galaxy_members

    print("Done adding stars.")

    # Save the data
    if savepath is not None:
        print("Saving catalog to JSON...")
        with open(savepath, 'w', encoding="utf-8") as f:
            print(type(output))
            json.dump(output, f, indent = 4)

    # Update column name configs
    with open("config.json", 'r', encoding="utf-8") as f:
        config_file = json.load(f)

    # del col_list['Members']
    save_to_json(config_file, "config.json")

    print("Done saving catalog.")
    return output

def merge_tables(catalog1_path, catalog2_path, catalog1_name, catalog2_name, alter_savepath = None):
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
        alter_savepath: The save path for the merged catalog. If None, the receiving catalog will be overwritten.
    
    Returns:
        The merged catalog
    """
    if catalog1_name != catalog2_name:
        raise ValueError("Why are you merging two different catalogs?")
    
    # Load the catalogs
    catalog1 = json.load(open(catalog1_path, encoding="utf-8"))
    catalog2 = json.load(open(catalog2_path, encoding="utf-8"))

    # Check each of their depths
    catalog1_depth = dict_depth(catalog1)
    catalog2_depth = dict_depth(catalog2)

    # Check if the catalogs are formatted correctly
    if catalog1_depth != catalog2_depth:
        raise ValueError("Catalogs must have the same depth.")
    if catalog1_depth != 5:
        msg = "Need 4-level catalogs. "
        msg += f"Catalog 1 has {catalog1_depth} levels. "
        msg += f"Catalog 2 has {catalog2_depth} levels."
        raise ValueError(msg)

    # Iterate through each galaxy in the incoming catalog
    for gal_name in catalog2:
        if gal_name == META_DATA_KEY: # If it's metadata, simply merge them
            catalog1_temp = catalog1[META_DATA_KEY]['Included titles'][catalog1_name]
            catalog2_temp = catalog2[META_DATA_KEY]['Included titles'][catalog2_name]

            # new_columns = catalog1_columns | catalog2_columns
            new_columns = list(set(catalog1_temp['Columns'] + catalog2_temp['Columns']))
            new_members = list(set(catalog1_temp['Members'] + catalog2_temp['Members']))

            catalog1_info = catalog1[META_DATA_KEY]['Included titles'][catalog1_name]["Table info"]
            catalog2_info = catalog2[META_DATA_KEY]['Included titles'][catalog2_name]["Table info"]

            catalog1_tables = catalog1_info["Tables used"]
            catalog2_tables = catalog2_info["Tables used"]

            new_tables = list(set(catalog1_tables + catalog2_tables))
            catalog1_info["Tables used"] = new_tables

            new_cat_metadata = {"Columns": new_columns, "Members": new_members,
                                "Table info": catalog1_info}

            catalog1[META_DATA_KEY]['Included titles'][catalog1_name] = new_cat_metadata

            continue

        gal = catalog2[gal_name]
        for star2_name in gal:
            star2 = gal[star2_name]
            try:
                catalog1_gal = catalog1[gal_name]
            except KeyError:
                # KeyError happens when gal_name doesn't exist in catalog1
                # When that happens, just copy the whole thing over.
                print(f"Galaxy {gal_name} not found in catalog1. Adding galaxy...")
                catalog1[gal_name] = catalog2[gal_name].copy()
                catalog1_gal = catalog1[gal_name]
            if star2_name in catalog1_gal.keys():
                new_star = catalog1_gal[star2_name].copy()

                for cat in star2:
                    if cat == catalog1_name:
                        new_star[cat].update(star2[cat])
                    else:
                        new_star[cat] = star2[cat]

                catalog1_gal[star2_name] = new_star
            else: # Automatically adds the star if it wasn't in catalog1
                print(f"Star {star2_name} not found in {catalog1_path}. Adding star...")
                catalog1_gal[star2_name] = star2

    # Save the merged catalog
    savepath = alter_savepath if alter_savepath is not None else catalog1_path
    save_to_json(catalog1, savepath)

    print(f"Done merging {catalog2_path} into {catalog1_path}. \n")
    return catalog1

# Write a function to sort new stars into galaxies
def add_catalog(catalog_path, catalog_name, cache_path, ref_catalog):
    """Takes a single-catalog JSON file and adds it to the cache.
    
    Args:
        catalog_path: The path to the catalog to be added.
        catalog_name: The name of the catalog to be added.
        cache_path: The path to the cache.
        ref_catalog: The name of the reference catalog. MUST ALREADY BE IN THE CACHE.

    Returns:
        The updated cache.
    """
    # Load the catalog
    catalog = json.load(open(catalog_path, encoding="utf-8"))
    catalog_contents = list(catalog[META_DATA_KEY]['Included titles'].keys())

    # Check if the incoming file is a single-catalog file
    if len(catalog_contents) > 1:
        raise ValueError("Incoming catalog must be a single-catalog file.")

    # Load the cache
    try:
        cache = json.load(open(cache_path, encoding="utf-8"))

    # If the cache is empty, copy the reference catalog to the cache
    except JSONDecodeError:
        print("Cache is empty. Copying reference catalog...")
        cache = {META_DATA_KEY: {"Included titles": {}}}

        for gal in catalog:
            gal_temp = catalog[gal]

            if gal == META_DATA_KEY:
                cache[META_DATA_KEY] = catalog[META_DATA_KEY]
                continue
            count = 0

            cache[gal] = {}
            for star in gal_temp:
                # Generate name for the star
                star_name = f"{gal}_{count}"
                count += 1

                # Add the star to the cache
                cache[gal][star_name] = gal_temp[star]
                cache[gal][star_name][catalog_name]["Star ID"] = star_name

        # Save the cache
        save_to_json(cache, cache_path)
        print("Done copying reference catalog to cache.")
        return cache

    # Do galaxy crossmatching
    print("Crossmatching galaxies...")
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
                status_msg = f"Galaxy {gal_name} found in cache. "
                status_msg += "Crosmatching and adding stars..."
                print(status_msg)

                # star crossmatch happens here, and add the stars to the cache
                incoming_gal = catalog[gal_name]
                ref_gal = cache[gal_name_]

                matched_gal, _ = star_crossmatch(incoming_gal, catalog_name,
                                              ref_gal, ref_catalog, matched_gal_name = gal_name_)

                # Add the updated galaxy to the cache
                cache[gal_name_] = matched_gal # BUG: this step somehow converted dict to tuple. Fixed for now, but need to investigate later
                break
            # If the galaxy is not in the cache, find the nearest galaxy
            # and add the stars to that galaxy
        if galaxy_found is False:
            status_msg = f"Galaxy {gal_name} not found in cache. "
            # status_msg += "Finding nearest galaxy..."
            print(status_msg)

            # Find the nearest galaxy
            nearest_galaxy = galaxy_crossmatch(catalog[gal_name],catalog_name,
                                               ref_catalog, cache)

            if nearest_galaxy is None:
                # This happens when the galaxy was not previously recorded
                msg = f"No galaxy in the cach matches galaxy {gal_name}. "
                msg += "Adding galaxy to cache..."
                print(msg)

                # Add the galaxy to the cache with updated star names
                count = 0
                new_gal = {}
                for star in catalog[gal_name]:
                    star_name = f"{gal_name}_{count}"
                    count += 1
                    new_gal[star_name] = catalog[gal_name][star]
                    new_gal[star_name][catalog_name]["Star ID"] = star_name

                cache[gal_name] = new_gal

                print(f"Done adding {gal_name} to cache.\n")
                continue
            else:
                # star crossmatch happens here, AGAIN
                # This happens when the galaxy was previously recorded, but under a different name
                msg = f"Found {nearest_galaxy} as the nearest galaxy to {gal_name}. "
                incoming_gal = catalog[gal_name]
                ref_gal = cache[nearest_galaxy]
                matched_gal, _ = star_crossmatch(incoming_gal, catalog_name,
                                              ref_gal, ref_catalog, matched_gal_name = gal_name_)

                # Add the updated galaxy to the cache
                cache[nearest_galaxy] = matched_gal

        print(f"Done adding {gal_name} to cache.\n")

    # Add catalog to the list of included titles
    catalog_metadata = catalog[META_DATA_KEY]['Included titles'][catalog_name]
    cache[META_DATA_KEY]['Included titles'][catalog_name] = catalog_metadata

    # Update member count
    cache, gal_count, star_count = member_count(cache)
    titles_count = len(cache[META_DATA_KEY]['Included titles'].keys())
    msg = f"The updated cache has {gal_count} galaxies and {star_count} stars"
    msg += f"from {titles_count} catalogs."
    print(msg)

    # Final pass to ensure all stars have unique names
    for gal in cache:
        if gal == "Meta data":
            continue

        for star in cache[gal]:
            for paper in cache[gal][star]:
                try:
                    cache[gal][star][paper]['Star ID']
                except KeyError:
                    cache[gal][star][paper]['Star ID'] = star

    # Save the cache
    print("Saving cache...")
    save_to_json(cache, cache_path)

    print(f"Added {catalog_name} to cache.\n")

    return cache




# Execution code
if __name__ == "__main__":
    # add_catalog("Temp/Kirby 2009.json", "Kirby+ 2010", "Data/cache.json", "Kirby+ 2010")
    # add_catalog("Temp/Theler 2020.json", "Theler+ 2020", "Data/cache.json", "Kirby+ 2010")
    # add_catalog("Temp/Reichert 2020.json", "Reichert+ 2020", "Data/cache.json", "Kirby+ 2010")
    # add_catalog("Temp/Kirby+, 2020.json", "Kirby+, 2020", "Data/cache.json", "Kirby+ 2010")
    # add_catalog("Temp/Kirby+ Carbon.json", "Kirby+ Carbon", "Data/cache.json", "Kirby+ 2010")
    # add_catalog("Temp/Kirby+, 2017.json", "Kirby+, 2017", "Data/cache.json", "Reichert+ 2020")
    # add_catalog("Temp/Kirby+ 2017 multigal.json", "Kirby+ 2017 multigal", "Data/cache.json", "Reichert+ 2020")
    # add_catalog("Temp/Duggan+, 2018.json", "Duggan+, 2018", "Data/cache.json", "Reichert+ 2020")
    # add_catalog("Temp/de los Reyes 2020.json", "de los Reyes 2020", "Data/cache.json", "Reichert+ 2020")
    retrieve_catalog("J/A+A/642/A176/", "Temp/Theler(8).json", UI="ID", paper_name="Theler 2020", single_galaxy=True, gal_name="Sex")
    merge_tables("Temp/Theler 2020.json", "Temp/Theler(8).json", 
                 "Theler 2020", "Theler 2020", alter_savepath="Temp/Theler 2020.json")
    # Sex, UMi, and Dra: show what catalogs I have about those galaxies,
    # which abundances are included, and how many stars.

    # Retreive the info as a panda file, or fits, or a data table format of some kind
    # Oh and the csv thing.

    # Make a "galaxy" called "field" that just records field stars.
