"""Reads data from a catalog file and stores it locally in a dictionary."""
from astroquery.vizier import Vizier
from astropy.table import Table
from astropy.io import fits
from config import COLUMNS
import config as c
import json

def retreive_catalog(prompt:str) -> dict:
    """Retreives a catalog from a file.
    
    Args:
        prompt: The prompt for searching on ViZier.
    
    Returns:
        A dict containing the catalog.
    """
    V = Vizier(
        columns = COLUMNS,
        row_limit = -1
    )
    catalog_list = V.get_catalogs(prompt)

    print(f"Found {len(catalog_list)} catalogs.")
    catalog_names = catalog_list.keys()

    catalogs_found = {}
    for count, name in enumerate(catalog_names):
        # Find the catalog with all the columns we want
        catalog_columns = catalog_list[name].keys()
        if len(catalog_columns) == len(COLUMNS):
            catalogs_found[name] = catalog_list[name]

    if count == len(catalog_names):
        raise NameError("No catalog found with all the columns we want.")
    
    print(len(catalogs_found), "catalogs found with all the columns we want.")
    return catalogs_found

def convert_catalog(catalog:Table, catalog_name:str, save_path = None) -> dict:
    """Converts a catalog to a dictionary, and saves it to a JSON file.
    
    Args:
        catalog: The catalog to convert.
        catalog_name: The name of the catalog.
        save_path: The path to save the catalog to.
    """
    print("Rearranging data for saving...")
    star_dict = {}
    catalog_columns = catalog.keys()

    for star in catalog:
        gal_name = star["dSph"]

        try: # Assign stars to galaxies. If no such galaxy, make one
            star_dict[gal_name][star["Name"]] = {catalog_name: {}}
        except KeyError:
            star_dict[gal_name] = {star["Name"]: {catalog_name: {}}}

        for column in catalog_columns:
            # Check if column name is mentioned in the alternative column names
            if column in c.ALT_COL_NAMES.keys():
                column = c.ALT_COL_NAMES[column]
            star_dict[gal_name][star["Name"]][catalog_name][column] = str(star[column])

            # If the colunns are RA and DEC, convert them to floats
            if column == "RAJ2000" or column == "DEJ2000":
                coord = str(star[column]).split(" ")
                coord = float(coord[0]) + float(coord[1]) /60 + float(coord[2]) /3600
                star_dict[gal_name][star["Name"]][catalog_name][column] = str(coord)


    print("Saving data to JSON...")
    if save_path:
        with open(save_path, "w+", encoding = "utf-8") as f:
            f.write(json.dumps(star_dict, indent=4))

    print(f"Catalog '{catalog_name}' saved to JSON.")
    return star_dict

# Execution code
if __name__ == "__main__":
    Kirby_2009 = retreive_catalog("Kirby 2009")

    for name, item in Kirby_2009.items():
        catalog = item
        convert_catalog(catalog, name, f"{c.TEMP_FOLDER}/Kirby 2009.json")
