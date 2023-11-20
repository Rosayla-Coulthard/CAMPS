"""Utility functions go here"""
import numpy as np
import json
from config import TEMP_FOLDER

def find_nearest_galaxy(star_coords, catalog_path = None, catalog_name = None):
    """Finds the nearest star in a catalog to a given star."""
    if catalog_name is None or catalog_path is None:
        raise ValueError("Catalog name must be specified.")
    
    # Load JSON file
    catalog = json.load(open(catalog_path, encoding="utf-8"))

    # Find the nearest galaxy
    gal_dist = np.inf
    for gal in catalog.keys():
        first_star = list(catalog[gal].keys())[0]
        first_star = catalog[gal][first_star][catalog_name]

        first_star_coord = first_star["RAJ2000"], first_star["DEJ2000"]
        
        gal_dist_new = (star_coords["RAJ2000"] - first_star_coord[0])**2 + (star_coords["DEJ2000"] - first_star_coord[1])**2

        if gal_dist_new < gal_dist:
            gal_dist = gal_dist_new
            nearest_gal = gal
        
    return nearest_gal, gal_dist

if __name__ == "__main__":
    test_star = {
        "RAJ2000": 17.35407777777778,
        "DEJ2000": 57.869527777777776,
    }
    # _ = find_nearest_star(test_star, f"{TEMP_FOLDER}/Kirby_2009.json", 'J/ApJS/191/352/abun')

