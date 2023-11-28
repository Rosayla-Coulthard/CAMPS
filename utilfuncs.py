"""Utility functions go here"""
import numpy as np
import json
from config import META_DATA_KEY, TEMP_FOLDER

# Catch errors
def catalog_type_check(catalog):
    """Checks if the catalog is a path to a JSON file or a dictionary."""
    try:
        catalog = json.load(open(catalog, encoding="utf-8"))
    except TypeError:
        if type(catalog) != dict:
            raise TypeError("Catalog must be a path to a JSON file or a dictionary.")

    return catalog
            
def catalog_name_check(catalog, catalog_name):
    """Checks if a specified catalog had already been compiled."""
    collection = list(catalog[META_DATA_KEY]['Included titles'].keys())
    
    if catalog_name not in collection:
        raise ValueError("Catalog was never added to the collection.")
    elif catalog_name is None:
        raise ValueError("Catalog name must be specified.")

    return catalog_name

def find_nearest_galaxy(star_coords, catalog = str or dict, catalog_name = None):
    """Finds the nearest star in a catalog to a given star.
    
    Args:
        star_coords: A dictionary containing the coordinates of the incoming star.
        catalog_path: The path to the catalog file.
        catalog_name: The name of the catalog in the catalog file.
    """
    # Load JSON file and check type
    catalog = catalog_type_check(catalog)

    catalog_name = catalog_name_check(catalog, catalog_name)

    # Find the nearest galaxy
    gal_dist = np.inf
    for gal in catalog.keys():
        if gal == META_DATA_KEY:
            continue

        first_star = list(catalog[gal].keys())[0]
        first_star = catalog[gal][first_star][catalog_name]

        first_star_coord = first_star["RAJ2000"], first_star["DEJ2000"]

        gal_dist_new = (star_coords["RAJ2000"] - first_star_coord[0])**2 + (star_coords["DEJ2000"] - first_star_coord[1])**2

        if gal_dist_new < gal_dist:
            gal_dist = gal_dist_new
            nearest_gal = gal

    # Find the nearest galaxy's centroid and boundary,
    # and check if the star is within the boundary
    gal_member_coords = []
    for star in catalog[nearest_gal].keys():
        star = catalog[nearest_gal][star][catalog_name]
        ra_temp = np.float64(star['RAJ2000'])
        dec_temp = np.float64(star['DEJ2000'])

        gal_member_coords.append([ra_temp, dec_temp])

    ra_, dec_ = np.transpose(gal_member_coords)

    ra_bound = (max(ra_) - min(ra_)) * 0.5
    dec_bound = (max(dec_) - min(dec_)) * 0.5

    # Check if the star is within the boundary
    star_ra, star_dec = star_coords["RAJ2000"], star_coords["DEJ2000"]
    star_dist = (star_ra - (max(ra_) + min(ra_)) * 0.5)**2
    star_dist += (star_dec - (max(dec_) + min(dec_)) * 0.5)**2

    if star_dist > ra_bound**2 + dec_bound**2:
        nearest_gal = None
        gal_dist = None

    return nearest_gal, gal_dist

# Test the functions
if __name__ == "__main__":
    test_star = {
        "RAJ2000": 1.9893277777777777,
        "DEJ2000": -31.28391666666666,
    }
    near_gal, near_dist = find_nearest_galaxy(test_star,
                                              catalog=f"{TEMP_FOLDER}/Kirby_2009.json",
                                              catalog_name = 'J/ApJS/191/352/abun')
