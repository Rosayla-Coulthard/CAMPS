"""Utility functions go here"""
from datetime import datetime
import json
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky

# Import configurations
config = json.load(open("config.json", encoding="utf-8"))
META_DATA_KEY = config["META_DATA_KEY"]
TEMP_FOLDER = config["TEMP_FOLDER"]

# Catch errors
def catalog_type_check(catalog):
    """Checks if the catalog is a path to a JSON file or a dictionary.
    
    Args:
        catalog: The catalog to be checked.

    Returns:
        The catalog as a dictionary.
    """
    try:
        catalog = json.load(open(catalog, encoding="utf-8"))
    except TypeError as e:
        if type(catalog) != dict:
            raise TypeError("Catalog must be a path to a JSON file or a dictionary.") from e

    return catalog
            
def catalog_name_check(catalog, catalog_name):
    """Checks if a specified catalog had already been compiled.
    
    Args:
        catalog: The catalog to be checked.
        catalog_name: The name of the catalog to be checked.

    Returns:
        The name of the catalog if it is in the collection.
    """
    collection = list(catalog[META_DATA_KEY]['Included titles'].keys())
    
    if catalog_name not in collection:
        raise ValueError("Catalog was never added to the collection.")
    elif catalog_name is None:
        raise ValueError("Catalog name must be specified.")

    return catalog_name

def dict_depth(d):
    """Finds the depth of a nested dictionary.
    
    Args:
        d: The dictionary to be checked.

    Returns:
        The depth of the dictionary.
    """
    if not isinstance(d, dict) or not d:
        # If not a dictionary or an empty dictionary, depth is 0
        return 0

    # Recursive call to find the depth of nested dictionaries
    depths = [dict_depth(value) for value in d.values()]

    # Return the maximum depth plus 1 for the current dictionary
    return max(depths) + 1

def save_to_json(data, savepath):
    """Saves data to a JSON file."""
    with open(savepath, 'w', encoding="utf-8") as f:
        json.dump(data, f, indent = 4)



def find_nearest_galaxy(star_coords, catalog = str or dict, ref_catalog = None):
    """Finds the nearest star in a catalog to a given star.
    
    Args:
        star_coords: A dictionary containing the coordinates of the incoming star.
        catalog_path: The path to the catalog file.
        ref_catalog: Reference catalog to cross match with.
    
    Returns:
        nearest_gal: The name of the nearest galaxy.
        gal_dist: The distance between the incoming star and the nearest galaxy.
    """
    # Load JSON file and check type
    print(f"[{datetime.now()}] Finding the nearest galaxy...")
    catalog = catalog_type_check(catalog)

    ref_catalog = catalog_name_check(catalog, ref_catalog)
    # Find the nearest galaxy
    gal_dist = np.inf

    # NOTE: this is where the problem is
    # dict_depth_ = dict_depth(catalog)
    for gal in catalog.keys():
        if gal == META_DATA_KEY:
            continue

        first_star = list(catalog[gal].keys())[0]
        first_star = catalog[gal][first_star][ref_catalog]

        first_star_coord = first_star["RAJ2000"], first_star["DEJ2000"]

        gal_dist_new = (star_coords["RAJ2000"] - first_star_coord[0])**2 
        gal_dist_new += (star_coords["DEJ2000"] - first_star_coord[1])**2

        if gal_dist_new < gal_dist:
            gal_dist = gal_dist_new
            nearest_gal = gal

    # Find the nearest galaxy's centroid and boundary,
    # and check if the star is within the boundary
    gal_member_coords = []
    for star in catalog[nearest_gal].keys():
        star = catalog[nearest_gal][star][ref_catalog]
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
        time_stamp = datetime.now()
        print(f"[{time_stamp}] Star is not within the boundary of the nearest galaxy, returning None.")

    return nearest_gal, gal_dist

def galaxy_crossmatch(gal1:dict, gal1_catalog:str, ref_catalog, cache):
    """See if gal1 is already in the cache. If so, return gal_name. If not,
        return None. Just provide a galaxy in dict format, and the catalog
        it is from.

        This is done by providing a star from gal1, and see if it is within 
        the boundary of any galaxy in the cache.
    
    Args:
        gal1: The galaxy to be crossmatched.
        gal1_catalog: The name of the catalog for gal1.
        ref_catalog: Reference catalog to cross match with.
        cache: The cache to be crossmatched against.
        
    Returns:
        gal_name if gal1 is already in the cache, None otherwise.
        """
    # Provides a star from gal1, and see if it is within the boundary 
    # of any galaxy in the cache
    first_star = list(gal1.keys())[0]
    first_star = gal1[first_star][gal1_catalog]
    print(first_star)

    first_star_coord = {"RAJ2000": first_star["RAJ2000"],
                        "DEJ2000": first_star["DEJ2000"]}
    gal_name, _ = find_nearest_galaxy(first_star_coord, cache, ref_catalog)

    return gal_name

def star_crossmatch(gal1:dict, gal1_catalog:str, gal2:dict, gal2_catalog:str):
    """Main function for crossmatching stars between two galaxies in two catalogs.
    
        Precondition: the two galaxies are already crossmatched.
    """
    # Get the coordinates of all stars in gal1
    gal1_coords = []
    for star in gal1.keys():
        star = gal1[star][gal1_catalog]
        ra_temp = np.float64(star['RAJ2000'])
        dec_temp = np.float64(star['DEJ2000'])

        gal1_coords.append([ra_temp, dec_temp])

    # Get the coordinates of all stars in gal2
    gal2_coords = []
    for star in gal2.keys():
        star = gal2[star][gal2_catalog]
        ra_temp = np.float64(star['RAJ2000'])
        dec_temp = np.float64(star['DEJ2000'])

        gal2_coords.append([ra_temp, dec_temp])

    # Crossmatch the two galaxies
    gal1_coords = SkyCoord(gal1_coords, unit="deg")
    gal2_coords = SkyCoord(gal2_coords, unit="deg")

    idx, d2d, d3d = match_coordinates_sky(gal1_coords, gal2_coords)

    # Find the closest star in gal2 for each star in gal1
    gal1_star_names = list(gal1.keys())
    gal2_star_names = list(gal2.keys())

    gal1_star_names = np.array(gal1_star_names)[idx]
    gal2_star_names = np.array(gal2_star_names)

    gal1_star_names = gal1_star_names.tolist()
    gal2_star_names = gal2_star_names.tolist()

    # Add the closest star in gal2 to gal1
    for i in range(len(gal1_star_names)):
        gal1[gal1_star_names[i]][gal1_catalog]["Closest star"] = gal2_star_names[i]

    # Add the closest star in gal1 to gal2
    for i in range(len(gal2_star_names)):
        gal2[gal2_star_names[i]][gal2_catalog]["Closest star"] = gal1_star_names[i]

    return gal1, gal2

# Test the functions
if __name__ == "__main__":
    # Load catalog and cache
    catalog = json.load(open("Temp/J_ApJ_838_83.json", encoding="utf-8"))
    cache = json.load(open("Data/cache.json", encoding="utf-8"))

    gal1 = cache['Scl']
    gal2 = cache['Scl']

    star_crossmatch(gal1, "J/ApJS/191/352/abun", gal2, "J/ApJS/191/352/abun")
