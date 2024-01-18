"""Makes plots."""
import numpy as np
from matplotlib import pyplot as plt
import json
from utilfuncs import catalog_name_check, META_DATA_KEY

def data_by_galaxy(cache_path:str, catalog_name:str, galaxies, columns):
    """Plots a pair of specified columns for each galaxy in the list, with
    all the data coming from one catalog.
    
    Args:
        cache_path (str): The path to the cache file.
        catalog_name (str): The name of the catalog to use.
        galaxies (list): A list of galaxy names.
        columns (list): A list of column names.
    
    Returns:
        output (dict): A dictionary of the data, organized by galaxy.
    """
    # Check inputs
    if len(columns) != 2:
        raise ValueError("columns must have exactly two elements")

    if META_DATA_KEY in columns:
        raise ValueError("columns cannot contain the meta data")

    # Load cache
    with open(cache_path, "r", encoding="utf-8") as cache_file:
        cache = json.load(cache_file)

    catalog_name = catalog_name_check(cache, catalog_name)

    # Read data
    col1_name, col2_name = columns
    output = {}
    for gal_name in galaxies:
        col1, col2 = [], []
        if gal_name is META_DATA_KEY:
            continue
        if gal_name not in cache.keys():
            raise KeyError("Galaxy not found in catalog")

        for star_name in cache[gal_name].keys():
            star = cache[gal_name][star_name]

            if catalog_name not in star.keys():
                continue
            if col1_name not in star[catalog_name].keys():
                continue
            if col2_name not in star[catalog_name].keys():
                continue

            col1.append(star[catalog_name][col1_name])
            col2.append(star[catalog_name][col2_name])

        if len(col1) == 0:
            print(f"{gal_name} not found in {catalog_name}, skipping")
            continue
        output[gal_name] = {col1_name: col1, col2_name: col2}


    return output

def data_by_catalog(cache_path:str, galaxy_name:str, catalogs:list, columns):
    """Plots a pair of specified columns for each catalog in the list,
    with all the data coming from one galaxy.

    Args:
        cache_path (str): The path to the cache file.
        galaxy_name (str): The name of the galaxy to use.
        catalogs (list): A list of catalog names.
        columns (list): A list of column names.

    Returns:
        ouptut (dict): A dictionary of the data, organized by catalog.
    """
    # Check inputs
    if len(columns) != 2:
        raise ValueError("columns must have exactly two elements")

    if galaxy_name is META_DATA_KEY:
        raise ValueError("galaxy_name cannot be the meta data")

    # Load cache
    with open(cache_path, "r", encoding="utf-8") as cache_file:
        cache = json.load(cache_file)

    # Read data
    galaxy = cache[galaxy_name]

    # Organize the data by cataglog
    col1, col2 = columns
    output = {}
    for cat in catalogs:
        if cat not in cache[META_DATA_KEY]["Included titles"].keys():
            print(f"{cat} not found in galaxy, skipping")
            continue

        col1, col2 = [], []
        col1_name, col2_name = columns
        output[cat] = {col1_name: col1, col2_name: col2}

        for star_name in galaxy.keys():
            star = galaxy[star_name]
            if cat not in star.keys():
                continue
            if col1_name not in star[cat].keys():
                continue
            if col2_name not in star[cat].keys():
                continue

            output[cat][col1_name].append(star[cat][col1_name])
            output[cat][col2_name].append(star[cat][col2_name])

    return output

def plot_data(data:dict, title:str = None, savepath:str=None, plot_size = (5, 5)):
    """Plots the data extracted from the cache file. Can be applied to the 
    output of either data_by_galaxy or data_by_catalog.

    Args:
        data (dict): The data to plot.
        title (str): The title of the plot.
        savepath (str): The path to save the plot to. If None, the plot is
            displayed instead of saved.

    Returns:
        None (plots the data and saves it)
    """
    marker_styles = ["o", "v", "^", "<", ">", "s", "p", "P", "*", "h", "H"]
    fig, ax = plt.subplots(1, 1, figsize=plot_size)

    count = 0
    for plot_num in data.keys():
        cols = data[plot_num].keys()
        col1_name, col2_name = cols

        col1 = data[plot_num][col1_name]
        col2 = data[plot_num][col2_name]

        ax.scatter(col1, col2, label = plot_num, s = 15,
                   marker = marker_styles[count], alpha=0.8)

        count += 1

    ax.set_xlabel(f'[{col1_name}]')
    ax.set_ylabel(f'[{col2_name}]')
    ax.set_title(title)
    ax.legend()

    if savepath is not None:
        fig.savefig(savepath)
        return

    plt.show()


if __name__ == "__main__":
    # data_galaxy = data_by_galaxy("Data/cache.json", "J/ApJS/191/352/abun",
                        #   ["Scl", "For", "TriII"], ["Fe_H", "Mg_Fe"])
    # data_catalog = data_by_catalog("Data/cache.json", "Scl",
                    # ["J/ApJS/191/352/abun", "J/ApJ/838/83", "J/A+A/641/A127"],
                    # ["Fe_H", "Mg_Fe"])

    # plot_data(data_galaxy, "Test", plot_size=(8, 8))

    Kirby_Mg_abun = data_by_galaxy("Data/cache.json", "J/ApJS/191/352/abun",
                          [ "Scl", "For", "LeoI", "Sex", "LeoII", "CVnI", "UMi", "Dra"], ["Fe_H", "Mg_Fe"])
    plot_data(Kirby_Mg_abun, savepath = "Output/Kirby 2009 Mg Abundance",
              plot_size=(8, 8))
    
    # catalog_cross_ref = data_by_catalog
