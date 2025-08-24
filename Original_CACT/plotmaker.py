"""Makes plots."""
import json
import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

from utilfuncs import catalog_name_check, META_DATA_KEY

def linear(x, m, b):
    """A linear function."""
    return m*x + b

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
        included_titles = cache[META_DATA_KEY]["Included titles"]
        cat_list = list(included_titles.keys())
        cat_col_list = list(included_titles[cat]["Columns"])
        cat_mem_list = list(included_titles[cat]["Members"])

        # print(included_titles)
        # print(cat_list)
        # print(cat_col_list)
        # print(cat_mem_list)

        # quit()

        if cat not in cat_list:
            print(f"{cat} not found in cache, skipping")
            continue
        elif col1 not in cat_col_list:
            print(f"{col1} not found in {cat}, skipping")
            continue
        elif col2 not in cat_col_list:
            print(f"{col2} not found in {cat}, skipping")
            continue
        elif galaxy_name not in cat_mem_list:
            print(f"{galaxy_name} not found in {cat}, skipping")
            continue

        col1_list, col2_list = [], []
        col1_name, col2_name = columns
        output[cat] = {col1_name: col1_list, col2_name: col2_list}

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

def data_for_crossmatch(cache_path:str, catalogs:list, column:str, gal_list=None):
    """Extracts data for crossmatching between catalogs. Only extracts data from 
    the same star.
    
    Args:
        cache_path (str): The path to the cache file.
        catalogs (list): A list of catalog names.
        column (str): The column to extract.
        gal_list (list): A list of galaxies to extract data from. If None, all
            galaxies are used.

    Returns:
        output (dict): A dictionary of the data, organized by catalog.
    """
    # Check inputs
    if len(catalogs) != 2:
        raise ValueError("catalogs must have exactly two elements")

    # Load cache
    with open(cache_path, "r", encoding="utf-8") as cache_file:
        cache = json.load(cache_file)

    # Read data
    cat1, cat2 = catalogs
    cat1_data, cat2_data = [], []

    if gal_list is None:
        gal_list = list(cache.keys())[1:]

    # Check if the requested column is in both catalogs
    included_titles = cache[META_DATA_KEY]["Included titles"]
    if column not in included_titles[cat1]["Columns"]:
        raise ValueError(f"{column} not found in {cat1}")
    if column not in included_titles[cat2]["Columns"]:
        raise ValueError(f"{column} not found in {cat2}")

    # Scan through the cache
    for gal_name in gal_list:
        gal = cache[gal_name]
        for star_name in gal.keys():
            star = gal[star_name]

            try:
                temp_1 = star[cat1][column]
                temp_2 = star[cat2][column]
            except KeyError:
                continue

            cat1_data.append(temp_1)
            cat2_data.append(temp_2)

    output = {column:{cat1: cat1_data, cat2: cat2_data}}
    return output


def plot_data(data:dict, title:str = None, savepath:str=None,
              plot_size = (5, 5), nbins = 0, units_xy = ("", "")):
    """Plots the data extracted from the cache file. Can be applied to the 
    output of either data_by_galaxy or data_by_catalog.

    Args:
        data (dict): The data to plot.
        title (str): The title of the plot.
        savepath (str): The path to save the plot to. If None, the plot is shown 
            instead of saved.
        plot_size (tuple): The size of the plot.
        nbins (int): The number of bins to use for the error bars. If 0, no error
            bars are plotted.
        units_xy (tuple): The units of the x and y axes.

    Returns:
        None (Plots the data and saves the plot if savepath is not None.)
    """
    marker_styles = ["o", "v", "^", "<", ">", "s", "p", "P", "*", "h", "H"]
    fig, ax = plt.subplots(1, 1, figsize=plot_size)

    count = 0
    for plot_num in data.keys():
        cols = data[plot_num].keys()
        col1_name, col2_name = cols

        col1 = data[plot_num][col1_name]
        col2 = data[plot_num][col2_name]

        # ax.scatter(col1, col2, label = plot_num, s = 15,
        #            marker = marker_styles[count], alpha=0.5)

        if nbins > 0:
            hist, edges = np.histogram(col1, bins = nbins)
            middle = (edges[1:] + edges[:-1]) / 2

            col_min = min(col1)
            bins = {str(i):[] for i in range(nbins)}
            
            for i, point in enumerate(col1):
                if math.isnan(col2[i]) is False:
                    # Calculate the index of the bin to put the point in
                    bin_indx = min(math.floor(nbins * (point / abs(col_min) + 1)), nbins - 1)
                    bins[str(bin_indx)].append(col2[i])

            bin_means = []
            bin_stds = []
            for bin_temp in bins.values():
                if len(bin_temp) != 0:
                    bin_means.append(np.mean(bin_temp))
                    bin_stds.append(np.std(bin_temp))
                else:
                    bin_means.append(np.nan)
                    bin_stds.append(np.nan)

            ax.errorbar(middle, bin_means, yerr = bin_stds, alpha = 0.6,
                        marker = marker_styles[count], label = plot_num, markersize = 10, capsize=3)
            ax.scatter(col1, col2, s = 15,
                       marker = marker_styles[count], alpha=0.1)

        else:
            ax.scatter(col1, col2, s = 15,
                       marker = marker_styles[count], alpha=0.5)

        count += 2
        print(f"Plotted {plot_num}")

    units_x, units_y = units_xy
    
    if units_x != "":
        units_x = f" [{units_xy[0]}]"

    if units_y != "":
        units_y = f" [{units_xy[1]}]"


    ax.set_xlabel(f'[{col1_name}] {units_x}')
    ax.set_ylabel(f'[{col2_name}] {units_y}')
    ax.set_title(title)
    ax.legend(ncol = min(int((count+1)/2), 5), bbox_to_anchor=(0, 1),
              loc='lower left')
    ax.grid()

    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath)
        return

    plt.show()

def plot_crossmatch(data:dict, title:str = None, savepath:str=None, plot_size = (5, 5),
                    units_xy = ("", ""), col_title = ''):
    """Plots the data extracted from the cache file. Can be applied to the 
    output of either data_by_galaxy or data_by_catalog.

    Args:
        data (dict): The data to plot.
        title (str): The title of the plot.
        savepath (str): The path to save the plot to. If None, the plot is shown 
            instead of saved.
        plot_size (tuple): The size of the plot.
        units_xy (tuple): The units of the x and y axes.
        col_title (str): The title of the column being plotted.

    Returns:
        None (Plots the data and saves the plot if savepath is not None.)
    """
    # Unpack data
    column = list(data.keys())[0]
    data = data[column]
    cat1, cat2 = list(data.keys())
    cat1_data, cat2_data = data[cat1], data[cat2]

    # Fit data to a linear function
    popt, pcov = curve_fit(linear, cat1_data, cat2_data)
    m, b = popt
    m_err, b_err = np.sqrt(np.diag(pcov))

    fit_curve = linear(np.array(cat1_data), m, b)

    # Plot data
    fig, (ax, residual) = plt.subplots(2, 1, figsize=plot_size)

    ax.scatter(cat1_data, cat2_data, label = f"[{column}] between {cat1} and {cat2}",
               s = 15, marker = "o", alpha=0.5)
    plot_label = f"Linear fit: y = ({m:.3f} +/- {round(m_err, 3)})x + ({b:.3f} +/- {round(b_err, 3)})"
    ax.plot(cat1_data, fit_curve, label = plot_label,
            color = "red", alpha = 0.5)
    ax.plot(cat1_data, linear(np.array(cat1_data), 1, 0), label = "y = x", linestyle = "--",
            color = "black", alpha = 0.5)

    fig.supxlabel(f'[{column}] {cat1} {units_xy[0]}')
    ax.set_ylabel(f'[{column}] {cat2} {units_xy[1]}')
    ax.set_title(title)
    ax.legend(ncol = 1, bbox_to_anchor=(0, 1), loc='lower left')
    ax.grid()

    resi = np.array(cat2_data) - fit_curve
    print(f"Residual mean: {np.mean(resi)}, std: {np.std(resi)}")

    residual.scatter(cat1_data, resi, label = "Residuals",
                        s = 15, marker = "o", alpha=0.5)
    residual.axhline(0, color = "red", alpha = 0.5, label = 'linear fit')
    residual.set_ylabel('Residuals')
    residual.legend(ncol = 2, bbox_to_anchor=(0, 1), loc='lower left')

    fig.tight_layout()

    fitting_params = np.array([m, m_err, b, b_err])
    fit_data = np.array([np.array(cat1_data), np.array(cat2_data)])

    if savepath is not None:
        fig.savefig(savepath)
        return fitting_params, fit_data, resi

    plt.show()

    return fitting_params, fit_data, resi

if __name__ == "__main__":
    gal_list = ["Scl", "For", "LeoI"]#, "Sex", "LeoII", "CVnI", "UMi", "Dra"]
    Kirby_Mg_abun = data_by_galaxy("Data/cache.json", "Kirby 2009",
                                    gal_list , ["Fe_H", "Mg_Fe"])
    plot_data(Kirby_Mg_abun, savepath = "Plots/Kirby 2009 Mg Abundance",
              nbins = 10, units_xy = ("log", "log"))

    # catalog_list = ["Kirby 2009", "Reichert 2020", "Theler 2020", "Skuladottir 2019"]
    # catalog_cross_ref_scl = data_by_catalog("Data/cache.json", "Scl",
                                            # catalog_list, ["Fe_H", "Mg_Fe"])
    # plot_data(catalog_cross_ref_scl, savepath = "Output/Scl Cross Reference",
            #   units_xy=("log", "log"))

    # catalog_cross_ref_sex = data_by_catalog("Data/cache.json", "Sex",
                                            # catalog_list, ["Fe_H", "Mg_Fe"])
    # plot_data(catalog_cross_ref_sex, savepath = "Output/Sex Cross Reference",
            #   units_xy=("log", "log"))
    
    # calibration_data = data_for_crossmatch("Data/cache.json",
                                            #  ["Reichert 2020", "Kirby 2009"],
                                            #  "Fe_H")
    # plot_crossmatch(calibration_data, savepath = "Output/Fe_H Calibration")
