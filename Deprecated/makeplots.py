"""Plot Spatial distribution, Color-Magnitude Diagram, 
Metallicity Distribution Function, and [X/Fe] vs. [Fe/H] distributions for each galaxy
"""
# NOTE: but first, separate the galaxies
import json
import numpy as np
import matplotlib.pyplot as plt

import config as c

def extract_columns(filepath, galaxy = None, catalog = None, columns = None):
    """Extracts specified data from the catalog.
    Three modes for each parameter:
        None: Extracts all data
        str: Extracts only specified data
        list: Extracts only specified data
    """
    # Load data
    data_file = json.load(open(filepath, "r", encoding = "utf-8"))

    # If any parameter is a string, turn that string into a single-element list
    # If any paramter is a list, use that list
    # If any parameter is None, turn that parameter into a list of all keys
    if isinstance(galaxy, str):
        galaxy = [galaxy]
    elif galaxy is None:
        galaxy = list(data_file.keys())

    if not isinstance(catalog, str):
        raise ValueError("Must specify catalog(s) to extract from.")

    if isinstance(columns, str):
        columns = [columns]
    elif columns is None:
        star0 = list(data_file[galaxy[0]].keys())[0]
        print(star0)
        columns = list(data_file[galaxy[0]][star0][catalog].keys())

    # Create output dictionary
    output_dict = {}

    for gal in galaxy:
        output_dict[gal] = {col: [] for col in columns}
        star_list = data_file[gal].keys()
        for star in star_list:
            star_data = data_file[gal][star][catalog]
            for col in columns:
                try:
                    output_dict[gal][col].append(float(star_data[col]))
                except ValueError:
                    output_dict[gal][col].append(np.nan)
        
        for col in columns:
            output_dict[gal][col] = np.array(output_dict[gal][col])

    # print(output_dict)
    return output_dict

def remove_nan_for_plot(x, xerr, y, yerr):
    """Cleans up the nan inputs for plotting later"""
    nan_mask = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(xerr) & ~np.isnan(yerr)
    x = x[nan_mask]
    y = y[nan_mask]
    xerr = xerr[nan_mask]
    yerr = yerr[nan_mask]

    return x, xerr, y, yerr

if __name__ == "__main__":
    filepath = "Temp/Kirby 2009.json"
    catalog_name = "J/ApJS/191/352/abun"

    data_column = ["__Fe_H_", "e__Fe_H_",
                "__Mg_Fe_", "e__Mg_Fe_",
                "__Si_Fe_", "e__Si_Fe_",
                "__Ca_Fe_", "e__Ca_Fe_"]

    # Extract the relevant data
    data = extract_columns(filepath, 
                        galaxy = None, 
                        catalog = catalog_name, 
                        columns = data_column)

    fig, ax = plt.subplots(3, 1, figsize = (9, 7))

    MgFe, CaFe, SiFe = ax

    FEH_LABEL = data_column[0]
    fig.supxlabel(c.ALT_COL_NAMES[FEH_LABEL])

    for gal in data.items():
        gal_name, gal_data = gal
        FeH, FeH_err = gal_data[data_column[0]], gal_data[data_column[1]]

        for i, _ in enumerate(data_column[2::2]):
            abun, abun_err = data_column[2 + i*2], data_column[3 + i*2]

            abun_label = c.ALT_COL_NAMES[abun]
            abun, abun_err = gal_data[abun], gal_data[abun_err]

            FeH_temp, FeH_err_temp, abun, abun_err = remove_nan_for_plot(FeH, FeH_err, abun, abun_err)

            ax[i].errorbar(FeH_temp, abun,
                        xerr = FeH_err_temp, yerr = abun_err,
                        fmt = '.', label = gal_name,
                        color = 'black', alpha = 0.4)
            ax[i].set_ylabel(abun_label)
            ax[i].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        break


    fig.suptitle("Kirby 2009 Chemical Abundance Plot")

    fig.tight_layout()
    fig.savefig(f"{c.OUTPUT_FOLDER}/Kirby 2009 Chem Abun Plot.png")


