"""Plots chemical abundance and metallicity distribution of galaxies"""
import numpy as np
import matplotlib.pyplot as plt
from makeplots import extract_columns

column_list = ["__Fe_H_",
               "e__Fe_H_",
               "__Mg_Fe_",
               "e__Mg_Fe_",
               "__Si_Fe_",
               "e__Si_Fe_",
               "__Ca_Fe_",
               "e__Ca_Fe_",
               "__Ti_Fe_",
               "e__Ti_Fe_"]

# Extract data
data = extract_columns("Temp/Kirby 2009.json",
                       "J/ApJS/191/352/abun",
                       ["dSph"])
                       # column_list)

print(data)
