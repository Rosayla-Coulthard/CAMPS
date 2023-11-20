"""Debug file, for testing purposes."""
from matplotlib import pyplot as plt
import numpy as np
from config import TEMP_FOLDER
import json

catalog_path = f"{TEMP_FOLDER}/Kirby_2009.json"
catalog = json.load(open(catalog_path, encoding="utf-8"))

# catalog = [catalog['Scl'][i]['J/ApJS/191/352/abun'] for i in list(catalog['Scl'].keys())]
star_coords = []

# Find the max and min RA and DEC
for i in list(catalog['Scl'].keys()):
    star = catalog['Scl'][i]['J/ApJS/191/352/abun']
    RA_temp = np.float64(star['RAJ2000'])
    DEC_temp = np.float64(star['DEJ2000'])

    star_coords.append([RA_temp, DEC_temp])

RA, DEC = np.transpose(star_coords)

RA_bound = (max(RA) - min(RA)) * 0.5
DEC_bound = (max(DEC) - min(DEC)) * 0.5

angle = np.linspace(0, 2*np.pi, 100)

bound_circle_x = 1.25 * RA_bound * np.cos(angle) + (max(RA) + min(RA)) * 0.5
bound_circle_y = 1.25 * DEC_bound * np.sin(angle) + (max(DEC) + min(DEC)) * 0.5
# Outside of this bound, start a new galaxy

plt.plot(bound_circle_x, bound_circle_y)
plt.scatter(RA, DEC, s=0.1)
# plt.Circle(np.average(RA), np.average(DEC), 0.5, color='r')
plt.show()
