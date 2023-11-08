"""Debug file, for testing purposes."""
from matplotlib import pyplot as plt
import numpy as np
import dataloader as dl

def catalog_statistics(prompt):
    """Prints the number of stars and galaxies in one catalog, given a Vizier prompt."""
    catalog = dl.retreive_catalog(prompt)
    catalog_name = list(catalog.keys())[0]
    catalog = catalog[catalog_name]

    count = 0
    for star in catalog:
        count += 1
    star_count = str(count) + " stars"

    galaxies = []
    for star in catalog:
        galaxies.append(star['dSph'])
    galaxies = set(galaxies)
    gal_count = str(len(galaxies)) + " galaxies"

    print(catalog_name, "has", star_count, "in", gal_count)

catalog_statistics("2011yCat..21910352K")

# dl.retreive_catalog("2011yCat..21910352K")
