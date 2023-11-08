"""Makes position plot of """
from matplotlib import pyplot as plt
from makeplots import extract_columns

filepath = "Temp/Kirby 2009.json"
catalog_name = "J/ApJS/191/352/abun"

if __name__ == "__main__":
    # Make position plot
    loc = extract_columns(filepath, 
                              catalog_name, 
                              ["RAJ2000", "DEJ2000"])

    fig, (pos, gal1, gal2, gal3) = plt.subplots(1, 4, figsize = (15, 5))

    pos.scatter(loc["RAJ2000"], loc["DEJ2000"], s = 1, color = 'black')

    fig.supxlabel("RA (deg)")
    fig.supylabel("Dec (deg)")
    pos.set_title("Spatial Distribution of all Stars in Kirby 2009")
    pos.grid()

    gal1.scatter(loc["RAJ2000"], loc["DEJ2000"], s = 1, color = 'red')
    gal1.set_title("Lower Left Section")
    gal1.set_xlim(0.5, 3)
    gal1.set_ylim(-34, -32)
    gal1.grid()

    gal2.scatter(loc["RAJ2000"], loc["DEJ2000"], s = 1, color = 'blue')
    gal2.set_title("Middle Section")
    gal2.set_xlim(10, 13.55)
    gal2.set_ylim(-2, 36)
    gal2.grid()

    gal3.scatter(loc["RAJ2000"], loc["DEJ2000"], s = 1, color = 'green')
    gal3.set_title("Upper Right Section")
    gal3.set_xlim(15, 17.5)
    gal3.set_ylim(57, 68)
    gal3.grid()

    fig.suptitle("Kirby 2009 Position Plot, total stars: " + str(len(loc["RAJ2000"])))
    fig.tight_layout()

    fig.savefig("Output/Kirby 2009 Position Plot.png")