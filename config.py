"""Configuration file."""

DATA_FOLDER = "Data"
TEMP_FOLDER = "Temp"
OUTPUT_FOLDER = "Output"

# The columns corresponding to values we concern
COLUMNS=[
        "RAJ2000",
        "DEJ2000",
        "dSph",
        "Name",
        "Vmag",
        "Imag",
        "Teff",
        "e_Teff",
        "logg",
        "e_logg",
        "[Fe/H]",
        "e_[Fe/H]",
        "[Mg/Fe]",
        "e_[Mg/Fe]",
        "[Si/Fe]",
        "e_[Si/Fe]",
        "[Ca/Fe]",
        "e_[Ca/Fe]",
        "[Ti/Fe]",
        "e_[Ti/Fe]",
    ]

# Reference list of alternative names
ALT_COL_NAMES = {
    '__Fe_H_': '[Fe/H]',
    'e__Fe_H_': 'e[Fe/H]',
    '__Mg_Fe_': '[Mg/Fe]',
    'e__Mg_Fe_': 'e[Mg/Fe]',
    '__Si_Fe_': '[Si/Fe]',
    'e__Si_Fe_': 'e[Si/Fe]',
    '__Ca_Fe_': '[Ca/Fe]',
    'e__Ca_Fe_': 'e[Ca/Fe]'
}