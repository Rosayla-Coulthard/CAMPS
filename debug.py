from astropy.table import Table

# Make an astropy table row by row
t = Table(names=('a', 'b', 'c'), dtype=('i4', 'f4', 'S2'))
for i in range(10):
    t.add_row((i, i**2, 'x' * i))
    print(t)

# Print the table
print(t)

