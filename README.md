# Catalog-Maker

## How to use
1. Retreive data from Vizier:
Run `retrieve_catalog()` following the instructions in the docstring and the prompts in the terminal. Save the output as a separate .json file. If `config.json` has no record of a column header contained in the retrieved table, it will prompt you to rename the column as it will appear in the cache file. If the column header is already in `config.json`, it will automatically rename the column as it appears in the cache file. Once all new column headers have been renamed and added, the script will terminate and you will need to run it again. 
2. Merge tables:
If the Vizier prompt retrieves multiple tables, save them as separate `.json` files and run `merge_tables()` to merge them into one table.
3. Create catalog compilation:
Copy one of the existing `.json` catalog files and make it the cache file. Run `add_catalog()` to add one of the already retrieved catalogs to the cache file.

