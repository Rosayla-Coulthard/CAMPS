import json
import pandas as pd
from astropy import table 
from astropy.table import Table, Column, MaskedColumn, vstack, unique, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from IPython.display import clear_output
from astropy.io import ascii
import numpy as np
from astroquery.vizier import Vizier
import collections
from collections import Counter
import sqlutilpy
import utilfuncs as uf
import dataloader as dl
import sys
sys.path.append('/home/jupyter-tingli/pythoncode/')
import crossmatcher
import warnings
warnings.filterwarnings('ignore')
config = json.load(open("config.json", encoding="utf-8"))
renames = config["RENAMES"]
gals_list = config["GALS_LIST"]
paper_renames = config["PAPER_RENAMES"]
objects = config["OBJECTS"]
saga_corrs = config["SAGA_CORRS"]
elem_order = config["ELEM_ORDER"]

def fix_types(lit_cat:Table):
    potential = ['Gaia_ID','Galaxy','Object Type','Paper','Bib','Star_ID', 'Name', 'Object Name', 'From', 'Member', 'Facility', 'Resolution', 'Solar Abundance']
    for col in potential:
        if col in lit_cat.columns:
            lit_cat[col] = np.array(lit_cat[col].data, dtype=str)
    for col in lit_cat.columns:
        if col not in potential:
            try:
                lit_cat[col] = np.array(lit_cat[col], dtype=float)
            except:
                lst2=[]
                for entry in lit_cat[col]:
                    if isinstance(entry,float):
                        lst2.append(entry)
                    elif entry is np.ma.masked or entry == '' or entry.isalpha():
                        lst2.append(np.nan)
                    else:
                        lst2.append(float(entry))
                lit_cat[col] = lst2

def retrieve_catalog_bulk(paper_tag, tables, paper_name, star_name, paper_type,
                          gal_title, gal_name = None, save_path = None):
    """Retrieves all tables from the same paper"""
    # Retrieve all tables
    path_list = []
    for table in tables:
        if len(tables) <= 1:
            temp_path = save_path
        else:
            temp_path = save_path.split(".")
            temp_path = f"{temp_path[0]}({table}).{temp_path[1]}"

        path_list.append(temp_path)

        if gal_title is None:
            IsSingleGalaxy = True
        else:
            IsSingleGalaxy = False
            
        output = dl.retrieve_catalog(paper_type, paper_tag, temp_path, UI=star_name, gal_title=gal_title,
                            table_num=table, paper_name=paper_name,
                            single_galaxy=IsSingleGalaxy, gal_name = gal_name)
        uf.save_to_json(output, temp_path)


    # Merge tables if there are more than one
    if len(path_list) > 1:
        for count, path in enumerate(path_list):
            if count == 0:
                first_table_path = path
            else:
                output = dl.merge_tables(first_table_path, path, paper_name, paper_name,
                                save_path)

    uf.save_to_json(output, save_path)

    # Save the output. If there is only one table, the output was determined earlier
    return output

def stack(lit_cat: Table, paper_list:dict, paper:str, entry_count: int):
    clear_output(wait=True)
    print('Merging ', paper, ' . . .')
    load = json.load(open(paper_list[paper]['path'], "r"))
    incoming = Table.from_pandas(uf.make_dataframe(load, False))
    print(incoming['Paper'])
    
    #Giving the incoming entries a unique star ID, then accounting for duplicates within the paper itself
    incoming['Star_ID'] = Column([f"{i:05.0f}" for i in range(entry_count + 1, len(incoming) + entry_count+1)])
    dupes = ([item for item, count in collections.Counter(incoming['RAJ2000']).items() if count > 1])
    swap_tuples = []
    for duplicate in dupes:
        indices = [i for i, x in enumerate(incoming['RAJ2000']) if x == duplicate]
        swap_tuples.append(indices)
    #swap_tuples is a list of indices that have the same star name. They need to be changes to have the same star_ID name
    for duplicate in swap_tuples:
        for i in range(1,len(duplicate)-1):
            incoming['Star_ID'][i] = incoming['Star_ID'][duplicate[0]]
    #slightly modify the star ID names so that there each star can be identified as unique
    entry_count = entry_count+len(incoming)
    
    c = SkyCoord(ra=incoming['RAJ2000']*u.degree, dec=incoming['DEJ2000']*u.degree)
    catalog = SkyCoord(ra=lit_cat['RAJ2000']*u.degree, dec=lit_cat['DEJ2000']*u.degree)

    match, d2d, d3d = c.match_to_catalog_sky(catalog)

    
    lit_inIncoming = lit_cat[match][d2d < 1*u.arcsec]
    Incoming_inlit = incoming[d2d < 1*u.arcsec]
    Incoming_notinlit = incoming[d2d >= 1*u.arcsec]
    
    
    #for stars in incoming, already in the catalogue change their name to what the catalogue has
    key = list(zip(lit_inIncoming['Star_ID'], Incoming_inlit['Star_ID']))
    key2 = list(zip(lit_inIncoming['Galaxy'], Incoming_inlit['Galaxy']))
    
    for idx in range(len(key)):
        incoming['Star_ID'][incoming['Star_ID']==key[idx][1]] = key[idx][0]
        incoming['Galaxy'][incoming['Galaxy']==key2[idx][1]] = key2[idx][0]
    

    incoming_mod = incoming['Star_ID'].tolist()        

    all_cols = set(lit_cat.columns) | set(incoming.columns)
    # Add missing columns to table1
    for col in all_cols:
        if col not in lit_cat.columns:
            lit_cat[col] = ['nan'] * len(lit_cat)
        if col not in incoming.columns:
            incoming[col] = ['nan'] * len(incoming)

    # Ensure columns are in the same order
    lit_cat = lit_cat[sorted(all_cols)]
    incoming = incoming[sorted(all_cols)]
    lit_cols = list(lit_cat.columns)
    
    for col in lit_cols:
        if col != 'RAJ2000' and col != 'DEJ2000':
            lit_cat[col] = lit_cat[col].astype(str)
            incoming[col] = incoming[col].astype(str)
    
    print(paper, ' merged!\n')
    return (vstack([lit_cat, incoming]),paper_list, paper, entry_count)

def galaxies(lit_cat: Table):
    curr_gals = lit_cat['Galaxy'].tolist()
    print(f"Found the following galaxies: {np.unique(curr_gals)}")
    new = []
    for gal in curr_gals:
        if gal not in renames:
            msg = f"Found galaxy '{gal}' not in convention file. Please specify a name (not shortened): "
            user_input1 = input(msg)
            if user_input1 == "":
                print("WARNING: Input empty. Operation aborted, no changes made.")
            msg = f"Please specify a short form (eg: Aquarius --> Aqr): "
            user_input2 = input(msg)
            if user_input2 == "":
                print("WARNING: Input empty. Operation aborted, no changes made.")
            msg = f"Please specify an object type (eg: UFD, GC, dIrr): "
            user_input3 = input(msg)
            if user_input3 == "":
                print("WARNING: Input empty. Operation aborted, no changes made.")
            config["RENAMES"][gal] = user_input1
            config["GALS_LIST"][user_input1] = user_input2
            config["OBJECTS"][user_input1] = user_input3
            uf.save_to_json(config, "config.json")
        
        new.append(renames[gal])
    lit_cat['Galaxy'] = new
    return lit_cat

def clean(lit_cat:Table):
    key4 = Table()
    key4['old'] = list(objects.keys())
    key4['new'] = list(objects.values())

    mapping4 = {row['old']: row['new'] for row in key4}
    new_names = []

    for old_name in lit_cat['Galaxy']:
        new_names.append(mapping4.get(old_name, None))
    lit_cat['Object Type'] = np.array(new_names, dtype=str)
    lit_cat.rename_column('Galaxy', 'Object Name')
    
    ordered = ['Gaia_ID','Star_ID','Name','Paper', 'Bib', 'From','RAJ2000','DEJ2000','Object Name', 'Object Type','log_g','er_log_g','T_eff',
                          'er_T_eff','RV','er_RV','V_mic','er_V_mic','Member','S_N','Facility', 'Resolution', 'Solar Abundance','Fe_H','er_Fe_H','FeII_H','er_FeII_H','FeI_H','er_FeI_H', 'M_H', 'A_C', 'A_Na', 'A_Mg', 'A_Ca', 'A_Ti', 'A_Fe', 'A_FeI', 'A_FeII', 'A_Sr', 'A_Ba']
    ordered2 = elem_order
    cols = list(lit_cat.columns)
    
    one = []
    two = []
    three = []
    for element in ordered:
        if element in cols:
            one.append(element)
    for element in ordered2:
        if element in cols:
            three.append(element)
    for element in cols:
        if element not in ordered:
            if element not in ordered2:
                two.append(element)
    one.extend(three)
    for col in two:
        lit_cat.remove_column(col)
    return(lit_cat[one])

def collapse(lit_cat: Table, Add_SAGA):
    names = ["H_Fe","He_Fe","Li_Fe","Be_Fe","B_Fe","C_Fe","N_Fe","O_Fe","F_Fe","Ne_Fe","Na_Fe","Mg_Fe",
                "Al_Fe","Si_Fe","P_Fe","S_Fe","Cl_Fe","Ar_Fe","K_Fe","Ca_Fe","Sc_Fe","Ti_Fe","V_Fe","Cr_Fe","Mn_Fe",
                "Co_Fe","Ni_Fe","Cu_Fe","Zn_Fe","Ga_Fe","Ge_Fe","As_Fe","Se_Fe","Br_Fe","Kr_Fe","Rb_Fe","Sr_Fe","Y_Fe",
                "Zr_Fe","Nb_Fe","Mo_Fe","Tc_Fe","Ru_Fe","Rh_Fe","Pd_Fe","Ag_Fe","Cd_Fe","In_Fe","Sn_Fe","Sb_Fe",
                "Te_Fe","I_Fe","Xe_Fe","Cs_Fe","Ba_Fe","La_Fe","Ce_Fe","Pr_Fe","Nd_Fe","Pm_Fe","Sm_Fe","Eu_Fe","Gd_Fe",
                "Tb_Fe","Dy_Fe","Ho_Fe","Er_Fe","Tm_Fe","Yb_Fe","Lu_Fe","Hf_Fe","Ta_Fe","W_Fe","Re_Fe","Os_Fe","Ir_Fe",
                "Pt_Fe","Au_Fe","Hg_Fe","Tl_Fe","Pb_Fe","Bi_Fe","Po_Fe","At_Fe","Rn_Fe","Fr_Fe","Ra_Fe","Ac_Fe","Th_Fe",
                "Pa_Fe","U_Fe","Np_Fe","Pu_Fe","Am_Fe","Cm_Fe","Bk_Fe","Cf_Fe","Es_Fe","Fm_Fe","Md_Fe","No_Fe","Lr_Fe",
                "Rf_Fe","Db_Fe","Sg_Fe","Bh_Fe","Hs_Fe","Mt_Fe","Ds_Fe","Rg_Fe","Cn_Fe","Nh_Fe","Fl_Fe","Mc_Fe","Lv_Fe",
                "Ts_Fe","Og_Fe","C2_Fe","CH_Fe","CN_Fe","HF_Fe","NH_Fe","OH_Fe", 'Fe_H']
    if Add_SAGA:
        lit_cat.rename_column('BI_FeI', 'B_Fe')
        lit_cat.rename_column('TcI_FeI', 'Tc_Fe')
        
    for name in names:
        er = 'er_' + name
        if name == 'Fe_H':
            var1 = 'FeI_H'
            var2 = 'FeII_H'
            er1 = 'er_FeI_H'
            er2 = 'er_FeII_H'
        else:
            split1 = name.split('_')
            var1 = split1[0] + 'I_' + split1[1] + 'I'
            er1 = 'er_' + split1[0] + 'I_' + split1[1] + 'I'
            var2 = split1[0] + 'II_' + split1[1] + 'II'
            er2 = 'er_' + split1[0] + 'II_' + split1[1] + 'II'
        try:
            mask1 = np.isnan(lit_cat[name])
            lit_cat[name][mask1] = lit_cat[var1][mask1] 
            lit_cat.remove_column(var1)
        except:
            pass
        try:
            mask1 = np.isnan(lit_cat[name])
            lit_cat[name][mask1] = lit_cat[var2][mask1]
            lit_cat.remove_column(var2)
        except:
            pass
        try:
            mask2 = np.isnan(lit_cat[er])
            lit_cat[er][mask2] = lit_cat[er1][mask2]
            lit_cat.remove_column(er1)
        except:
            pass
        try:
            mask2 = np.isnan(lit_cat[er])
            lit_cat[er][mask2] = lit_cat[er2][mask2]
            lit_cat.remove_column(er2)
        except:
            pass

def fix_names(lit_cat: Table):
    # reduces the table to only have unique values
    curr_gals = np.unique(lit_cat['Galaxy']).tolist()
    # all of the current galaxies
    curr_name = []
    new_name = []

    for gal in curr_gals:
        mask = lit_cat['Galaxy'] == gal
        clear_output(wait=True)
        print('Sorting ', gal, '...')
        lit_cat_gal = lit_cat[mask]
        form = gals_list[gal]
        # only examining one table at a time

        names = list(set(lit_cat_gal['Star_ID'].tolist()))
        corresponding = [f"{form}_{i:05.0f}" for i in range(1, len(names) + 1)]
        # names is the list of unique star names (accounting for possible duplicate paper names)
        # corresponding is the list of what their names should be
        curr_name = curr_name + names
        new_name = new_name + corresponding
        #for idx in range(len(names)):
         #   lit_cat['Star_ID'][lit_cat['Star_ID']== names[idx]]= corresponding[idx]
          #  print(names[idx], corresponding[idx])
            # Change each instance of the name in the real catalogue,
        print(gal +' complete')

    name_map = dict([(curr_name[i], new_name[i]) for i in range(0, len(curr_name))])
    new_names = [name_map[old] for old in lit_cat['Star_ID']]
    lit_cat['Star_ID'] = new_names
    lit_cat = lit_cat[np.argsort(lit_cat['Star_ID'])]
    clear_output(wait=True)
    
    return lit_cat

def lit_bibs(lit_cat, tag_coll):
    lit_cat['From'] = np.full(len(lit_cat), 'Manual Compilation')

    lit_bibs = []
    papers = []
    for i in range(len(tag_coll)):
        papers.append(tag_coll[i][0])
    for bib in tag_coll:
        if '.' in bib[1]:
            lit_bibs.append(bib[1])
            continue
        else:
            split1 = bib[1].split('/')
            author = bib[0].split(' ')[-2]
            year = bib[0].split(' ')[-1]
            year = year.split('-')[0]
            journal = split1[1].strip()
            if journal == 'A+A':
                journal = 'A&A'
            volume = split1[2].strip()
            page = split1[3].strip()

            if not volume[-1].isnumeric():
                p = volume[-1]
                volume = volume[:-1]
            else:
                p = '.'
            if not (page[0]).isnumeric():
                p = page[0]
                page = page[1:]
            else:
                p = '.'
            lit_bibs.append(f"{year}{journal:{'.'}<5}{volume:{'.'}>4}{p}{page:{'.'}>4}{author[0]}")
    adds = []
    for entry in lit_cat['Paper']:
        adds.append(lit_bibs[papers.index(entry)])
    lit_cat['Bib'] = adds
    for i in range(len(lit_cat)):
        if lit_cat['Paper'][i] == 'Yoon 2016':
            lit_cat['Paper'][i] == 'nan'
            lit_cat['Bib'][i] == 'nan'
            lit_cat['From'][i] = 'Yoon 2016'
    masky = (lit_cat['Paper'] == 'Yoon 2016')
    lit_cat['Paper'][masky] = lit_cat['Paper_Y'][masky] 
    lit_cat['Bib'][masky] = lit_cat['Bib_Y'][masky] 
    lit_cat.remove_column('Bib_Y')
    lit_cat.remove_column('Paper_Y')

            
def saga(lit_cat: Table, curr_count: int):
    SAGAv2 = ascii.read('/home/jupyter-tingli/Tutorials/saga_cleaned_catalog.tsv')
    SAGAv2.rename_column('RAdeg', 'RA2000')
    SAGAv2.rename_column('DECdeg', 'DEC2000')
    SAGAv2 = SAGAv2[~np.isnan(SAGAv2['RA2000'])]
    
    SAGAv2['Star_ID'] = Column([f"{i:05.0f}" for i in range(1, len(SAGAv2)+1)])
    SAGAv2['From'] = np.full(len(SAGAv2), 'SAGA')

    coords = SkyCoord(ra=SAGAv2['RA2000'] * u.deg, dec=SAGAv2['DEC2000'] * u.deg)
    idx, d2d, d3d = coords.match_to_catalog_sky(coords)
    matches = (d2d < 1.0 * u.arcsec) & (idx != range(len(SAGAv2)))
    matched_sources = SAGAv2[matches]
    matched_indices = np.where(matches)[0]

    match_tuples = []
    for i, match in enumerate(matched_sources):
        matched_star_idx = matched_indices[i]
        og_star_idx = idx[matches][i]  # Get index of the matched star
        match_tuples.append((og_star_idx, matched_star_idx))

    for pair in match_tuples:
        SAGAv2['Star_ID'][pair[1]] = SAGAv2['Star_ID'][pair[0]]
                                                                                                       
    c = SkyCoord(ra=SAGAv2['RA2000']*u.degree, dec=SAGAv2['DEC2000']*u.degree)
    catalog = SkyCoord(ra=lit_cat['RAJ2000']*u.degree, dec=lit_cat['DEJ2000']*u.degree)

    match, d2d, d3d = c.match_to_catalog_sky(catalog)

    lit_inSAGA = lit_cat[match][d2d < 1*u.arcsec]
    SAGA_inlit = SAGAv2[d2d < 1*u.arcsec]
    SAGA_notinlit = SAGAv2[d2d >= 1*u.arcsec]

    key = Table()
    key['New'] = lit_inSAGA['Star_ID']
    key['Old'] = SAGA_inlit['Star_ID']
    key['Galaxy'] = lit_inSAGA['Galaxy']
    
    appends = np.unique(SAGA_notinlit['Star_ID'])

    corresponding = [f"{'Star'}_{i:05.0f}" for i in range(curr_count + 1, curr_count + len(appends)+ 1)]
    curr_count = curr_count + len(appends)
    key2 = Table()
    key2['New'] = corresponding
    key2['Old'] = appends
    key2['Galaxy'] = np.full(len(corresponding), 'nan')

    key = vstack([key, key2])
    mapping = {row['Old']: row['New'] for row in key}
    mapping2 = {row['Old']: row['Galaxy'] for row in key}
    new_names = []
    new_gals = []
    for old_name in SAGAv2['Star_ID']:
        new_name = mapping.get(old_name, None) 
        new_names.append(new_name)
        new_gal = mapping2.get(old_name, None)
        new_gals.append(new_gal)
    SAGAv2['Star_ID'] = new_names
    SAGAv2['Galaxy'] = new_gals
    SAGAv2 = SAGAv2[np.argsort(SAGAv2['Star_ID'])]
    

    #THIS MIGHT NEED TO BE REMOVED
    refers = []
    for entry in SAGAv2['Reference'].data:
        split1 = entry.split('.')
        short1 = split1[-1]
        split2 = short1.split(',')
        year = split2[-1]
        name = split2[0].strip('+')
        refers.append(name + year)
    SAGAv2['Paper'] = refers
    
    # Adding bibs
    bibs = []
    for i in range(len(SAGAv2['Reference'])):
        entry = SAGAv2['Reference'][i]
        split1 = entry.split(',')
        author= SAGAv2['Paper'][i].split(' ')[-2]
        year = split1[4].strip()
        journal = split1[1].strip()
        if journal == 'A+A':
            journal = 'A&A'
        volume = split1[2].strip()
        page = split1[3].strip()
        if not volume[-1].isnumeric():
            p = volume[-1]
            volume = volume[:-1]
        else:
            p = '.'
        if not (page[0]).isnumeric():
            p = page[0]
            page = page[1:]
        else:
            p = '.'
        bibs.append(f"{year}{journal:{'.'}<5}{volume:{'.'}>4}{p}{page:{'.'}>4}{author[0]}")
    SAGAv2['Bib'] = bibs


    # 2. Rename columns to match literature catalogue
    SAGAv2.rename_column('[Fe/H]', 'Fe_H')
    SAGAv2.rename_column('[Fe I/H]', 'FeI_H')
    SAGAv2.rename_column('[Fe II/H]', 'FeII_H')


    config = json.load(open("config.json", encoding="utf-8"))
    ALT_COL_NAMES = config["ALT_COL_NAMES"]
    META_DATA_KEY = config["META_DATA_KEY"]

    # Refer to config file for the column name conversion
    catalog.colnames = SAGAv2.columns
    col_list = {}
    new_column_found = False
    for col in catalog.colnames:
        try:
            col_list[col] = ALT_COL_NAMES[col]
        except KeyError:
            # If the column is not in the config file, ask the user to specify a name
            new_column_found = True
            msg_ = f"(enter _ to not include column '{col}' in the catalog)"
            msg = f"Found column '{col}' not in config file. Please specify a name {msg_}: "
            user_input = input(msg)
            if user_input == "":
                print("WARNING: Input empty. Operation aborted, no changes made.")
            col_list[col] = user_input
            config["ALT_COL_NAMES"][col] = col_list[col]
    curr_names = list(col_list.keys())
    new_names = list(col_list.values())
    for i in range(len(curr_names)):
        if new_names[i] == "_":
            SAGAv2.remove_column(curr_names[i])
        else:
            SAGAv2.rename_column(curr_names[i], new_names[i])


    # 4. Add The /Fe columns
    columns = []
    SAGAv2.columns
    for col in SAGAv2.columns:
        if 'H' in col:
            columns.append(col)

    iron_cols = ['Fe_H', 'FeI_H', 'FeII_H']
    columns.remove('Fe_H')
    columns.remove('FeI_H')
    columns.remove('FeII_H')
    columns.remove('[OI]_H')
    columns.remove('[SI]_H')
    columns.remove('M_H')

    given = list(saga_corrs.keys())
    subs = list(saga_corrs.values())   
    for col in columns:
        corr = subs[given.index(col)][0]
        backup1 = subs[given.index(col)][1]
        backup2 = subs[given.index(col)][2]
        element = col.split('_')[0]
        new_name = element + '_' + corr.split('_')[0]
        new_data = []
        for idx in range(len(SAGAv2)):
            if isinstance(SAGAv2[corr][idx], np.float64):
                new_data.append(SAGAv2[col][idx] - SAGAv2[corr][idx])
            else:
                if isinstance(SAGAv2[backup1][idx], np.float64):
                    new_data.append(SAGAv2[col][idx] - SAGAv2[backup1][idx])
                else:
                    if isinstance(SAGAv2[backup2][idx], np.float64):
                        new_data.append(SAGAv2[col][idx] - SAGAv2[backup2][idx])
                    else:
                        new_data.append(np.nan)
      
        SAGAv2.add_column(Column(new_data, name=new_name))
    
    # 5. Combine SAGA with Literature

    all_cols = set(lit_cat.columns) | set(SAGAv2.columns)

    # Add missing columns to table1
    for col in all_cols:
        if col not in lit_cat.columns:
            lit_cat[col] = ['nan'] * len(lit_cat)

    # Add missing columns to table2
    for col in all_cols:
        if col not in SAGAv2.columns:
            SAGAv2[col] = ['nan'] * len(SAGAv2)

    # Ensure columns are in the same order
    lit_cat = lit_cat[sorted(all_cols)]
    SAGAv2 = SAGAv2[sorted(all_cols)]
    lit_cols = list(lit_cat.columns)
    for col in lit_cols:
        if isinstance(SAGAv2[col],MaskedColumn):          
            if col == 'S_N' or col == 'binary' or col == 'solar_scale':
                new_col = SAGAv2[col].filled('nan')
            elif col == 'res':
                new_col = SAGAv2[col].filled(-9999)
            else:
                new_col = SAGAv2[col].filled(np.nan)

            SAGAv2[col] = Column(new_col)
        lit_cat[col] = lit_cat[col].astype(str)
        SAGAv2[col] = SAGAv2[col].astype(str)
    big_lit = vstack([lit_cat, SAGAv2])
    for col in big_lit.colnames:
        big_lit[col] = [value if value != '' else 'nan' for value in big_lit[col]]
        try:
            big_lit[col] = big_lit[col].astype(float)

        except:
            pass

    big_lit = big_lit[np.argsort(big_lit['Star_ID'])]
    data = big_lit['S_N'].data
    new = []
    for entry in data:
        
        # If isinstance(entry, string):
            #if it has > or < in it just remove the > or <
            # If it has - in it, then take the lower value
        # Also update this in the overleaf doc
        
        if not entry.isnumeric():
            if '>' in entry or '<' in entry:
                entry = entry[1:]
                entry = float(entry)
            elif '-' in entry:
                entry = entry.split('-')[1].strip()
            else:
                entry = np.nan
            new.append(np.nan)
        else:
            new.append(entry)
    big_lit['S_N'] = new
    big_lit['S_N'] = big_lit['S_N'].astype(float)
    big_lit.remove_column('ref')
    fix_types(big_lit)
    return(big_lit, curr_count)

def ji(lit_cat:Table, curr_count: int): 
    alex_df = pd.read_csv('dwarf_lit_all_coo.org', delimiter = '|', header = 1)
    alex = Table.from_pandas(alex_df)
    alex['From'] = np.full(len(alex), 'Alex Ji')
    alex.remove_column("Unnamed: 0")
    alex.remove_column('Unnamed: 32')
    alex.remove_column('        RAstr ')
    alex.remove_column('        Decstr ')
    alex.remove_column(' Full ')
    old_names = [' Star         ',' Galaxy    ','                 RA ','                 Dec ','      C ','      N ','     O ','     Na ','    Mg ','     Al ','     Si ','      K ','   CaI ','   ScII ','  TiII ','   CrI ','     Mn ','    Fe ','     Co ','     Ni ','     Zn ','     Sr ','     Ba ',' Eu     ',' Teff ',' logg ','   vt ',' Source                ']
    new_names = ['Name', 'Galaxy','RAJ2000','DEJ2000','C_Fe', 'N_Fe', 'O_Fe', 'Na_Fe', 'Mg_Fe', 'Al_Fe', 'Si_Fe', 'K_Fe', 'CaI_FeI', 'ScII_FeII', 'TiII_FeII', 'CrI_FeI', 'Mn_Fe', 'Fe_H', 'Co_Fe', 'Ni_Fe', 'Zn_Fe', 'Sr_Fe', 'Ba_Fe', 'Eu_Fe', 'T_eff', 'log_g', 'V_mic','Source']
    
    refs = ['Frebel14', 'Simon10', 'Francois16', 'Frebel10', 'Ji16a', 'FEL09', 'NOR10', 'ISH14', 'GIL13(GM)', 'LAI11', 'Norris10b', 'Koch08', 'Koch13', 'Roederer14', 'Ji16b', 'Roed16', 'Chiti18', 'Chiti22', 'Ji19', 'Venn17', 'Kirby17', 'Hansen17', 'Marshall19', 'Nagasawa17', 'Spite18', 'Ji20', 'Hansen20', 'Waller22', 'Webber23', 'Hansen24', 'FRE14']
    bibs = ['2014ApJ...786...74F', '2010ApJ...716..446S', '2016A&A...588A...7F', '2010IAUS..265..237F', '2016ApJ...817...41J', '2009A&A...508L...1F', '2010ApJ...711..350N', '2014A&A...562A.146I' ,'2013ApJ...763...61G', '2011ApJ...738...51L', '2010ApJ...711..350N', '2008ApJ...688L..13K', '2013A&A...554A...5K', '2014AJ....147..136R', '2016ApJ...830...93J', '2016AJ....151...82R', '2018ApJ...857...74C', '2022ApJ...939...41C', '2019ApJ...870...83J', '2017MNRAS.466.3741V', '2017ApJ...838...83K', '2017ApJ...838...44H', '2019ApJ...882..177M', '2018ApJ...852...99N', '2018A&A...617A..56S', '2020ApJ...889...27J', '2020ApJ...897..183H', '2023MNRAS.519.1349W', '2023ApJ...959..141W', '2024ApJ...968...21H','2014ApJ...786...74F']
    
    #new_names = []
    #banned = ['Star', 'Galaxy', 'Source', 'RA', 'Dec', 'Teff', 'logg', 'vt', 'From']
    #for col in alex.colnames:
    #    if col.strip() == 'Fe':
    #        new_names.append(col.strip() + '_H')
    #    elif col.strip() in banned:
    #        new_names.append(col.strip())
    #    else:
    #        new_names.append(col.strip() + '_Fe')

    abuns = ['C_Fe', 'N_Fe', 'O_Fe', 'Na_Fe', 'Mg_Fe', 'Al_Fe', 'Si_Fe', 'K_Fe', 'CaI_FeI', 'ScII_FeII', 'TiII_FeII', 'CrI_FeI', 'Mn_Fe',
             'Fe_H', 'Co_Fe', 'Ni_Fe', 'Zn_Fe', 'Sr_Fe', 'Ba_Fe', 'Eu_Fe', 'T_eff', 'log_g', 'V_mic']

    for idx in range(len(new_names)):
        alex.rename_column(old_names[idx], new_names[idx])
    #alex.rename_column('Star', 'Name')
    #alex.rename_column('RA', 'RAJ2000')
    #alex.rename_column('Dec', 'DEJ2000')
    for abun in abuns:
        alex.add_column(Column(np.full(len(alex), np.nan), name='er_' + abun), index=alex.colnames.index(abun) + 1)

    alex_bibs = []
    for entry in alex['Source']:
        if ',' in entry:
            temp = ''
            split1 = entry.split(',')
            for sub in split1:
                temp = temp + bibs[refs.index(sub.strip())] + ', '
            alex_bibs.append(temp[:(len(temp)-2)])
        else:
            alex_bibs.append(bibs[refs.index(entry.strip())])
    alex['Bib'] = alex_bibs
    # If value is an upper bound, the error is 9999
    # if value is a lower bound, the error is -9999
    #alex.rename_column('Teff', 'T_eff')
    #alex.rename_column('logg', 'log_g')
    #alex.rename_column('vt', 'V_mic')
    #alex.rename_column('er_Teff', 'er_T_eff')
    #alex.rename_column('er_logg', 'er_log_g')
    #alex.rename_column('er_vt', 'er_V_mic')

    for abun in abuns:
        new = []
        err = []
        data = alex[abun].data
        for entry in data:
            if np.ma.is_masked(entry):
                entry = ''
            entry = str(entry)
            if '<' in entry:
                err.append(-9999)
                flt = entry.replace('<', '').strip()
            elif '>' in entry:
                err.append(9999)
                flt = entry.replace('>', '').strip()
            else:
                err.append(np.nan)
                flt = entry.strip()
            if flt == '':
                new.append(np.nan)
            else:
                new.append(float(flt))
        alex[abun] = new
        alex['er_' + abun] = err


    # Change the format of the source column to fit the rest of the catalogue
    refers = []
    for ref in alex['Source'].data:
        new = ''
        split1 = ref.split(',')[0]
        chars = list(split1)
        new = new + (chars[1]).upper()
        for idx in range(len(chars[1:])):
            if chars[2+idx].isalpha():
                new = new + chars[2+idx].lower()
            else:
                new = new + ' 20' + chars[2+idx] + chars[3+idx]
                break
        refers.append(new)
    alex['Paper'] = refers
    alex.remove_column('Source')

    from collections import Counter

    alex['Star_ID'] = Column([f"{i:05.0f}" for i in range(1, len(alex)+1)])

    coords = SkyCoord(ra=alex['RAJ2000'] * u.deg, dec=alex['DEJ2000'] * u.deg)
    idx, d2d, d3d = coords.match_to_catalog_sky(coords)
    matches = (d2d < 1.0 * u.arcsec) & (idx != range(len(alex)))
    matched_sources = alex[matches]
    matched_indices = np.where(matches)[0]

    match_tuples = []
    for i, match in enumerate(matched_sources):
        matched_star_idx = matched_indices[i]
        og_star_idx = idx[matches][i]  # Get index of the matched star
        match_tuples.append((og_star_idx, matched_star_idx))

    for pair in match_tuples:
        alex['Star_ID'][pair[1]] = alex['Star_ID'][pair[0]]
   
    c = SkyCoord(ra=alex['RAJ2000']*u.degree, dec=alex['DEJ2000']*u.degree)
    catalog = SkyCoord(ra=lit_cat['RAJ2000']*u.degree, dec=lit_cat['DEJ2000']*u.degree)

    match, d2d, d3d = c.match_to_catalog_sky(catalog)

    lit_inalex = lit_cat[match][d2d < 1*u.arcsec]
    alex_inlit = alex[d2d < 1*u.arcsec]
    alex_notinlit = alex[d2d >= 1*u.arcsec]

    key = Table()
    key['New'] = lit_inalex['Star_ID']
    key['Old'] = alex_inlit['Star_ID']
    key['New_gals'] = lit_inalex['Galaxy']

    appends = np.unique(alex_notinlit['Star_ID'])

    corresponding = [f"{'Star'}_{i:05.0f}" for i in range(curr_count + 1, curr_count + len(appends)+ 1)]
    curr_count = curr_count + len(appends)
    
    key2 = Table()
    key2['New'] = corresponding
    key2['Old'] = appends
    key2['New_gals'] = alex_notinlit['Galaxy']

    key = vstack([key, key2])
    mapping = {row['Old']: row['New'] for row in key}
    mapping2 = {row['Old']: row['New_gals'] for row in key}
    new_names = []
    new_gals = []
    for old_name in alex['Star_ID']:
        new_name = mapping.get(old_name, None) 
        new_names.append(new_name)
        new_gal = mapping2.get(old_name, None)
        new_gals.append(new_gal)
    alex['Star_ID'] = new_names
    alex['Galaxy'] = new_gals
    alex = alex[np.argsort(alex['Star_ID'])]
    
    alex['Galaxy'] = alex['Galaxy'].astype('S15')
    fix_types(alex)
    fix_types(lit_cat)
    curr_gals = np.unique(alex['Galaxy']).tolist()
    for gal in curr_gals:
        alex['Galaxy'][alex['Galaxy']== gal]= renames[gal]
    key5 = Table()
    key5['old'] = list(objects.keys())
    key5['new'] = list(objects.values())

    mapping5 = {row['old']: row['new'] for row in key5}
    new_names = []

    for old_name in alex['Galaxy']:
        new_names.append(mapping5.get(old_name, None))
    alex['Object Type'] = np.array(new_names, dtype=str)
    #########
    
    all_cols = set(lit_cat.columns) | set(alex.columns)

    # Add missing columns to table1
    for col in all_cols:
        if col not in lit_cat.columns:
            lit_cat[col] = [np.nan] * len(lit_cat)

    # Add missing columns to table2
    for col in all_cols:
        if col not in alex.columns:
            alex[col] = [np.nan] * len(alex)

    fix_types(alex)
    fix_types(lit_cat)
    # Ensure columns are in the same order
    lit_cat = lit_cat[sorted(all_cols)]
    alex = alex[sorted(all_cols)]
    bigger_lit = vstack([lit_cat, alex])
    
    return(bigger_lit, curr_count)

def jina(lit_cat:Table, curr_count: int):
    JINA = ascii.read('/raid/users/heigerm/catalogues/jina.csv')
    JINA['From'] = np.full(len(JINA), 'JINA')
    # Remove excess columns
    remove = ['﻿JINA_ID','Simbad_Identifier','RA','DEC','Vel_bibcode','U_mag','B_mag','V_mag','R_mag','I_mag','J_mag',
              'H_mag','K_mag']
    for col in remove:
        JINA.remove_column(col)

    #temp1 = JINA['ra']
    #temp2 = JINA['dec']
    #JINA.rename_column('ra', 'RAJ2000')
    #JINA.rename_column('dec', 'DEJ2000')    
    # Rename remaining columns
    config = json.load(open("config.json", encoding="utf-8"))
    ALT_COL_NAMES = config["ALT_COL_NAMES"]
    META_DATA_KEY = config["META_DATA_KEY"]

    # Refer to config file for the column name conversion
    colnames = JINA.columns
    col_list = {}
    new_column_found = False
    for col in colnames:
        try:
            col_list[col] = ALT_COL_NAMES[col]
        except KeyError:
            # If the column is not in the config file, ask the user to specify a name
            new_column_found = True
            msg_ = f"(enter _ to not include column '{col}' in the catalog)"
            msg = f"Found column '{col}' not in config file. Please specify a name {msg_}: "
            user_input = input(msg)
            if user_input == "":
                print("WARNING: Input empty. Operation aborted, no changes made.")
            col_list[col] = user_input
            config["ALT_COL_NAMES"][col] = col_list[col]
    curr_names = list(col_list.keys())
    new_names = list(col_list.values())
    for i in range(len(curr_names)):
        if new_names[i] == "_":
            JINA.remove_column(curr_names[i])
        else:
            JINA.rename_column(curr_names[i], new_names[i])

    floats=['er_Fe_H','er_Ba_Fe','er_Li_Fe','er_Be_Fe','er_C_Fe','er_N_Fe','er_O_Fe','er_F_Fe','er_Na_Fe','er_Mg_Fe','er_Al_Fe',
            'er_Si_Fe','er_P_Fe','er_S_Fe','er_K_Fe','er_Ca_Fe','er_Sc_Fe','er_Ti_Fe','er_V_Fe','er_Cr_Fe','er_Mn_Fe','er_Co_Fe',
            'er_Ni_Fe','er_Cu_Fe','er_Zn_Fe','er_Ga_Fe','er_Ge_Fe','er_As_Fe','er_Se_Fe','er_Rb_Fe','er_Sr_Fe','er_Y_Fe','er_Zr_Fe',
            'er_Nb_Fe','er_Mo_Fe','er_Ru_Fe','er_Rh_Fe','er_Pd_Fe','er_Ag_Fe','er_Cd_Fe','er_Sn_Fe','er_Te_Fe','er_La_Fe',
            'er_Ce_Fe','er_Pr_Fe','er_Nd_Fe','er_Sm_Fe','er_Eu_Fe','er_Gd_Fe','er_Tb_Fe','er_Dy_Fe','er_Ho_Fe','er_Er_Fe',
            'er_Tm_Fe','er_Yb_Fe','er_Lu_Fe','er_Hf_Fe','er_W_Fe','er_Os_Fe','er_Ir_Fe','er_Pt_Fe','er_Au_Fe','er_Hg_Fe','er_Pb_Fe',
            'er_Bi_Fe','er_Th_Fe','er_U_Fe','er_CaII_FeII','er_TiII_FeII','er_VII_FeII','er_CrII_FeII','er_MnII_FeII', 'T_eff', 'F_Fe',
            'As_Fe','Se_Fe', 'W_Fe', 'Hg_Fe', 'Bi_Fe']

    for col in floats:
        JINA[col] = JINA[col].astype(float)

    for col in JINA.columns:
        if isinstance(JINA[col],MaskedColumn):          
            if col == 'Paper' or col == 'Name':
                new_col = JINA[col].filled('nan')
            else:
                new_col = JINA[col].filled(np.nan)

            JINA[col] = Column(new_col)

    lims=['er_Fe_H','er_Ba_Fe','er_Li_Fe','er_Be_Fe','er_C_Fe','er_N_Fe','er_O_Fe','er_F_Fe','er_Na_Fe','er_Mg_Fe','er_Al_Fe',
            'er_Si_Fe','er_P_Fe','er_S_Fe','er_K_Fe','er_Ca_Fe','er_Sc_Fe','er_Ti_Fe','er_V_Fe','er_Cr_Fe','er_Mn_Fe','er_Co_Fe',
            'er_Ni_Fe','er_Cu_Fe','er_Zn_Fe','er_Ga_Fe','er_Ge_Fe','er_As_Fe','er_Se_Fe','er_Rb_Fe','er_Sr_Fe','er_Y_Fe','er_Zr_Fe',
            'er_Nb_Fe','er_Mo_Fe','er_Ru_Fe','er_Rh_Fe','er_Pd_Fe','er_Ag_Fe','er_Cd_Fe','er_Sn_Fe','er_Te_Fe','er_La_Fe',
            'er_Ce_Fe','er_Pr_Fe','er_Nd_Fe','er_Sm_Fe','er_Eu_Fe','er_Gd_Fe','er_Tb_Fe','er_Dy_Fe','er_Ho_Fe','er_Er_Fe',
            'er_Tm_Fe','er_Yb_Fe','er_Lu_Fe','er_Hf_Fe','er_W_Fe','er_Os_Fe','er_Ir_Fe','er_Pt_Fe','er_Au_Fe','er_Hg_Fe','er_Pb_Fe',
            'er_Bi_Fe','er_Th_Fe','er_U_Fe','er_CaII_FeII','er_TiII_FeII','er_VII_FeII','er_CrII_FeII','er_MnII_FeII']

    for col in lims:
        data = JINA[col].data
        for idx in range(len(data)):
            if data[idx] == 1:
                data[idx] = 9999
            elif data[idx] == 0:
                data[idx] = -9999
        JINA[col] = data

    JINA['Galaxy'] = np.full(len(JINA), 'nan')
    
    JINA['Star_ID'] = Column([f"{i:05.0f}" for i in range(1, len(JINA)+1)])

    coords = SkyCoord(ra=JINA['RAJ2000'] * u.deg, dec=JINA['DEJ2000'] * u.deg)
    idx, d2d, d3d = coords.match_to_catalog_sky(coords)
    matches = (d2d < 1.0 * u.arcsec) & (idx != range(len(JINA)))
    matched_sources = JINA[matches]
    matched_indices = np.where(matches)[0]

    match_tuples = []
    for i, match in enumerate(matched_sources):
        matched_star_idx = matched_indices[i]
        og_star_idx = idx[matches][i]  # Get index of the matched star
        match_tuples.append((og_star_idx, matched_star_idx))

    for pair in match_tuples:
        JINA['Star_ID'][pair[1]] = JINA['Star_ID'][pair[0]]

    c = SkyCoord(ra=JINA['RAJ2000']*u.degree, dec=JINA['DEJ2000']*u.degree)
    catalog = SkyCoord(ra=lit_cat['RAJ2000']*u.degree, dec=lit_cat['DEJ2000']*u.degree)

    match, d2d, d3d = c.match_to_catalog_sky(catalog)

    lit_inJINA = lit_cat[match][d2d < 1*u.arcsec]
    JINA_inlit = JINA[d2d < 1*u.arcsec]
    JINA_notinlit = JINA[d2d >= 1*u.arcsec]

    key = Table()
    key['New'] = lit_inJINA['Star_ID']
    key['Old'] = JINA_inlit['Star_ID']
    key['Galaxy'] = lit_inJINA['Galaxy']

    appends = np.unique(JINA_notinlit['Star_ID'])

    corresponding = [f"{'Star'}_{i:05.0f}" for i in range(curr_count + 1, curr_count + len(appends)+ 1)]
    curr_count = curr_count + len(appends)
    
    key2 = Table()
    key2['New'] = corresponding
    key2['Old'] = appends
    key2['Galaxy'] = np.full(len(corresponding), 'nan')

    key = vstack([key, key2])
    mapping = {row['Old']: row['New'] for row in key}
    mapping2 = {row['Old']: row['Galaxy'] for row in key}
    new_names = []
    new_gals = []
    for old_name in JINA['Star_ID']:
        new_name = mapping.get(old_name, None) 
        new_names.append(new_name)
        new_gal = mapping2.get(old_name, None)
        new_gals.append(new_gal)
    JINA['Star_ID'] = new_names
    JINA['Galaxy'] = new_gals
    JINA = JINA[np.argsort(JINA['Star_ID'])]

    # Fix paper naming conventions
    key3 = Table()
    key3['New'] = [value[0] for value in paper_renames.values()]
    key3['bibs'] = [value[1] for value in paper_renames.values()]
    key3['Old'] = list(paper_renames.keys())
    mapping3 = {row['Old']: row['New'] for row in key3}
    mapping4 = {row['Old']: row['bibs'] for row in key3}
    new_names = []
    new_bibs = []
    for old_name in JINA['Paper']:
        new_name = mapping3.get(old_name, None)
        new_bib = mapping4.get(old_name, None)
        new_names.append(new_name)
        new_bibs.append(new_bib)
    JINA['Paper'] = new_names
    JINA['Bib'] = new_bibs

    all_cols = set(lit_cat.columns) | set(JINA.columns)

    # Add missing columns to table1
    for col in all_cols:
        if col not in lit_cat.columns:
            lit_cat[col] = [np.nan] * len(lit_cat)

    # Add missing columns to table2
    for col in all_cols:
        if col not in JINA.columns:
            JINA[col] = [np.nan] * len(JINA)

    # Ensure columns are in the same order
    lit_cat = lit_cat[sorted(all_cols)]
    JINA = JINA[sorted(all_cols)]
    fix_types(JINA)
    fix_types(lit_cat)
    biggest_lit = vstack([lit_cat, JINA])
    return(biggest_lit, curr_count)

def cust_stack(incoming: Table, lit_cat: Table, Name: str, curr_count:int):
    if 'Star_ID' not in incoming.columns:
        incoming['Star_ID'] =  Column([f"{i:05.0f}" for i in range(len(incoming))])
    if 'Galaxy' not in incoming.columns:
        incoming['Galaxy'] = ['nan' for i in range (len(incoming))]
    galaxies(incoming)    

    coords = SkyCoord(ra=incoming['RAJ2000'] * u.deg, dec=incoming['DEJ2000'] * u.deg)
    idx, d2d, d3d = coords.match_to_catalog_sky(coords)
    matches = (d2d < 1.0 * u.arcsec) & (idx != range(len(incoming)))
    matched_sources = incoming[matches]
    matched_indices = np.where(matches)[0]

    match_tuples = []
    for i, match in enumerate(matched_sources):
        matched_star_idx = matched_indices[i]
        og_star_idx = idx[matches][i]  # Get index of the matched star
        match_tuples.append((og_star_idx, matched_star_idx))
    
    for pair in match_tuples:
        incoming['Star_ID'][pair[1]] = incoming['Star_ID'][pair[0]]      
            
    c = SkyCoord(ra=incoming['RAJ2000']*u.degree, dec=incoming['DEJ2000']*u.degree)
    catalog = SkyCoord(ra=lit_cat['RAJ2000']*u.degree, dec=lit_cat['DEJ2000']*u.degree)

    match, d2d, d3d = c.match_to_catalog_sky(catalog)

    lit_inInc = lit_cat[match][d2d < 1*u.arcsec]
    Inc_inlit = Inc[d2d < 1*u.arcsec]
    Inc_notinlit = Inc[d2d >= 1*u.arcsec]

    key = Table()
    key['New'] = lit_inInc['Star_ID']
    key['Old'] = Inc_inlit['Star_ID']
    key['Galaxy'] = lit_inInc['Galaxy']
 
    key2 = Table()
    key2['New'] = [f"{'Star'}_{i:05.0f}" for i in range(curr_count + 1, curr_count + len(appends)+ 1)]
    key2['Old'] = Inc_notinlit['Star_ID']
    key2['Galaxy'] = Inc_notinlit['Galaxy']
    curr_count = curr_count + len(appends)

    key = vstack([key, key2])
    mapping = {row['Old']: row['New'] for row in key}
    mapping2 = {row['Old']: row['Galaxy'] for row in key}
    new_names = []
    new_gals = []
    for old_name in incoming['Star_ID']:
        new_name = mapping.get(old_name, None) 
        new_names.append(new_name)
        new_gal = mapping2.get(old_name, None)
        new_gals.append(new_gal)
    incoming['Star_ID'] = new_names
    incoming['Galaxy'] = new_gals
    incoming = incoming[np.argsort(incoming['Star_ID'])]
    
    all_cols = set(lit_cat.columns) | set(incoming.columns)
    # Add missing columns to table1
    for col in all_cols:
        if col not in lit_cat.columns:
            lit_cat[col] = ['nan'] * len(lit_cat)
        if col not in incoming.columns:
            incoming[col] = ['nan'] * len(incoming)

    # Ensure columns are in the same order
    lit_cat = lit_cat[sorted(all_cols)]
    incoming = incoming[sorted(all_cols)]
    lit_cols = list(lit_cat.columns)
    
    for col in lit_cols:
        if col != 'RAJ2000' and col != 'DEJ2000':
            lit_cat[col] = lit_cat[col].astype(str)
            incoming[col] = incoming[col].astype(str)
            
    print(Name, ' merged!\n')
    return (vstack([lit_cat, incoming]), curr_count)

def gaia(lit_cat:Table):
    tbl = lit_cat
    gaia_table = 'gaia_dr3.gaia_source'
    gaia_data = crossmatcher.doit(gaia_table, tbl['RAJ2000'], tbl['DEJ2000'], '*', radeccols=('ra', 'dec'), rad=1.,\
                                        extra=None,yourradeccols=('ra_gaia', 'dec_gaia'), asDict=True)
    tmp = table.Table(gaia_data)
    tmp['source_id'] = tmp['source_id'].astype(int)
    lit_cat.add_column(tmp['source_id'], index=2, name='Gaia_ID')
    lit_cat['Gaia_ID'] = lit_cat['Gaia_ID'].astype(int)
    
