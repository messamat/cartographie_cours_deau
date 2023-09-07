from setup_classement import * #Get directory structure and packages

overwrite = False #Whether to overwrite outputs or skip them if done

#Read compilation of metadata
geometadata_path = list(datdir.glob('metadonnes_cartographie_cours_deau*xlsx'))[-1] #Get most recent table
geometadata_pd = pd.read_excel(geometadata_path, sheet_name = 'Métadonnées_réseau_SIG')

#Create output directory
out_dir = Path(resdir, 'reseaux_departementaux_copies')

#Copy all layers for renaming
net_copylist = getfilelist(out_dir, repattern='.*_copie[.](TAB|shp)')

# in_net_path = net_editlist[0]
# in_geometadata_pd = geometadata_pd
def format_net(in_net_path, in_geometadata_pd):
    print(in_net_path)

    #Get row in metadata tab
    depnum = int(re.sub("[A-zÀ-ÿ_']", '', Path(in_net_path).parent.name))
    row_metadata = in_geometadata_pd.loc[in_geometadata_pd['Numéro']==depnum].iloc[0,:]

    #Export to shapefile
    out_file = re.sub(
        '[.](TAB|shp)$', '.gpkg',
        re.sub('_copie(?=[.][a-zA-Z]{2,3})', '_edit', in_net_path)
    )

    if (not Path(out_file).exists()) or overwrite:
        print(out_file)
        net = gpd.read_file(in_net_path, encodings=row_metadata['Encodage'])
        #gpkg requires fid column to be int - sometimes read in as float from shp
        if 'fid' in net.columns:
            if not pd.api.types.is_integer_dtype(net.fid):
                net['fid'] = net['fid'].astype(int)
            if any(net.fid.duplicated()):
                net.rename(columns={"fid":"fid_orig"}, inplace=True)
        try:
            net.to_file(out_file, encoding='utf-8')
        #Some layers have mixed encoding
        except ValueError:
            print('Probable mixed encoding. Removing problematic records...')
            for col in np.where(net.dtypes=='object')[0]:
                for i in range(len(net)):
                    record = net.iloc[i, col]
                    if not pd.isnull(record):
                        if isinstance(record, bytes):
                            #print(col)
                            net.iloc[i, col] = np.nan
            net.to_file(out_file, encoding='utf-8')

        # def detect_encoding_series(in_series, in_origencoding):
        #     for record in in_series:
        #         if not pd.isnull(record):
        #             print(record)
        #             if not isinstance(record, bytes):
        #                 record = record.encode(in_origencoding)
        #             check = chardet.detect(record)
        #             if check['confidence'] > 0.99:
        #                 break
        #     return(check['encoding'])
        # net.select_dtypes(include='object').apply(
        #     lambda col: col.\
        #     str.encode(row_metadata['Encodage']).\
        #     str.decode(detect_encoding_series(col, row_metadata['Encodage'])))


    #Read in editable layer
    net = gpd.read_file(out_file)
    net.columns = [col.lower() for col in net.columns]

    #Rename attributes
    renaming_dict = {
        "Nom de l'attribut désignant le type d'écoulement": "type_stand",
        "Nom de l'attribut auxiliaire désignant le type d'écoulement": "type_aux",
        "Nom de l'attribut désignant le régime hydrologique": "regime",
        "Nom de l'attribut désignant la méthode d'identification de l'écoulement": "nat_id",
        "Nom de l'attribut désignant la source de la modification, de la suppression du tronçon BD TOPO®, ou de l’ajout d’un nouveau tronçon": "orig_mo",
        "Nom de l'attribut désignant la date de l'identification du type d'écoulement": "date_id"
    }

    for k,v in renaming_dict.items():
        if not pd.isnull(row_metadata[k].iloc[0]):
            colnames_orig_list = split_strip(row_metadata[k].iloc[0])
            for i in range(len(colnames_orig_list)):
                if i == 0:
                    net[v] = net[colnames_orig_list[i].lower()]
                else:
                    net['{0}{1}'.format(v, i+1)] = net[colnames_orig_list[i].lower()]

    #Recategorize TYPE_ECOUL
    #row_metadata.columns
    dict_recat_type_ecoul = {'CE_recatégorisé': "Cours d'eau",
                             'NCE_recatégorisé': "Non cours d'eau",
                             'Indéterminé_ou_autre_recatégorisé': "Indéterminé",
                             'Inexistant_recatégorisé': "Inexistant",
                             'Hors_département_recatégorisé': "Hors département"}


    # in_row_net = net.iloc[0,:]
    # in_row_metadata = row_metadata
    # in_dict_recat = dict_recat_type_ecoul
    def recat_gpdcol(in_row_net, in_row_metadata, in_dict_recat):
        for cats_orig_colname, cat_new in in_dict_recat.items():
            #print(in_row_net['TYPE_ECOUL'])
            if not pd.isnull(in_row_metadata[cats_orig_colname].iloc[0]):
                cats_orig = [
                    cat if cat != 'NULL' else np.nan
                    for cat in split_strip(
                        str(in_row_metadata[cats_orig_colname].iloc[0]),
                        sep=';')
                ]
                #print(cats_orig)
                if in_row_net['type_stand'] in cats_orig:
                    #print(cat_new)
                    return (cat_new)
                    break


    net['type_stand'] = net.apply(recat_gpdcol,
                                  in_row_metadata=row_metadata,
                                  in_dict_recat=dict_recat_type_ecoul,
                                  axis = 1)

    # Reproject all to the same layer 2154
    if net.crs.to_epsg() != 2514:
        net_proj = net.to_crs("epsg:2154")
    else:
        net_proj = net

    #Write out results
    net_proj.to_file(out_file)

for net_path in net_copylist:
    format_net(in_net_path=net_path,
               in_geometadata_pd=geometadata_pd
               )

#Merge into one layer


#Intersect with administrative boundaries
#Intersect with hydrographic units boundaries