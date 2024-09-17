from setup_classement import * #Get directory structure and packages

overwrite = True #Whether to overwrite outputs or skip them if done
datdir = Path(datdir)

#Read compilation of metadata
geometadata_path = list(datdir.glob('metadonnes_cartographie_cours_deau*xlsx'))[-1] #Get most recent table
geometadata_pd = pd.read_excel(geometadata_path, sheet_name = 'Métadonnées_réseau_SIG',
                               na_values = 'NA', keep_default_na = False)
#Remove duplicate layers (for Ile-de-France)
geometadata_pd = geometadata_pd.loc[~(geometadata_pd['Lien local données'].duplicated())]

#Create output directory
out_dir = Path(resdir, 'reseaux_departementaux_copies')
if not out_dir.exists():
    out_dir.mkdir()

#Output
out_QCtab = Path(resdir, 'QC_nets.csv')

yonne_out_mergedfile = Path(rootdir,
                            geometadata_pd[geometadata_pd['Numéro'] == 89]['Lien local données'].iloc[0])
yonne_in_dir = Path(rootdir,
                    Path(geometadata_pd[geometadata_pd['Numéro'] == 89]['Lien local données'].iloc[0]).parent)

#in_row = geometadata_pd.iloc[92,:]
def copy_net(in_row, in_rootdir, out_dir, quiet, overwrite):
    if not quiet:
        print(in_row['Numéro'])
    depdir = Path(out_dir, "D{0}_{1}".format(
        in_row.Numéro,
        re.sub('[-]', '_', in_row.Département)
    ))

    if not depdir.exists():
        os.mkdir(depdir)

    if (in_row["Lien local données"] != 'TBD') and (in_row['Titre du fichier'] != 'TBD'):
        if in_rootdir is not None:
            orig_dir = Path(in_rootdir, Path(in_row["Lien local données"]).parent)
        else:
            orig_dir = Path(in_row["Lien local données"]).parent

        for in_file in getfilelist(dir=orig_dir,
                                   repattern='^{0}[.][a-zA-Z]{{2,3}}$'.format(Path(in_row["Lien local données"]).stem)):
            out_file = Path(depdir,
                            '{0}_copie{1}'.format(Path(in_file).stem,
                                                  Path(in_file).suffix)
                            )
            if not out_file.exists() or overwrite:
                if not quiet:
                    print(out_file)
                shutil.copy(in_file, out_file)
            else:
                if not quiet:
                    print('{} already exists.'.format(out_file))

        return(
            Path(depdir,
                 '{0}{1}'.format(
                     out_file.stem,
                     Path(in_row["Lien local données"]).suffix)
                 )
        )
    else:
        return(np.nan)

#Define functions
def remove_invalid_geom(ds):
    if isinstance(ds, Path):
        ds = str(ds)

    if os.path.splitext(ds)[1] == '.shp':
        net_ogr = ogr.Open(ds, update=1)
        lyr = net_ogr.GetLayer()

    elif os.path.splitext(ds)[1] == '.TAB':
        driver = ogr.GetDriverByName("MapInfo File")
        net_ogr = driver.Open(ds, 1)
        lyr = net_ogr.GetLayer(0)

    for feature in lyr:
        geom = feature.GetGeometryRef()
        if geom is None:
            print('Deleting FID {}...'.format(feature.GetFID()))
            lyr.DeleteFeature(feature.GetFID())
        elif geom.IsValid() == False:
            print('Deleting FID {}...'.format(feature.GetFID()))
            lyr.DeleteFeature(feature.GetFID())

    net_ogr = None

def check_colmatch(in_colrecord, net):
    if pd.isnull(in_colrecord):
        return(True)
    else:
        return(all([(col in net.columns)
                    for col in split_strip(in_colrecord)]))

#in_row = geometadata_pd.iloc[55,:]

def QC_row_metadata(in_row, in_dict):
    print(in_row['Numéro'])

    if not pd.isnull(in_row['data_copy_path']):
        #Read in network feature class
        try:
            net = gpd.read_file(in_row["data_copy_path"], encodings=in_row['Encodage'])
        except:
            # Clean network of invalid geometries first
            remove_invalid_geom(in_row['data_copy_path'])
            # Then read it
            net = gpd.read_file(in_row["data_copy_path"], encodings=in_row['Encodage'])

        #Format categories of ecoulement from metadata table
        ce_cats_raw = in_row["Catégories de l'attribut désignant le type d'écoulement"]
        if not pd.isnull(ce_cats_raw):
            if ";" in ce_cats_raw:
                unique_typecoul_tab = split_strip(ce_cats_raw, ';')
            else:
                unique_typecoul_tab = split_strip(ce_cats_raw)

        type_ecoul_colname_tab = in_row["Nom de l'attribut désignant le type d'écoulement"]


        in_dict[in_row.Numéro] = [
            #Number of features match - n_match
            "{0}: {1}-{2}".format(
                in_row["Nombre total d'écoulements référencés"] == len(net),
                in_row["Nombre total d'écoulements référencés"],
                len(net)
            ),
            # all columns in network are recorded - col_match
            [col for col in net.columns if ((col != 'geometry') and (not col in in_row["Nom d'attributs"]))],
            # cours d'eau attribute name exists - ce_colname_match
            type_ecoul_colname_tab in net.columns or pd.isnull(type_ecoul_colname_tab),
            # cours d'eau categories match - ce_cats_match
            [cat for cat in net[type_ecoul_colname_tab].unique() if (not cat in unique_typecoul_tab)] \
                if not pd.isnull(in_row["Catégories de l'attribut désignant le type d'écoulement"]) \
                else True,

            check_colmatch(in_row["Nom de l'attribut auxiliaire désignant le type d'écoulement"], net), #auxi_colname_match
            check_colmatch(in_row["Nom de l'attribut désignant le régime hydrologique"], net), #'regime_colname_match'
            check_colmatch(in_row["Nom de l'attribut désignant la méthode d'identification de l'écoulement"], net), #natident_colname_match
            check_colmatch(in_row["Nom de l'attribut désignant la source de la modification; de la suppression du tronçon BD TOPO; ou de l’ajout d’un nouveau tronçon"], net), #origmodif_colname_match
            check_colmatch(in_row["Nom de l'attribut désignant la date de l'identification du type d'écoulement"], net) #dateident_colname_match
        ]
    else:
        return(np.nan)

#Format Yonne — merge layers
if (not Path(yonne_out_mergedfile).exists()) or overwrite:
    if Path(yonne_out_mergedfile).exists():
        for fname in getfilelist(yonne_in_dir, '{}*'.format(os.path.splitext(yonne_out_mergedfile.name)[0])):
            os.remove(fname)
    yonne_in_filepaths = getfilelist(yonne_in_dir, '.*[.](TAB|shp)$')
    pd.concat(
        [convert_bytes_to_na(
            gpd.read_file(file_path,
                          encoding=geometadata_pd[geometadata_pd['Numéro'] == 89]['Encodage'].iloc[0]). \
                assign(status_lyr = Path(file_path).stem)
        )
            for file_path in yonne_in_filepaths]
    ).pipe(gpd.GeoDataFrame). \
        to_file(Path(yonne_out_mergedfile))

#Copy data to results folder for manipulation
geometadata_pd['data_copy_path'] = geometadata_pd.apply(copy_net,
                                                        in_rootdir=rootdir,
                                                        out_dir=out_dir,
                                                        quiet=False,
                                                        overwrite=overwrite,
                                                        axis=1)

#QC metadata
if not out_QCtab.exists() or overwrite:
    QC_dict = defaultdict(list)
    geometadata_QCed_pd = geometadata_pd.apply(QC_row_metadata,
                                               in_dict = QC_dict,
                                               axis=1)

    QC_pd = pd.DataFrame.from_dict(QC_dict, orient='index')
    QC_pd.columns = ['n_match',
                     'col_match',
                     'ce_colname_match',
                     'ce_cats_match',
                     'auxi_colname_match',
                     'regime_colname_match',
                     'natident_colname_match',
                     'origmodif_colname_match',
                     'dateident_colname_match'
                     ]
    QC_pd.to_csv(out_QCtab)

#Check for typos, etc.
uvals = geometadata_pd.apply(lambda col: col.unique())
# in_row = geometadata_pd.loc[geometadata_pd.Numéro==25]
# net = gpd.read_file(in_row["data_copy_path"].values[0], encodings=in_row['Encodage'].values[0])
# net.columns
# geometadata_pd.loc[geometadata_pd.Numéro==25]["Nom de l'attribut auxiliaire désignant le type d'écoulement"].values[0].split(',').
#
# [col for col in split_strip(geometadata_pd.loc[geometadata_pd.Numéro==25]["Nom de l'attribut auxiliaire désignant le type d'écoulement"].values[0])
#  if col not in net.columns]
# #uvals['Commentaire']






