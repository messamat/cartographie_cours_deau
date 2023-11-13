import fiona
import shapely
from setup_classement import * #Get directory structure and packages

overwrite = True #Whether to overwrite outputs or skip them if done

datdir = Path(datdir)
resdir = Path(resdir)

#Read compilation of metadata
geometadata_path = list(datdir.glob('metadonnes_cartographie_cours_deau*xlsx'))[-1] #Get most recent table
geometadata_pd = pd.read_excel(geometadata_path, sheet_name = 'Métadonnées_réseau_SIG')

#Create output directory
out_dir = Path(resdir, 'reseaux_departementaux_copies')
out_net_path = Path(resdir, 'carto_loi_eau_france.gpkg')

#Copy all layers for renaming
net_copylist = getfilelist(out_dir, repattern='.*_copie[.](TAB|shp|gpkg)$')

# in_net_path = net_copylist[93]
# in_geometadata_pd = geometadata_pd
def create_editable_lyr(in_copylyr, out_editlyr, in_encoding, overwrite):
    if (not Path(out_editlyr).exists()) or overwrite:
        print(out_editlyr)
        try:
            net = gpd.read_file(in_copylyr,
                                encoding=in_encoding.lower(),
                                )
        #If one of the records has single-point LineStrings, remove these individual LineStrings and re-read
        except shapely.errors.GEOSException: #If one of the records has invalid geometry that prevents reading
            print("Removing Linestrings with a single point...")
            out_fixed_net = "{0}_fixed{1}".format(os.path.splitext(in_copylyr)[0],
                                                  os.path.splitext(in_copylyr)[1])
            with fiona.open(in_copylyr) as input:
                with fiona.open(out_fixed_net, 'w', 'ESRI Shapefile',
                                input.schema.copy(), input.crs) as output:
                    #x=0
                    for elem in input:
                        # GeoJSON to shapely geometry
                        if elem['geometry'] is not None:
                            if elem['geometry'].type == 'MultiLineString':
                                out_coords = [l for l in elem['geometry'].coordinates if len(l) > 1]
                                if len(out_coords) > 0:
                                    out_type = 'MultiLineString'
                                else:
                                    out_type = 'LineString'
                                out_elem=  {
                                    'geometry': {'type':out_type,
                                                 'coordinates': out_coords
                                                 },
                                    'properties': elem.properties
                                }
                                output.write(out_elem)

                            elif elem['geometry'].type == 'LineString':
                                if len(elem['geometry'].coordinates) > 1:
                                    output.write(elem)
                        #x += 1
            net = gpd.read_file(out_fixed_net,
                                encoding=in_encoding.lower(),
                                )

        #gpkg requires fid column to be int - sometimes read in as float from shp
        if 'fid' in net.columns:
            if any(net.fid.duplicated()):
                net.rename(columns={"fid":"fid_orig"}, inplace=True)
            elif not pd.api.types.is_integer_dtype(net.fid):
                net['fid'] = net['fid'].astype(int)
        try:
            net.to_file(out_editlyr, encoding='utf-8')
        #Some layers have mixed encoding
        except ValueError:
            print('Probable mixed encoding. Removing problematic records...')
            convert_bytes_to_na(net).to_file(out_editlyr, encoding='utf-8')

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

def regex_dir_and_file(in_path):
    return(r".*{0}.*{1}$".format(Path(in_path).parent.stem,
                                 re.sub('_copie[.]', '.', Path(in_path).name)
                                 )
           )

def format_net(in_net_path, in_geometadata_pd, overwrite):
    print(in_net_path)

    #Get row in metadata tab
    #depnum = int(re.sub("[A-zÀ-ÿ_']", '', Path(in_net_path).parent.name))
    if len(np.where(in_geometadata_pd['Lien local données'].str.match(regex_dir_and_file(in_net_path)))[0]) != 0:
        row_metadata = in_geometadata_pd.loc[
                       in_geometadata_pd['Lien local données'].str.match(regex_dir_and_file(in_net_path))].\
                           iloc[0,:]

        #Export to shapefile
        out_file = re.sub(
            '[.](TAB|shp)$', '.gpkg',
            re.sub('_copie(?=[.][a-zA-Z]{2,3})', '_edit', in_net_path)
        )

        create_editable_lyr(in_copylyr=in_net_path,
                            out_editlyr=out_file,
                            in_encoding=row_metadata['Encodage'],
                            overwrite=overwrite)

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
            if not pd.isnull(row_metadata[k]):
                colnames_orig_list = split_strip(row_metadata[k])
                for i in range(len(colnames_orig_list)):
                    if i == 0:
                        net[v] = net[colnames_orig_list[i].lower()]
                    else:
                        net['{0}{1}'.format(v, i+1)] = net[colnames_orig_list[i].lower()]

        #Delete duplicate lines, keeping the ones with the most information
        net.drop_duplicates(subset=[col for col in net.columns if col!='status_lyr'],
                            keep='first',
                            inplace=True)
        if len(net.duplicated(subset='geometry')) > 0:
            print('Removing duplicate geometries...')
            net['number_of_nonna_attricols'] = net.apply(lambda in_row:
                                                         sum([(not (in_row.loc[col]==np.nan))
                                                              for col in renaming_dict.values()
                                                              if col in net.columns]),
                                                         axis=1)
            net = net.sort_values('number_of_nonna_attricols', ascending=False).\
                drop_duplicates('geometry', keep='first').\
                sort_index().\
                drop('number_of_nonna_attricols', axis=1)

        #Recategorize TYPE_ECOUL
        #row_metadata.columns
        dict_recat_type_ecoul = {'CE_recatégorisé': "Cours d'eau",
                                 'NCE_recatégorisé': "Non cours d'eau",
                                 'Indéterminé_ou_autre_recatégorisé': "Indéterminé",
                                 'Inexistant_recatégorisé': "Inexistant",
                                 'Hors_département_recatégorisé': "Hors département"}


        def recat_gpdcol(in_row_net, in_row_metadata, in_dict_recat):
            if pd.isnull(in_row_metadata["Nom de l'attribut désignant le type d'écoulement"]):
                return("Cours d'eau")
            else:
                for cats_orig_colname, cat_new in in_dict_recat.items():
                    #print(in_row_net['type_ecoul'])
                    if in_row_metadata[cats_orig_colname] != "Pas de catégorie correspondante":
                        cats_orig = split_strip(
                            str(in_row_metadata[cats_orig_colname]),
                            sep=';')
                        if 'NULL' in cats_orig:
                            cats_orig.append(str(np.nan))

                        if not pd.isnull(in_row_net['type_stand']):
                            if isinstance(in_row_net['type_stand'], float):
                                in_row_net.loc['type_stand'] = int(in_row_net['type_stand'])

                        if str(in_row_net['type_stand']) in cats_orig:
                            #print(cat_new)
                            return(cat_new)
                            break


        net.loc[:, 'type_stand'] = net.apply(recat_gpdcol,
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

    else:
        print('No metadata associated with this layer. Skipping.')

for net_path in net_copylist[50:]:
    #print(net_path)
    format_net(in_net_path=net_path,
               in_geometadata_pd=geometadata_pd,
               overwrite=overwrite
               )

#Format Aisne -------------------------
print('Formatting data for the Aisne department')
aisne_editfile = getfilelist(Path(out_dir, 'D2_Aisne'), '.*_edit[.]gpkg')[0]
aisne_gpd = gpd.read_file(aisne_editfile)

#Merge fields indicating status
aisne_gpd.loc[aisne_gpd['type_ecoul'].isna(), 'type_ecoul'] = aisne_gpd.loc[aisne_gpd['type_ecoul'].isna(),'conclusion']
aisne_cat_dict = {
    "Cours d'eau": "Cours d'eau",
    'Indéterminé': 'Indéterminé',
    'CE': "Cours d'eau",
    'IndÃ©terminÃ©': 'Indéterminé',
    'NON': "Non cours d'eau",
    'SUPPRESSION':"Non cours d'eau",
    'MAINTIENT':"Cours d'eau",
    "AJOUT":"Cours d'eau",
    "MODIF":"Cours d'eau",
    "MAINTIEN":"Cours d'eau",
    "MODIFICATION":"Cours d'eau",
    "Permanent":"Cours d'eau"
}

aisne_gpd.loc[:,'type_stand'] = aisne_gpd.loc[:,'type_ecoul'].apply(
    lambda cat_orig: aisne_cat_dict[cat_orig] if not pd.isnull(cat_orig) else "Cours d'eau" #Cours d'eau by default
)

#Convert buffer polygons to centerlines
def gpd_centerline(in_row):
    if in_row.geometry.geom_type in ['Polygon', 'MultiPolygon']:
        polygon = in_row.geometry
        centerline = pygeoops.centerline(polygon, simplifytolerance=0) #https://pygeoops.readthedocs.io/en/latest/api/pygeoops.centerline.html#pygeoops.centerline
        in_row.loc['geometry'] = centerline
    return(in_row)
if aisne_gpd.geom_type.isin(['Polygon', 'MultiPolygon']).any():
    aisne_gpd = aisne_gpd.apply(gpd_centerline, axis=1)
aisne_gpd.to_file(aisne_editfile)

#Merge into one layer
sel_cols = ["type_stand","type_aux","regime","regime2", "nat_id", "nat_id2", "orig_mo","orig_mo2", "orig_mo3",
            "date_id","fid","id", "id_loc", "code_hydro","orig_layer", "geometry"
            ]
if not Path(out_net_path).exists() or overwrite:
    print('Merging all layers')
    net_editlist = getfilelist(out_dir, '.*_edit[.]gpkg$')
    net_gpdlist = []
    for net_path in net_editlist:
        print(net_path)
        net_gpdlist.append(gpd.read_file(net_path).assign(orig_layer=Path(net_path).parent.stem))

    net_merged = pd.concat(
        net_path_sub.loc[:,net_path_sub.columns.isin(sel_cols)]
        for net_path_sub in net_gpdlist).\
        pipe(gpd.GeoDataFrame)
    net_merged.to_file(Path(out_net_path))

#Create stable UID - #######################find way to integrate in open source code ###################################
ce_net = os.path.join(resdir, "carto_loi_eau_france.gpkg", "main.carto_loi_eau_france")
if 'UID_CE' not in [f.name for f in arcpy.ListFields(ce_net)]:
    arcpy.AddField_management(in_table=ce_net, field_name='UID_CE', field_type='LONG')
    with arcpy.da.UpdateCursor(ce_net, ['UID_CE', 'OID@']) as cursor:
        for row in cursor:
            row[0] = row[1]
            cursor.updateRow(row)

