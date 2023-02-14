import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) #TO GET READ OF ITERITEMS SPAM FROM PANDAS

from classement_setup import * #Get directory structure and packages

#Create folder to store geo-ide
geoide_dir = Path(resdir, 'geoide')
if not geoide_dir.exists():
    os.mkdir(geoide_dir)

#Read compilation of data sources by French department
datasources_path = list(datdir.glob('metadonnes_cartographie_cours_deau*xlsx'))[-1] #Get most recent table
datasources_pd = pd.read_excel(datasources_path)

#Format as long table
re_geoide = re.compile('Catalogue interministériel de données géographiques')
geoide_cols = list(filter(re_geoide.match, datasources_pd.columns))
datasources_pdlong = pd.melt(datasources_pd,
                         id_vars=list(itertools.filterfalse(re_geoide.match, datasources_pd.columns)),
                         value_vars=geoide_cols,
                         var_name='sitenum',
                         value_name='geoide_url',
                         ignore_index=True).\
    dropna(subset='geoide_url').\
    reset_index(drop=True)

#Format base URL to URL to XML metadata
datasources_pdlong['metadata_xml_url'] = datasources_pdlong['geoide_url'].apply(
    lambda geoide_url: f"http://catalogue.geo-ide.developpement-durable.gouv.fr/catalogue/srv/api/records/{geoide_url.rsplit('/', 1)[-1]}/formatters/xml"
)

#Alternative source
# xml_urlroot = "https://ogc.geo-ide.developpement-durable.gouv.fr/csw/all-dataset?REQUEST=GetRecordById&SERVICE=CSW&VERSION=2.0.2&RESULTTYPE=results&elementSetName=full&TYPENAMES=gmd:MD_Metadata&OUTPUTSCHEMA=http://www.isotc211.org/2005/gmd&ID="
# datasources_pdlong['metadata_xml_url'] = datasources_pdlong['geoide_url'].apply(
#     lambda x: f"{xml_urlroot}{x.rsplit('/', 1)[-1]}"
# )

#Parse/Extract metadata from XML
metadata_dict_list = datasources_pdlong.apply(get_geoide_metadata_tab, args=(geoide_dir, False), axis=1)

metadata_xml_pd = pd.concat(
    [pd.DataFrame.from_dict(d, orient='index').transpose() for d in metadata_dict_list if isinstance(d, dict)],
    axis=0
)

#Write table
metadata_xml_pd.\
    sort_values(by=['Numéro', 'revision_date']).\
    to_csv(Path(resdir, "metadonnes_cartographie_cours_deau_geoidemerge.csv"), encoding='utf-8')

#~~~~~~~~~~~~~~~~~~~~~~~~~ DOWNLOAD DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Download atom archive
metadata_xml_pd['Atom_archive_download_output'] = metadata_xml_pd.apply(download_atomarchive_fromgeoidemetadata,
                                                                        args=(geoide_dir,),
                                                                        axis=1)
#Check which departments have no atom archive
metadata_xml_pd[~metadata_xml_pd['Département'].isin(metadata_xml_pd[(metadata_xml_pd['Atom_archive_download_output'] != 'No atom archive')]['Département'].unique())]['Département'].unique()

#Unzip downloaded data
for zipf in getfilelist(dir=Path(resdir, 'geoide'), repattern=".*[.]zip"):
    unzip(zipf)

#Download WFS data
#For charentes
#https://ogc.geo-ide.developpement-durable.gouv.fr/wxs?map=/opt/data/carto/geoide-catalogue/1.4/org_37972/bef85666-f513-4b4f-9f9a-051f6e157674.internet.map





