from classement_setup import * #Get directory structure and packages
import geopandas as gpd
from datetime import date

#Folder where geo-ide are stored
geoide_dir = Path(datdir, 'geoide')
other_dir = Path(datdir, 'autres_téléchargements')

#List datasets
lyrs_geoide = getfilelist(geoide_dir, '.+[.](shp)|(tab)|(TAB)$')
lyrs_other = getfilelist(other_dir, '.+[.](shp)|(tab)|(TAB)$')

lyrs = lyrs_geoide + lyrs_other

lyrmetadata_dict = {}

for i in lyrs:
    dept_name = re.sub(r".*[\\]D[0-9]{1,2}_", "", i).split(os.path.sep)[0]
    print(dept_name)
    try:
        gdf = gpd.read_file(i) #Errors with some layers that contain Linestrings that are built on less than three coordinates or points
        lyrmetadata_dict[i] = [dept_name, len(gdf), ', '.join(list(gdf.columns)[:-1])]

    # Return any other type of error
    except:
        # By default any other errors will be caught here
        e = sys.exc_info()[1]
        print(e.args[0])

pd.DataFrame.from_dict(lyrmetadata_dict, orient='index').to_csv(
    Path(datdir, f"lyrsmetadata_{date.today().strftime('%Y%m%d')}.csv")
)