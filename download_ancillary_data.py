from classement_setup import  *


admin_express_url = "https://wxs.ign.fr/x02uy2aiwjo9bm8ce5plwqmr/telechargement/prepackage/ADMINEXPRESS_SHP_WGS84G_PACK_2023-05-04$ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03/file/ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03.7z"
admin_express_zip_path = os.path.join(datdir, 'ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03.7z')
admin_express_unzipped_path = os.path.join(datdir, 'ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03')
if not os.path.exists(admin_express_unzipped_path):
    with open(admin_express_zip_path, "wb") as file:
        # get request
        print(f"Downloading {Path(admin_express_url).name}")
        response = requests.get(admin_express_url, verify=False)
        file.write(response.content)