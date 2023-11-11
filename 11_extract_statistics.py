from setup_classement import  *

anci_dir = os.path.join(datdir, 'données_auxiliaires')#Ancillary data directory
pregdb = os.path.join(resdir, 'preprocessing_ancillary_data.gdb')

cats_hybasdeps = os.path.join(pregdb, 'BV_hybas0809_depsinters') #BV joined to HydroBASINS and intersected with Departements

#Compiled DDT networks
ce_net = os.path.join(resdir, "carto_loi_eau_france.gpkg", "main.carto_loi_eau_france")
#Pourpoints of compiled DDT networks
ce_net_prpts = os.path.join(resdir, "carto_loi_eau_france.gpkg", "main.carto_loi_eau_france_prpts")

#-------------------------------- INPUTS -------------------------------------------------------------------------------
#Land cover - theia oso
lc_dir = os.path.join(anci_dir, 'oso')
lc_filedict = {yr: getfilelist(os.path.join(lc_dir, "oso_{}".format(yr)),"Classif_Seed_.*[.]tif$")
               for yr in [2018, 2019, 2020, 2021]}

#AWC
gsm_awc_dir = os.path.join(anci_dir, "gsm")
gsm_awc_filelist = getfilelist(gsm_awc_dir, 'awc_mm_.*[.]tif$')

#-------------------------------- OUTPUTS ------------------------------------------------------------------------------
statsgdb = os.path.join(resdir, 'env_stats.gdb')
if not arcpy.Exists(statsgdb):
    arcpy.CreateFileGDB_management(out_folder_path=resdir, out_name=os.path.split(statsgdb)[1])

#-------------------------------- ANALYSIS -----------------------------------------------------------------------------
#Rasterize compiled DDT networks at 25 m resolution
ce_net_ras25m = os.path.join(pregdb, 'carto_loi_eau_france_ras25m')

#Rasterize cats_hybasdeps at the resolution of land cover data
#Create stable UID
# if 'UID_BV' not in [f.name for f in arcpy.ListFields(cats_hybasdeps)]:
#     arcpy.AddField_management(in_table=cats_hybasdeps, field_name='UID_BV', field_type='LONG')
#     with arcpy.da.UpdateCursor(cats_hybasdeps, ['UID_BV', 'OID@']) as cursor:
#         for row in cursor:
#             row[0] = row[1]
#             cursor.updateRow(row)

#USe 'UID_BV'

#Buffer 10 m on each side (only laterally) and rasterize buffer for streams at 10 m


#Analyze ONDE-----------------------------------------------------------------------------------------------------------
#statistics abbrevation: data source, spatial extent, statistics
{'slo_dg_sav': ['zonal', os.path.join(pregdb, "bdalti_25m_slope"), cats_hybasdeps, 'AVERAGE']
,'slo_dg_uav': ['cell', os.path.join(pregdb, "bdalti_25m_slope"), ce_net_prpts, 'AVERAGE']
,'elv_mt_rmx': ['zonal', os.path.join(pregdb, "bdalti_25m_preconditioned"), ce_net_ras25m, 'MAXIMUM'] #To compute stream gradient
,'elv_mt_rmn': ['zonal', os.path.join(pregdb, "bdalti_25m_preconditioned"), ce_net_ras25m, 'MINIMUM'] #To compute stream gradient
,'tcv_ix_rav': ['zonal', os.path.join(pregdb, "bdalti_25m_curvature_tangential"), ce_net_ras25m, 'AVERAGE']
,'tcv_ix_uav': ['cell', os.path.join(pregdb, "bdalti_25m_curvature_tangential"), ce_net_prpts, 'AVERAGE']
,'pcv_ix_rav': ['zonal', os.path.join(pregdb, "bdalti_25m_curvature_profile"), ce_net_ras25m, 'AVERAGE']
,'pcv_ix_uav': ['cell', os.path.join(pregdb, "bdalti_25m_curvature_profile"), ce_net_prpts, 'AVERAGE']
,'ari_ix_syr': ['zonal', os.path.join(pregdb, "ai_v3_yrav"), cats_hybasdeps, 'AVERAGE']
,'ari_ix_ssu': ['zonal', os.path.join(pregdb, "ai_v3_summerav"), cats_hybasdeps, 'AVERAGE']
,'ari_ix_uyr': ['cell', os.path.join(pregdb, "ai_v3_yrav"), ce_net_prpts, 'AVERAGE']
,'ari_ix_usu': ['cell', os.path.join(pregdb, "ai_v3_summerav"), ce_net_prpts, 'AVERAGE']
,'lc_pc_s18': ['tabulate', lc_filedict[2018], cats_hybasdeps]
,'lc_pc_s19': ['tabulate', lc_filedict[2018], cats_hybasdeps]
,'lc_pc_s20': ['tabulate', lc_filedict[2018], cats_hybasdeps]
,'lc_pc_s21': ['tabulate', lc_filedict[2018], cats_hybasdeps]
,'lc_pc_b18': []
,'lc_pc_b19': []
,'lc_pc_b20': []
,'lc_pc_b21': []
,'veg_pc_use': []
,'wet_pc_use': []
,'gla_pc_use': []
,'luc_cl_umj': []
,'urb_pc_use': []
,'agr_pc_use': []
,'crp_pc_use': []
,'scr_pc_use': []
,'vny_pc_use': []
,'awc_mm_sav': []
,'awc_mm_uav': []
,'cly_pc_sav': []
,'cly_pc_uav': []
,'slt_pc_sav': []
,'slt_pc_uav': []
,'snd_pc_sav': []
,'snd_pc_uav': []
,'ppd_pk_sav': []
,'ppd_pk_uav': []
}
