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
lc_filedict = {yr: getfilelist(os.path.join(lc_dir, "oso_{}".format(yr)),"Classif_Seed_.*[.]tif$")[0]
               for yr in [2018, 2019, 2020, 2021]}
lcav_dir = os.path.join(resdir, 'oso_lc_stats')

#AWC
gsm_dir = os.path.join(anci_dir, "gsm")

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
ce_buf = os.path.join(pregdb, 'carto_loi_eau_france_buf10m')
ce_buf_ras10m = os.path.join(pregdb, 'carto_loi_eau_france_bufras10m')

#Accumulate predictors-----------------------------------------------------------------------------------------------------------

os.path.join(lcav_dir, 'oso_veg')
os.path.join(lcav_dir, 'oso_veg_acc')

#Analyze ONDE-----------------------------------------------------------------------------------------------------------


#statistics abbrevation: data source, spatial extent, statistics
config_dict = {
    'slo_dg_sav': ['zonal', os.path.join(pregdb, "bdalti_25m_slope"), cats_hybasdeps, 'MEAN']
    ,'slo_dg_uav': ['cell', os.path.join(pregdb, "bdalti_25m_slope_acc"), ce_net_prpts, 'MEAN']
    ,'elv_mt_rmx': ['zonal', os.path.join(pregdb, "bdalti_25m_preconditioned"), ce_net_ras25m, 'MAXIMUM'] #To compute stream gradient
    ,'elv_mt_rmn': ['zonal', os.path.join(pregdb, "bdalti_25m_preconditioned"), ce_net_ras25m, 'MINIMUM'] #To compute stream gradient
    ,'tcv_ix_rav': ['zonal', os.path.join(pregdb, "bdalti_25m_curvature_tangential"), ce_net_ras25m, 'MEAN']
    ,'tcv_ix_uav': ['cell', os.path.join(pregdb, "bdalti_25m_curvature_tangential_acc"), ce_net_prpts, 'MEAN']
    ,'pcv_ix_rav': ['zonal', os.path.join(pregdb, "bdalti_25m_curvature_profile"), ce_net_ras25m, 'MEAN']
    ,'pcv_ix_uav': ['cell', os.path.join(pregdb, "bdalti_25m_curvature_profile_acc"), ce_net_prpts, 'MEAN']
    ,'ari_ix_syr': ['zonal', os.path.join(pregdb, "ai_v3_yrav"), cats_hybasdeps, 'MEAN']
    ,'ari_ix_ssu': ['zonal', os.path.join(pregdb, "ai_v3_summerav"), cats_hybasdeps, 'MEAN']
    ,'ari_ix_uyr': ['cell', os.path.join(pregdb, "ai_v3_yrav_acc"), ce_net_prpts, 'MEAN']
    ,'ari_ix_usu': ['cell', os.path.join(pregdb, "ai_v3_summerav_acc"), ce_net_prpts, 'MEAN']
    ,'lc_pc_s18': ['tabulate', lc_filedict[2018], cats_hybasdeps]
    ,'lc_pc_s19': ['tabulate', lc_filedict[2019], cats_hybasdeps]
    ,'lc_pc_s20': ['tabulate', lc_filedict[2020], cats_hybasdeps]
    ,'lc_pc_s21': ['tabulate', lc_filedict[2021], cats_hybasdeps]
    ,'lc_pc_b18': ['tabulate', lc_filedict[2018], ce_buf_ras10m]
    ,'lc_pc_b19': ['tabulate', lc_filedict[2019], ce_buf_ras10m]
    ,'lc_pc_b20': ['tabulate', lc_filedict[2020], ce_buf_ras10m]
    ,'lc_pc_b21': ['tabulate', lc_filedict[2021], ce_buf_ras10m]
    ,'veg_pc_use': ['cell', os.path.join(lcav_dir, 'oso_veg_acc'), ce_net_prpts, 'MEAN']
    ,'wet_pc_use': ['cell', os.path.join(lcav_dir, 'oso_wet_acc'), ce_net_prpts, 'MEAN']
    ,'gla_pc_use': ['cell', os.path.join(lcav_dir, 'oso_gla_acc'), ce_net_prpts, 'MEAN']
    ,'imp_pc_use': ['cell', os.path.join(lcav_dir, 'oso_imp_acc'), ce_net_prpts, 'MEAN']
    ,'agr_pc_use': ['cell', os.path.join(lcav_dir, 'oso_agr_acc'), ce_net_prpts, 'MEAN']
    ,'crp_pc_use': ['cell', os.path.join(lcav_dir, 'oso_crp_acc'), ce_net_prpts, 'MEAN']
    ,'scr_pc_use': ['cell', os.path.join(lcav_dir, 'oso_scr_acc'), ce_net_prpts, 'MEAN']
    ,'vny_pc_use': ['cell', os.path.join(lcav_dir, 'oso_vny_acc'), ce_net_prpts, 'MEAN']
    ,'ppd_pk_sav': ['zonal', os.path.join(pregdb, "ppd_pk"), cats_hybasdeps, 'MEAN']
    ,'ppd_pk_uav': ['cell', os.path.join(pregdb, "ppd_pk_acc"), ce_net_prpts, 'MEAN']
}

horizon_lims = [0, 5, 15, 30, 60, 100, 200]
horizons = ["{0}_{1}".format(horizon_lims[:-1][i], horizon_lims[1:][i])
            for i in range(len(horizon_lims)-1)]
for ho in horizons:
    config_dict['awc_mm_sav{}'.format(ho)] = ['zonal', os.path.join(gsm_dir, "awc_mm_{}.tif".format(ho)),
                                              cats_hybasdeps, 'MEAN']
    config_dict['cly_pc_sav{}'.format(ho)] = ['zonal', os.path.join(gsm_dir, "argile.{}.tif".format(ho)),
                                              cats_hybasdeps, 'MEAN']
    config_dict['slt_pc_sav{}'.format(ho)] = ['zonal', os.path.join(gsm_dir, "limon{}.tif".format(ho)),
                                              cats_hybasdeps, 'MEAN']
    config_dict['snd_pc_sav{}'.format(ho)] = ['zonal', os.path.join(gsm_dir, "sable{}.tif".format(ho)),
                                              cats_hybasdeps, 'MEAN']
    config_dict['awc_mm_uav{}'.format(ho)] = ['cell', os.path.join(gsm_dir, "awc_mm_{}_acc.tif".format(ho)),
                                              ce_net_prpts, 'MEAN']
    config_dict['cly_pc_uav{}'.format(ho)] = ['cell', os.path.join(gsm_dir, "argile.{}_acc.tif".format(ho)),
                                              ce_net_prpts, 'MEAN']
    config_dict['slt_pc_uav{}'.format(ho)] = ['cell', os.path.join(gsm_dir, "limon{}_acc.tif".format(ho)),
                                              ce_net_prpts, 'MEAN']
    config_dict['snd_pc_uav{}'.format(ho)] = ['cell', os.path.join(gsm_dir, "sable{}_acc.tif".format(ho)),
                                              ce_net_prpts, 'MEAN']


for k,v in config_dict.items():
    if v[0] == 'zonal':
        if arcpy.Exists(v[1]):
            if arcpy.Exists(v[2]):
                print("Processing {}...".format(k))
                ZonalStatisticsAsTable(in_zone_data=v[2], zone_field='UID_BV', in_value_raster=v[1],
                                       out_table=os.path.join(statsgdb, k), statistics_type=v[3])

    elif v[0] == 'tabulate':
        if arcpy.Exists(v[1]):
            if arcpy.Exists(v[2]):
                print("Processing {}...".format(k))
                TabulateArea(in_zone_data=v[2], zone_field='UID_BV', in_class_data=v[1], class_field='Value',
                             out_table=os.path.join(statsgdb, k), classes_as_rows='CLASSES_AS_FIELDS')
