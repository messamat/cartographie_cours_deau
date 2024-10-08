import arcpy.management

from setup_classement import  *

anci_dir = os.path.join(datdir, 'données_auxiliaires')#Ancillary data directory
pregdb = os.path.join(resdir, 'preprocessing_ancillary_data.gdb')
tempgdb = os.path.join(resdir, "scratch.gdb")

cats_hybasdeps = os.path.join(pregdb, 'BV_hybas0809_depsinters') #BV joined to HydroBASINS and intersected with Departements

#Compiled DDT networks
ce_net = os.path.join(resdir, "carto_loi_eau_france.gpkg", "main.carto_loi_eau_france")
#Pourpoints of compiled DDT networks
ce_net_prpts = os.path.join(resdir, "carto_loi_eau_france.gpkg", "main.carto_loi_eau_france_prpts")

#-------------------------------- INPUTS -------------------------------------------------------------------------------
#Land cover - theia oso
lc_dir = os.path.join(anci_dir, 'oso')
lc_filedict = {yr: getfilelist(os.path.join(lc_dir, "oso_{}".format(yr)),"Classif_Seed_.*[.]tif$")[0]
               for yr in [2019, 2020, 2021]}
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
#USe 'UID_BV'

#Buffer 10 m on each side (only laterally) and rasterize buffer for streams at 10 m
ce_buf = os.path.join(pregdb, 'carto_loi_eau_france_buf10m')
ce_buf_ras10m = os.path.join(pregdb, 'carto_loi_eau_france_bufras10m')

#Project cats_hybasdeps to extract population counts
cats_hybasdeps_LAEA = os.path.join(tempgdb, 'cats_hybasdeps_LAEA')
pop_ras_200m = os.path.join(pregdb, 'insee_pop_interp200m')
if not arcpy.Exists(cats_hybasdeps_LAEA):
    arcpy.management.Project(in_dataset=cats_hybasdeps,
                             out_dataset=cats_hybasdeps_LAEA,
                             out_coor_system=pop_ras_200m)

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
    , 'lc_pc_s19': ['tabulate', lc_filedict[2019], cats_hybasdeps]
    , 'lc_pc_s20': ['tabulate', lc_filedict[2020], cats_hybasdeps]
    , 'lc_pc_s21': ['tabulate', lc_filedict[2021], cats_hybasdeps]
    , 'lc_pc_b19': ['tabulate', lc_filedict[2019], ce_buf_ras10m]
    , 'lc_pc_b20': ['tabulate', lc_filedict[2020], ce_buf_ras10m]
    , 'lc_pc_b21': ['tabulate', lc_filedict[2021], ce_buf_ras10m]
    ,'veg_pc_use': ['cell', os.path.join(lcav_dir, 'oso_veg_acc'), ce_net_prpts, 'MEAN']
    ,'wet_pc_use': ['cell', os.path.join(lcav_dir, 'oso_wet_acc'), ce_net_prpts, 'MEAN']
    ,'gla_pc_use': ['cell', os.path.join(lcav_dir, 'oso_gla_acc'), ce_net_prpts, 'MEAN']
    ,'imp_pc_use': ['cell', os.path.join(lcav_dir, 'oso_imp_acc'), ce_net_prpts, 'MEAN']
    ,'agr_pc_use': ['cell', os.path.join(lcav_dir, 'oso_agr_acc'), ce_net_prpts, 'MEAN']
    ,'crp_pc_use': ['cell', os.path.join(lcav_dir, 'oso_crp_acc'), ce_net_prpts, 'MEAN']
    ,'scr_pc_use': ['cell', os.path.join(lcav_dir, 'oso_scr_acc'), ce_net_prpts, 'MEAN']
    ,'vny_pc_use': ['cell', os.path.join(lcav_dir, 'oso_vny_acc'), ce_net_prpts, 'MEAN']
    ,'veg_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_veg'), cats_hybasdeps, 'MEAN']
    ,'for_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_for'), cats_hybasdeps, 'MEAN']
    ,'wet_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_wet'), cats_hybasdeps, 'MEAN']
    ,'gla_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_cl22.tif'), cats_hybasdeps, 'MEAN']
    ,'rck_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_cl20.tif'), cats_hybasdeps, 'MEAN']
    ,'imp_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_imp'), cats_hybasdeps, 'MEAN']
    ,'agr_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_agr'), cats_hybasdeps, 'MEAN']
    ,'crp_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_crp'), cats_hybasdeps, 'MEAN']
    ,'scr_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_scr'), cats_hybasdeps, 'MEAN']
    ,'vny_pc_sse': ['zonal', os.path.join(lcav_dir, 'oso_cl15.tif'), cats_hybasdeps, 'MEAN']
    ,'ppc_in_sav': ['zonal', os.path.join(pregdb, 'insee_pop_interp200m'), cats_hybasdeps_LAEA, 'SUM']
    ,'ppc_in_uav': ['cell', os.path.join(pregdb, 'insee_pop_interp200m_acc'), cats_hybasdeps_LAEA, 'SUM']
}

lc_class_dict = {1:'urba1', 2:'urba2', 3:'indus', 4:'roads', 5:'wioil', 6:'straw', 7:'spoil', 8:"soy", 9:"sunfl",
              10:"corn", 11:"rice", 12:"roots", 13:"pastu", 14:"orchd", 15:"vinyd", 16:"forbr", 17:"forco",
              18:"grass", 19:"heath", 20:"rocks", 21:"beach", 22:'glasn', 23:'water'}

horizon_lims = [0, 5, 15, 30, 60, 100, 200]
horizons = ["{0}_{1}".format(horizon_lims[:-1][i], horizon_lims[1:][i])
            for i in range(len(horizon_lims)-1)]
for ho in horizons:
    config_dict['awc_mm_sav{}'.format(ho)] = ['zonal', os.path.join(gsm_dir, "awc_mm_{}.tif".format(ho)),
                                              cats_hybasdeps, 'MEAN']
    config_dict['cly_pc_sav{}'.format(ho)] = ['zonal', os.path.join(gsm_dir, "argile.{}.tif".format(ho)),
                                              cats_hybasdeps, 'MEAN']
    config_dict['slt_pc_sav{}'.format(ho)] = ['zonal', os.path.join(gsm_dir, "limon.{}.tif".format(ho)),
                                              cats_hybasdeps, 'MEAN']
    config_dict['snd_pc_sav{}'.format(ho)] = ['zonal', os.path.join(gsm_dir, "sable.{}.tif".format(ho)),
                                              cats_hybasdeps, 'MEAN']
    config_dict['awc_mm_uav{}'.format(ho)] = ['cell', os.path.join(gsm_dir, "awc_mm_{}_acc.tif".format(ho)),
                                              ce_net_prpts, 'MEAN']
    config_dict['cly_pc_uav{}'.format(ho)] = ['cell', os.path.join(gsm_dir, "argile.{}_acc.tif".format(ho)),
                                              ce_net_prpts, 'MEAN']
    config_dict['slt_pc_uav{}'.format(ho)] = ['cell', os.path.join(gsm_dir, "limon.{}_acc.tif".format(ho)),
                                              ce_net_prpts, 'MEAN']
    config_dict['snd_pc_uav{}'.format(ho)] = ['cell', os.path.join(gsm_dir, "sable.{}_acc.tif".format(ho)),
                                              ce_net_prpts, 'MEAN']

for k,v in config_dict.items():
    #print(k)
    if v[0] == 'zonal':
        if arcpy.Exists(v[1]):
            if arcpy.Exists(v[2]):
                out_tab = os.path.join(statsgdb, k)
                if not arcpy.Exists(out_tab):
                    print("Processing {}...".format(k))
                    ZonalStatisticsAsTable(in_zone_data=v[2], zone_field='UID_BV', in_value_raster=v[1],
                                           out_table=out_tab, statistics_type=v[3])
                else:
                    print("{} already exists. Skipping...".format(out_tab))

    elif v[0] == 'tabulate':
        if arcpy.Exists(v[1]):
            if arcpy.Exists(v[2]):
                out_tab = os.path.join(statsgdb, k)
                if not arcpy.Exists(out_tab):
                    print("Processing {}...".format(k))
                    TabulateArea(in_zone_data=v[2], zone_field='UID_BV', in_class_data=v[1], class_field='Value',
                                 out_table=out_tab, classes_as_rows='CLASSES_AS_FIELDS')
                else:
                    print("{} already exists. Skipping...".format(out_tab))

############### extra stuff
#for yr in [2018, 2019, 2020, 2021]:
#     lcyr_gdb = os.path.join(lcav_dir, 'lc{}_tiles.gdb'.format(yr))
#     for tile in lcyr_gdb:
#         tile_num = re.sub('lc{}_'.format(yr), '', tile)
#         config_dict['lc_pc_s{0}_{1}'.format(yr, tile_num)] = ['tabulate', tile, cats_hybasdeps]
#         config_dict['lc_pc_s{0}_{1}'.format(yr, tile_num)] = ['tabulate', tile, cats_hybasdeps]