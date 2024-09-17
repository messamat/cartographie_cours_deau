import arcpy
from setup_classement import *

overwrite=True

outputs_gdb = os.path.join(resdir, 'analysis_outputs.gdb')
pregdb = os.path.join(resdir, "preprocessing_ancillary_data.gdb")

#Inputs
env_dd_merged_bv_tab = os.path.join(resdir, 'env_dd_merged_bv.csv')
bv_deps_inters = os.path.join(pregdb, 'BV_hybas0809_depsinters')

#Outputs
bv_depinters_env_dd = os.path.join(outputs_gdb, 'bv_depinters_env_dd')

[f.name for f in arcpy.ListFields(bv_depinters_env_dd)]
[f.name for f in arcpy.ListFields(env_dd_merged_bv_tab)]
#Analysis
if (not arcpy.Exists(bv_depinters_env_dd)) or overwrite:
    arcpy.MakeFeatureLayer_management(bv_deps_inters, 'bd_deps_layer')
    arcpy.AddJoin_management('bd_deps_layer',
                             in_field='UID_BV',
                             join_table=env_dd_merged_bv_tab,
                             join_field='UID_BV'
                             )
    arcpy.CopyFeatures_management('bd_deps_layer', bv_depinters_env_dd)
