#To install new module in venv: .\pyenv\python.exe -m pip install

import chardet
from datetime import date
from collections import defaultdict
import geopandas as gpd
from inspect import getsourcefile
import itertools
import numpy as np
import os
from osgeo import gdal, osr, ogr #
import pandas as pd #Rpd.read_excel equires optional dependency openpyxl
from pathlib import Path
#import pip_system_certs.wrapt_requests
#import pygeoops #problem with numpy
import re
import requests
import shapefile
import shutil
import sys
import urllib3
import xmltodict
from flatten_dict import flatten #INstall from https://github.com/ianlini/flatten-dict #github connection solved with https://stackoverflow.com/questions/72486457/fatal-unable-to-connect-to-github-com-github-com0-140-82-121-4-errno-unkno #Then .\pyenv\python.exe -m pip install git+git://github.com/ianlini/flatten-dict.git
import zipfile

import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = False
sys.stdout.reconfigure(encoding='utf-8')

#Utility functions
def get_root_fromsrcdir():
    return(os.path.dirname(os.path.abspath(
        getsourcefile(lambda:0)))).split('\\src')[0]

#Folder structure~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rootdir = get_root_fromsrcdir()
datdir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

if not os.path.exists(datdir):
    os.mkdir(datdir)
if not os.path.exists(resdir):
    os.mkdir(resdir)

#Get all files in a ArcGIS workspace (file or personal GDB)
def getwkspfiles(dir, repattern=None):
    arcpy.env.workspace = dir
    filenames_list = (arcpy.ListDatasets() or []) +\
                     (arcpy.ListTables() or []) +\
                     (arcpy.ListFeatureClasses() or []) # Either LisDatsets or ListTables may return None so need to create empty list alternative
    if not repattern == None:
        filenames_list = [filen for filen in filenames_list if re.search(repattern, filen)]

    return ([os.path.join(dir, f) for f in filenames_list])
    arcpy.ClearEnvironment('workspace')

def getfilelist(dir, repattern=None, gdbf=True, nongdbf=True, fullpath=False):
    """Function to iteratively go through all subdirectories inside 'dir' path
    and retrieve path for each file that matches "repattern"
    gdbf and nongdbf allows the user to choose whether to consider ArcGIS workspaces (GDBs) or not or exclusively"""

    try:
        if arcpy.Describe(dir).dataType == 'Workspace':
            if gdbf == True:
                print('{} is ArcGIS workspace...'.format(dir))
                filenames_list = getwkspfiles(dir, repattern)
            else:
                raise ValueError(
                    "A gdb workspace was given for dir but gdbf=False... either change dir or set gdbf to True")
        else:
            filenames_list = []

            if gdbf == True:
                for (dirpath, dirnames, filenames) in os.walk(dir):
                    for in_dir in dirnames:
                        fpath = os.path.join(dirpath, in_dir)
                        if arcpy.Describe(fpath).dataType == 'Workspace':
                            print('{} is ArcGIS workspace...'.format(fpath))
                            filenames_list.extend(getwkspfiles(dir=fpath, repattern=repattern))
            if nongdbf == True:
                for (dirpath, dirnames, filenames) in os.walk(dir):
                    for file in filenames:
                        if repattern is None:
                            filenames_list.append(os.path.join(dirpath, file))
                        else:
                            if re.search(repattern, file):
                                filenames_list.append(os.path.join(dirpath, file))
        return (filenames_list)

    # Return geoprocessing specific errors
    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages(2))
    # Return any other type of error
    except:
        # By default any other errors will be caught here
        e = sys.exc_info()[1]
        print(e.args[0])

def divbb(bbox, res, divratio):
    box_lc_x, box_lc_y, box_rc_x, box_rc_y = bbox
    coln = (box_rc_x - box_lc_x) / float(res)
    rown = (box_rc_y - box_lc_y) / float(res)

    xbblist = np.arange(box_lc_x, box_rc_x + (float(res) * coln / divratio),
                        float(res) * coln / divratio)
    ybblist = np.arange(box_lc_y, box_rc_y + (float(res) * rown / divratio),
                        float(res) * rown / divratio)

    if abs(xbblist[-1]) > abs(box_rc_x):
        xbblist[-1] = box_rc_x

    if abs(ybblist[-1]) > abs(box_rc_y):
        ybblist[-1] = box_rc_y

    xbblist = np.unique(xbblist)
    ybblist = np.unique(ybblist)

    fullbblist = []
    for pairx in zip(xbblist[:-1], xbblist[1:]):
        for pairy in zip(ybblist[:-1], ybblist[1:]):
            fullbblist.append((pairx[0], pairy[0], pairx[1], pairy[1]))

    return (fullbblist)

def unzip(infile):
    # Unzip folder
    if zipfile.is_zipfile(infile):
        print('Unzipping {}...'.format(os.path.split(infile)[1]))
        outdir = Path(os.path.splitext(infile)[0])
        if not outdir.exists():
            outdir.mkdir()

        with zipfile.ZipFile(infile) as zipf:
            zipfilelist = [info.filename for info in zipf.infolist()]
            listcheck = [f for f in zipfilelist if Path(outdir, f).exists()]
            if len(listcheck) > 0:
                print('Overwriting {}...'.format(', '.join(listcheck)))
            for name in zipf.namelist():
                zipf.extract(name, outdir)
        del zipf
    else:
        raise ValueError('Not a zip file')

#Workflow functions
###### MODIFIED FUNCTION FROM https://stackoverflow.com/questions/52081545/python-3-flattening-nested-dictionaries-and-lists-within-dictionaries
# TO MAKE SURE EACH COMPRESSED KEY IS UNIQUE.
def flatten(d, uid=0):
    out = {}
    for key, val in d.items():
        uid += 1
        if isinstance(val, dict):
            val = [val]
        if isinstance(val, list):
            for subdict in val:
                uid += 1
                deeper = flatten(subdict, uid=uid).items()
                out.update({key + str(uid) + '_' + key2 + str(uid): val2 for key2, val2 in deeper})
        else:
            out[key + str(uid) + '_'] = val
    return out

def split_strip(in_record, sep=','):
    if isinstance(in_record.split(sep), list):
        return([cat.strip() for cat in in_record.split(sep)])
    else:
        in_record.strip()

def convert_bytes_to_na(in_gpd):
    for col in np.where(in_gpd.dtypes == 'object')[0]:
        for i in range(len(in_gpd)):
            record =in_gpd.iloc[i, col]
            if not pd.isnull(record):
                if isinstance(record, bytes):
                    #print(col)
                    in_gpd.iloc[i, col] = np.nan
    return(in_gpd)

# Regex search dictionary keys and return values associated with matches
# max_i limits the number of matches
# in_pattern can be either as a simple string or as a list where the first element is the pattern and the second is max_i
def re_search_dict(in_dict, in_pattern, max_i=0):
    if isinstance(in_pattern, list):
        max_i = in_pattern[1]
        in_pattern = in_pattern[0]

    out_list = [in_dict[k] for k in in_dict if re.search(in_pattern, k)]
    if out_list:
        if max_i == 0:
            return (out_list[0])
        else:
            return (out_list[:max_i])

#Export attribute table from an ESRI-compatible feature class to csv
def CopyRows_pd(in_table, out_table, fields_to_copy):
    #Make sure fields_to_copy is a list
    if type(fields_to_copy) == str:
        fields_to_copy = [fields_to_copy]

    if type(fields_to_copy) == dict:
        dict_for_renaming = fields_to_copy
        fields_to_copy = list(fields_to_copy.keys())

    fields_to_copy_valid = []
    intable_flist = [f2.name for f2 in arcpy.ListFields(in_table)]
    for f1 in fields_to_copy:
        if f1 in intable_flist:
            fields_to_copy_valid.append(f1)
        else:
            print("{0} field is not present in {1}".format(f1, in_table))

    rows_to_copy_dict = defaultdict(float)
    with arcpy.da.SearchCursor(in_table, ['OID@']+fields_to_copy_valid) as cursor: #Other fields are badly entered
        for row in cursor:
            rows_to_copy_dict[row[0]] = list(row[1:])

    out_pd = pd.DataFrame.from_dict(data=rows_to_copy_dict, orient='index')
    out_pd.columns = fields_to_copy_valid
    if 'dict_for_renaming' in locals():
        out_pd.rename(columns={k:v for k,v in dict_for_renaming.items() if k in fields_to_copy_valid},
                      inplace=True)
    out_pd.to_csv(out_table, index=False)

def get_geoide_metadata_fromxmldict(in_xmldict):
    flatdict = flatten(in_xmldict['gmd:MD_Metadata']) #If use alternative xml source, need to edit to #flatten(in_xmldict['csw:GetRecordByIdResponse']['gmd:MD_Metadata'])

    pattern_dict = {'id': 'fileIdentifier',
                    "title": ".*MD_DataIdentification.*gmd:citation.*gmd:CI_Citation.*gmd:title",
                    "title_alternate": ".*MD_DataIdentification.*gmd:citation.*gmd:CI_Citation.*gmd:alternateTitle",
                    "title_schema": ".*MD_ApplicationSchemaInformation.*gmd:CI_Citation.*gmd:title",
                    "status": 'gmd:status',
                    'organisation_name': 'organisationName',
                    'org_phone': 'CI_Telephone',
                    'org_address': 'CI_Address',
                    'org_city': 'city',
                    'org_postalcode': 'postalCode',
                    'org_email': 'electronicMailAddress',
                    'org_all': ['CI_ResponsibleParty', 40],
                    'online_resource': 'CI_OnlineResource',
                    "abstract": "gmd:abstract",
                    "purpose": "gmd:purpose",
                    "contact_all": "gmd:pointOfContact",
                    'use_limitations': ['(gmd:useLimitation)|(gmd:otherConstraints)', 40],
                    'keywords': ["gmd:keyword", 10],
                    'dates': ["gmd:CI_Date", 10],
                    'maintenance_code': 'gmd:MD_MaintenanceFrequencyCode',
                    'maintenance_note': 'gmd:maintenanceNote',
                    'lineage': 'gmd:LI_Lineage',
                    "extent": "gmd:EX_Extent",
                    "format": "gmd:distributionFormat",
                    'character_set': 'gmd:MD_CharacterSetCode',
                    'crs': 'referenceSystemIdentifier',
                    'scale': 'gmd:spatialResolution',
                    }

    metadata_dict = {}

    re_search_dict(flatdict, in_pattern='gmd:MD_CharacterSetCode')

    for id, re_pat in pattern_dict.items():
        metadata_dict[id] = re_search_dict(in_dict=flatdict, in_pattern=re_pat)

    transferopts = re_search_dict(flatdict, in_pattern='gmd:transferOptions', max_i=40)
    metadata_dict['atom_internet_url'] = list(
        filter(re.compile('.*atom[.]geo[-]ide[.]developpement[-]durable[.]gouv[.]fr.*').match,
               transferopts))
    metadata_dict['wfs_internet_url'] = list(filter(re.compile(".*internet[.]map$").match,
                                                    transferopts))

    bbox = re_search_dict(flatdict, in_pattern='gmd:EX_GeographicBoundingBox', max_i=10)
    metadata_dict['bbox_west'] = bbox[0]
    metadata_dict['bbox_east'] = bbox[1]
    metadata_dict['bbox_south'] = bbox[2]
    metadata_dict['bbox_north'] = bbox[3]

    if 'revision' in metadata_dict['dates']:
        metadata_dict['revision_date'] = metadata_dict['dates'][metadata_dict['dates'].index('revision')-1] #Get revision date, always before first occurrence of "revision" string in list
    metadata_dict['publication_date'] = metadata_dict['dates'][metadata_dict['dates'].index('publication')-1] #Get revision date, always before first occurrence of "revision" string in list

    return(metadata_dict)

def get_geoide_metadata_tab(row, in_dir, overwrite=False):
    print(f'Processing {row.Département} - {row.sitenum}')

    depdir = Path(in_dir, "D{0}_{1}".format(
        row.Numéro,
        re.sub('[-]', '_', row.Département)
    ))
    if not depdir.exists():
        os.mkdir(depdir)

    exurl = row.metadata_xml_url
    depxml = Path(depdir, f"metadata_{re.sub(' ', '_', row.sitenum)}.xml")

    if ((not depxml.exists()) or (overwrite==True)):
        response = requests.get(exurl)

        with open(depxml, 'wb') as file:
            file.write(response.content)
        xmldict = xmltodict.parse(response.text, encoding='utf-8')
    else:
        with open(depxml, 'rb') as file:
            xmldict = xmltodict.parse(file, encoding='utf-8')

    if 'gmd:MD_Metadata' not in xmldict:
        return({'Numéro': row.Numéro,
                'Département': row.Département,
                'geoide_url': row.geoide_url,
                'xml_url': exurl,
                'id': 'XML link exists but there are no data available'}
        )
    else:
        metadata_dict = {}
        metadata_dict['Numéro'] = row.Numéro
        metadata_dict['Département'] = row.Département
        metadata_dict.update(get_geoide_metadata_fromxmldict(xmldict))
        metadata_dict['geoide_url'] = row.geoide_url
        metadata_dict['xml_url'] = exurl
    return(metadata_dict)

def download_atomarchive_fromgeoidemetadata(row, in_dir, verify_SSLcertificate=True):
    if (not pd.isna(row.atom_internet_url)):
        if (len(row.atom_internet_url) > 0):
            exurl = row.atom_internet_url[0]

            depdir = Path(in_dir, "D{0}_{1}".format(
                row.Numéro,
                re.sub('[-]', '_', row.Département)
            ))

            out_atom = Path(depdir, f"{''.join([x if x.isalnum() else '_' for x in row.title][0:30])}_"
                                    f"{re.sub('[-]', '_', row.id)}.zip")

            #Format SSL certificate verification
            if (verify_SSLcertificate == False):
                in_cert_reqs = 'CERT_NONE'

            http = urllib3.PoolManager(cert_reqs = in_cert_reqs)
            if not out_atom.exists():
                try:
                    print(f'Attempting atom archive download for {row.Département}')
                    f = http.request('GET', exurl, preload_content=False)
                    ftype = f.headers.get('content-type').lower()
                    if ftype in ['application/zip', 'application/x-zip-compressed']:
                        with open(out_atom, 'wb') as out_file:
                            shutil.copyfileobj(f, out_file)
                    else:
                        f'Atom archive is of type {ftype}... did not download'
                    return(out_atom)
                except urllib3.exceptions.HTTPError as e:
                    print('Request failed:', e.reason)
                    return ('Download failed')
            else:
                print(f"{out_atom} already exists. Skipping...")
                return ('Record not in zip format')
        else:
            print('No atom archive')
            return ('No atom archive')
    else:
        print('No atom archive')
        return ('No atom archive')

