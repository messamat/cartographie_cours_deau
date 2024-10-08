from setup_classement import  *

import os, sys, re
import requests
from urllib.parse import urljoin
from bs4 import BeautifulSoup


ddtsites_datdir = Path(datdir, 'sitesDDT_metadata')

if not ddtsites_datdir.exists():
    ddtsites_datdir.mkdir()

#Read compilation of data sources by French department
datasources_path = list(datdir.glob('metadonnes_cartographie_cours_deau*xlsx'))[-1] #Get most recent table
datasources_pd = pd.read_excel(datasources_path)

def savePage(url, pagepath='page'):
    def savenRename(soup, pagefolder, session, url, tag, inner):
        if not os.path.exists(pagefolder): # create only once
            os.mkdir(pagefolder)
        for res in soup.findAll(tag):   # images, css, etc..
            if res.has_attr(inner): # check inner tag (file object) MUST exists
                try:
                    filename, ext = os.path.splitext(os.path.basename(res[inner])) # get name and extension
                    filename = re.sub('\W+', '', filename) + ext # clean special chars from name
                    fileurl = urljoin(url, res.get(inner))
                    filepath = os.path.join(pagefolder, filename)
                    # rename html ref so can move html and folder of files anywhere
                    res[inner] = os.path.join(os.path.basename(pagefolder), filename)
                    if not os.path.isfile(filepath): # was not downloaded
                        with open(filepath, 'wb') as file:
                            filebin = session.get(fileurl)
                            file.write(filebin.content)
                except Exception as exc:
                    print(exc, file=sys.stderr)
    session = requests.Session()
    #... whatever other requests config you need here
    response = session.get(url)
    soup = BeautifulSoup(response.text, "html.parser")
    path, _ = os.path.splitext(pagepath)
    pagefolder = path+'_files' # page contents folder
    tags_inner = {'img': 'src', 'link': 'href', 'script': 'src'} # tag&inner tags to grab
    for tag, inner in tags_inner.items(): # saves resource files and rename refs
        savenRename(soup, pagefolder, session, url, tag, inner)
    with open(path+'.html', 'wb') as file: # saves modified html doc
        file.write(soup.prettify('utf-8'))


def save_DDT_page(row, in_dir):
    print(f'Processing {row.Département} - {row.Numéro}')

    depdir = Path(in_dir, re.sub('[-]', '_', row.Département))
    if not depdir.exists():
        os.mkdir(depdir)

    depdocdir = Path(depdir, 'documents')
    if not depdocdir.exists():
        os.mkdir(depdocdir)

    savePage(row['Site préfecture 1'], Path(depdir, 'siteDDT'))

datasources_pd.apply(save_DDT_page,
                     args=(ddtsites_datdir,),
                     axis=1)


