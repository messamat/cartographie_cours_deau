import requests
from lxml import html
from googlesearch import search
from bs4 import BeautifulSoup
from setup_classement import  *


geoide_dir = Path(resdir, 'geoide')

#Read compilation of data sources by French department
datasources_path = list(datdir.glob('metadonnes_cartographie_cours_deau*xlsx'))[-1] #Get most recent table
datasources_pd = pd.read_excel(datasources_path)

dept_carto_urls = {}
for d in datasources_pd.Département:
    print(d)
    query = f"Cartographie des cours d'eau {d}"
    ## Google Search query results as a Python List of URLs
    dept_carto_urls[d] = list(search(query, tld="co.in", num=10, stop=3, pause=1))[0]


for k, v in dept_carto_urls.items():
    print(f"{k} : {v}")