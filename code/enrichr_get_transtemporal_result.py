	
import json
import requests


ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/export'
query_string = '?userListId=%s&filename=%s&backgroundType=%s'
user_list_id = 56767500
filename = 'transtemporal_DEGs'
gene_set_library = 'Tabula_Muris'

url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
response = requests.get(url, stream=True)

with open('./outputs/'+filename + '.txt', 'wb') as f:
    for chunk in response.iter_content(chunk_size=1024): 
        if chunk:
            f.write(chunk)