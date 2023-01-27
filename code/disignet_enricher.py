import json
import requests
import csv
#### gene_list_submission

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'

#all
genes_str = open('./data/zh.LD_DEG.txt', 'r', encoding="utf-8")

description = 'LD-DEGs'
payload = {
    'list': (None, genes_str),
    'description': (None, description)
}

response = requests.post(ENRICHR_URL, files=payload)
if not response.ok:
    raise Exception('Error analyzing gene list')

data = json.loads(response.text)
print(data)

user_lid= data['userListId']

#########DisGeNET

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'
user_list_id = user_lid
gene_set_library = 'DisGeNET'
response = requests.get(
    ENRICHR_URL + query_string % (user_list_id, gene_set_library)
 )
if not response.ok:
    raise Exception('Error fetching enrichment results')

data = json.loads(response.text)
print(data)

#dowload

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/export'
query_string = '?userListId=%s&filename=%s&backgroundType=%s'
user_list_id = user_lid
filename = 'DisGeNET_enrichment'
gene_set_library = 'DisGeNET'

url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
response = requests.get(url, stream=True)

with open('./data/LD_DEG_'+filename + '.txt', 'wb') as f:
    for chunk in response.iter_content(chunk_size=1024): 
        if chunk:
            f.write(chunk)
            
            
#down
genes_str = open('./data/zh.LD_DEG.dn.txt', 'r', encoding="utf-8")

description = 'LD-DEGs_down-regulated'
payload = {
    'list': (None, genes_str),
    'description': (None, description)
}

response = requests.post(ENRICHR_URL, files=payload)
if not response.ok:
    raise Exception('Error analyzing gene list')

data = json.loads(response.text)
#print(data)

user_lid= data['userListId']

#########DisGeNET

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'
user_list_id = user_lid
gene_set_library = 'DisGeNET'
response = requests.get(
    ENRICHR_URL + query_string % (user_list_id, gene_set_library)
 )
if not response.ok:
    raise Exception('Error fetching enrichment results')

data = json.loads(response.text)
#print(data)

#dowload

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/export'
query_string = '?userListId=%s&filename=%s&backgroundType=%s'
user_list_id = user_lid
filename = 'DisGeNET_enrichment'
gene_set_library = 'DisGeNET'

url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
response = requests.get(url, stream=True)

with open('./data/LD_DEG_dn_'+filename + '.txt', 'wb') as f:
    for chunk in response.iter_content(chunk_size=1024): 
        if chunk:
            f.write(chunk)
