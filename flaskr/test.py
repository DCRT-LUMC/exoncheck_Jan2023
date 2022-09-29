import requests
import json
import xmltodict
import urllib.parse
import urllib.request
import pandas as pd



def get_uniprot_info(ENSG_gene_id):
    # This function provides the UniProt link (TO DO: add credits) and general domain information
    # (note: not exon to be skipped focused)
    # Input: ENSG gene identifier (without version number)
    # Output: UniProt gene link and exon corresponding domain information
    try:
        req = requests.get(f'https://mygene.info/v3/gene/{ENSG_gene_id}?fields=uniprot')
        data = json.loads(req.content)
        uniprot_id = data['uniprot']['Swiss-Prot']
        # Get UniProt gene link
        uniprot_link = f'https://www.uniprot.org/uniprotkb/{uniprot_id}/entry'
    except:
        uniprot_id = 'N/A'
        uniprot_link = 'N/A'

    # Get general domain info
    try:
        file = urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.xml')
        data = file.read()
        file.close()
        # Convert xml to dict
        uniprot_dict = xmltodict.parse(data)
        print(uniprot_dict)
        # Retrieve domain information
        domain_info = uniprot_dict['uniprot']['entry']['comment'][7]['text']['#text']
    except:
        domain_info = 'N/A'

    return uniprot_link, domain_info

get_uniprot_info("ENSG00000196998")