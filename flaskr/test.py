import json
import pandas as pd
import numpy as np

import requests

try:
    genes_containing_exact_hits = []
    url = 'http://lovd.nl/search.php?build=hg19&position=chr13:32936733'
    url_data = pd.read_table(url, sep = '\t')

    for gene in url_data['gene_id']:
        genes_containing_exact_hits.append(gene)

    # Remove duplicates if any
    genes_containing_exact_hits = list(dict.fromkeys(genes_containing_exact_hits))


    print(genes_containing_exact_hits)

except:
    data_exact = 'N/A'
