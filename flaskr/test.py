import json
import pandas as pd
import numpy as np

import requests

def reformat_hg38_positions(NC_variant):
    # This function converts the hg38 genomic description (i.e. NC_000023.11:g.49075445_49075447del) to a format
    # that is accepted by LOVD (i.e. g.49075445_49075447). Note that chromosome information is not needed since
    # LOVD already known to which chromosome the gene belongs
    # Input: hg38 genomic description
    # Output: LOVD conform hg38 description
    try:
        hg38_coordinates = NC_variant.split('.')[-1].split('_')

        if len(hg38_coordinates) == 1:
            hg38_coordinates_for_gene_based_lovd = 'g.'.join([i for i in hg38_coordinates[0] if i.isdigit()])
        else:
            coordinate1 = hg38_coordinates[0]
            coordinate2 = ''.join([i for i in hg38_coordinates[1] if i.isdigit()])
            hg38_coordinates_for_gene_based_lovd = 'g.' + coordinate1 + '_' + coordinate2
    except:
        hg38_coordinates_for_gene_based_lovd = 'N/A'

    return hg38_coordinates_for_gene_based_lovd


def get_lovd_info(hg38_variant, NC_variant):
    # This function checks if there are exact matches available in the LOVD database (TO DO: add credits)
    # Currently it's only checking the hg19 based database, this should be elaborated with hg18 and hg17
    # Input: the variant in hg38 format, the variant with NC as reference
    # Output: A string containing the number of hits per found gene and the corresponding link
    # TO DO: IMPROVE THIS OUTPUT FORMAT, THINK OF A WAY TO CONVENIENTLY STORE THIS IN A SQL FORMAT

    exact_lovd_match_link = 'N/A'

    final_dict = dict()

    # First get the coordinates in the right format
    chromosome_of_variant = hg38_variant.split('-')[0]
    hg38_coordinates_for_gene_lovd = reformat_hg38_positions(NC_variant)
    hg38_coordinates_for_general_lovd = 'chr' + chromosome_of_variant + ':' + hg38_coordinates_for_gene_lovd[2:]

    # Find in which genes exact matches are found
    try:

        genes_containing_exact_hits = []
        dna_containing_exact_hits = []
        url = f'http://lovd.nl/search.php?build=hg38&position={hg38_coordinates_for_general_lovd}'
        url_data = pd.read_table(url, sep='\t')

        for gene in url_data['gene_id']:
            genes_containing_exact_hits.append(gene)

        for dna in url_data['DNA']:
            dna_containing_exact_hits.append(dna)

        for i in range(len(genes_containing_exact_hits)):
            final_dict[genes_containing_exact_hits[i]] = dna_containing_exact_hits[i]

        # Remove duplicates if any
        genes_containing_exact_hits = list(dict.fromkeys(genes_containing_exact_hits))
        dna_containing_exact_hits = list(dict.fromkeys(dna_containing_exact_hits))

    except:
        genes_containing_exact_hits = 'N/A'
        dna_containing_exact_hits = 'N/A'

    # Loop over the genes containing hits
    if genes_containing_exact_hits != 'N/A':
        output_exact_hits = ''

        for gene in genes_containing_exact_hits:
            # Start counting from zero again when querying a new gene
            number_exact_lovd_matches = 0

            # Check if variant position EXACTLY matches other variants
            try:
                req_exact = requests.get(
                    f'http://databases.lovd.nl/shared/api/rest.php/variants/{gene}?search_position={final_dict[gene]}&format=application/json')

                data_exact = json.loads(req_exact.content)

                for variant in data_exact:
                    number_exact_lovd_matches += 1
                    lovd_DBID = variant["Variant/DBID"]
                    exact_lovd_match_link = f'https://databases.lovd.nl/shared/view/{gene}' \
                                            f'?search_VariantOnGenome%2FDBID=%22{lovd_DBID}%22'
            except:
                exact_lovd_match_link = 'N/A'

            output_exact_hits += gene + ': ' + str(
                number_exact_lovd_matches) + ' hit(s), link: ' + exact_lovd_match_link + ', '

            # https://databases.lovd.nl/shared/view/{gene}?search_VariantOnGenome%2FDBID=%22BEST1_000022%22

    else:
        output_exact_hits = 'N/A'

    return output_exact_hits

get_lovd_info('X-49075439-AAGG-A', 'NC_000023.11:g.49075445_49075447del')


# file1 = open('data/normal_data_retina.csv', 'r')
# lines = file1.readlines()
#
# outputfile = open('data/selected_normal_data_retina.csv', 'w')
#
# for line in lines[1:]: # Skip header
#     # print(line.split(',')[5])
#     gene = line.split(',')[1]
#     retina_periphery_s1 = float(line.split(',')[-6])
#     retina_periphery_s2 = float(line.split(',')[-5])
#     retina_periphery_s3 = float(line.split(',')[-4])
#     retina_center_s1 = float(line.split(',')[-3])
#     retina_center_s2 = float(line.split(',')[-2])
#     retina_center_s3 = float(line.split(',')[-1])
#
#     # Check if all periphery are non zero
#     if retina_periphery_s1 != 0.0 and retina_periphery_s2 != 0.0:
#         retina_periphery = 'True'
#     elif retina_periphery_s1 != 0.0 and retina_periphery_s3 != 0.0:
#         retina_periphery = 'True'
#     elif retina_periphery_s2 != 0.0 and retina_periphery_s3 != 0.0:
#         retina_periphery = 'True'
#     else:
#         retina_periphery = 'False'
#
#     # Check if all periphery are non zero
#     if retina_center_s1 != 0.0 and retina_center_s2 != 0.0:
#         retina_center = 'True'
#     elif retina_center_s1 != 0.0 and retina_center_s3 != 0.0:
#         retina_center = 'True'
#     elif retina_center_s2 != 0.0 and retina_center_s3 != 0.0:
#         retina_center = 'True'
#     else:
#         retina_center = 'False'
#
#     outputfile.write(gene + ',' + retina_periphery +  ',' + retina_center + '\n')
