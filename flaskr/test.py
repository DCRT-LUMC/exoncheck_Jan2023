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
        print('NC_variant', NC_variant)
        hg38_coordinates = NC_variant.split('.')[-1].split('_')

        if len(hg38_coordinates) == 1:
            coordinate1 = ''.join([i for i in hg38_coordinates[0] if i.isdigit()])
            hg38_coordinates_for_gene_based_lovd = 'g.' + coordinate1
        else:
            coordinate1 = hg38_coordinates[0]
            coordinate2 = ''.join([i for i in hg38_coordinates[1] if i.isdigit()])
            hg38_coordinates_for_gene_based_lovd = 'g.' + coordinate1 + '_' + coordinate2
    except:
        hg38_coordinates_for_gene_based_lovd = 'N/A'

    return hg38_coordinates_for_gene_based_lovd


def exploit_variant_validator(MANE_select_NM_variant):
    """
    This function retrieves all VariantValidator information
    Input: The variant with the MANE select as reference
    Output: TO DO
    """

    NM_id = MANE_select_NM_variant.split(':')[0]

    req_variantvalidator = requests.get(f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/'
                                        f'{MANE_select_NM_variant}/{NM_id}?content-type=application%2Fjson')
    data_variantvalidator = json.loads(req_variantvalidator.content)

    req_gene2transcripts = requests.get(f'https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts_v2/'
                                        f'{NM_id}/{NM_id}?content-type=application%2Fjson')
    data_gene2transcripts = json.loads(req_gene2transcripts.content)

    # Get ENSG identifier
    try:
        ENSG_gene = data_variantvalidator[MANE_select_NM_variant]["gene_ids"]["ensembl_gene_id"]
    except:
        ENSG_gene = 'N/A'

    # Get gene symbol
    try:
        gene_symbol = data_variantvalidator[MANE_select_NM_variant]["gene_symbol"]
    except:
        gene_symbol = 'N/A'

    # Get hg38 variant
    try:
        NC_variant = data_variantvalidator[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"] \
            ["hgvs_genomic_description"]
    except:
        NC_variant = 'N/A'

    # Get hg38 variant position information
    try:
        hg38_chr = data_variantvalidator[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["vcf"]["chr"][3:]
        hg38_pos = data_variantvalidator[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["vcf"]["pos"]
        hg38_ref = data_variantvalidator[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["vcf"]["ref"]
        hg38_alt = data_variantvalidator[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["vcf"]["alt"]
        hg38_variant = hg38_chr + '-' + hg38_pos + '-' + hg38_ref + '-' + hg38_alt
    except:
        hg38_variant = 'N/A'

    # Get consequence of variant at protein level
    try:
        consequence_variant = data_variantvalidator[MANE_select_NM_variant]["hgvs_predicted_protein_consequence"]["tlr"]
        consequence_variant = consequence_variant.split('.')[-1][1:-1]  # Extract only mutation type part and
        # remove the brackets
    except:
        consequence_variant = 'N/A'

    # Get the latest NC reference sequence, which is later used for getting exon information
    try:
        latest_reference_sequence = data_variantvalidator[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["hgvs_genomic_description"].split(':')[0]
    except:
        latest_reference_sequence = 'N/A'

    # Get exon number
    try:
        start_exon_number = data_variantvalidator[MANE_select_NM_variant]["variant_exonic_positions"] \
            [latest_reference_sequence]["start_exon"]
        end_exon_number = data_variantvalidator[MANE_select_NM_variant]["variant_exonic_positions"] \
            [latest_reference_sequence]["end_exon"]

        total_exons = str(data_gene2transcripts["transcripts"][0]["genomic_spans"] \
                              [latest_reference_sequence]["total_exons"])

        # If variant covers multiple exons, save the first and last involved exons in variable exon_number
        if start_exon_number == end_exon_number:
            exon_number = start_exon_number
        else:
            exon_number = start_exon_number + 'till' + end_exon_number
    except:
        exon_number = 'N/A'
        total_exons = 'N/A'

    # Get total protein length
    try:
        coding_end = data_gene2transcripts["transcripts"][0]["coding_end"]
        coding_start = data_gene2transcripts["transcripts"][0]["coding_start"]
        total_protein_length = round((abs(coding_end - coding_start) + 1) / 3)  # Note we have to do +1 because of the
        # way of counting
    except:
        total_protein_length = 'N/A'

    # First and last exons can't be skipped. If the variant concerns the first or last exon, show message. Else leave
    # this message empty. !!!This needs to be designed better!!!
    exon_number_interpretation = ''

    # Get the exon skip in NC format and save the exon length
    # Besides, calculate distance to nearest splice site
    hg38_coordinates = reformat_hg38_positions(NC_variant)
    print('hg38_coordinates', hg38_coordinates)
    lower_limit_variant_hg38 = hg38_coordinates.split('.')[-1].split('_')[0]
    try:
        upper_limit_variant_hg38 = hg38_coordinates.split('.')[-1].split('_')[1]
    except:
        upper_limit_variant_hg38 = lower_limit_variant_hg38

    try:
        for exon in data_gene2transcripts["transcripts"][0]["genomic_spans"][latest_reference_sequence] \
                ["exon_structure"]:
            if str(exon["exon_number"]) == exon_number:
                print('exon: ', exon)
                exon_end = exon["genomic_end"]
                exon_start = exon["genomic_start"]
                exon_length = int(exon["cigar"][:-1])
                NC_exon_NC_format = latest_reference_sequence + ':g.' + str(exon_start) + '_' + str(exon_end) + 'del'

                # Include only the coding part for the exon length
                if exon_number == '1':
                    exon_length = exon_length - int(data_gene2transcripts["transcripts"][0]["coding_start"]) + 1
                    exon_number_interpretation = "First exon can't be skipped"
                elif exon_number == total_exons:
                    exon_length = exon_length - int(data_gene2transcripts["transcripts"][0]["coding_end"]) + 1
                    exon_number_interpretation = "Last exon can't be skipped"

                # Get nearest splice site
                print('abs(int(lower_limit_variant_hg38))', abs(int(lower_limit_variant_hg38)))
                print('abs(int(exon_start))', abs(int(exon_start)))
                distance_1 = abs(int(lower_limit_variant_hg38)) - abs(int(exon_start))
                distance_2 = abs(int(exon_end)) - abs(int(upper_limit_variant_hg38))
                if distance_1 < distance_2:
                    nearest_splice_distant = distance_1
                else:
                    nearest_splice_distant = distance_2
    except:
        nearest_splice_distant = 'N/A'

    print(nearest_splice_distant)

MANE_select_NM_variant = 'NM_000391.4:c.1397T>G'
exploit_variant_validator(MANE_select_NM_variant)


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
