import json

import requests

MANE_select_NM_variant = 'NM_001029896.2:c.749_751del'

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
    reference_sequences = data_variantvalidator[MANE_select_NM_variant]["variant_exonic_positions"].keys()
    latest_reference_sequence = ''
    for reference_sequence in reference_sequences:
        if reference_sequence.startswith('NC'):
            if reference_sequence > latest_reference_sequence:
                latest_reference_sequence = reference_sequence
except:
    latest_reference_sequence = 'N/A'
print(latest_reference_sequence)

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
try:
    for exon in data_gene2transcripts["transcripts"][0]["genomic_spans"][latest_reference_sequence] \
            ["exon_structure"]:
        if str(exon["exon_number"]) == exon_number:
            genomic_end = str(exon["genomic_end"])
            genomic_start = str(exon["genomic_start"])
            exon_length = int(exon["cigar"][:-1])
            NC_exon_NC_format = latest_reference_sequence + ':g.' + genomic_start + '_' + genomic_end + 'del'
            # Include only the coding part for the exon length
            if exon_number == '1':
                exon_length = exon_length - int(data_gene2transcripts["transcripts"][0]["coding_start"]) + 1
                exon_number_interpretation = "First exon can't be skipped"
            elif exon_number == total_exons:
                exon_length = exon_length - int(data_gene2transcripts["transcripts"][0]["coding_end"]) + 1
                exon_number_interpretation = "Last exon can't be skipped"
except:
    NC_exon_NC_format = 'N/A'
    exon_length = 'N/A'

# Convert exon length from nucleotides to amino acids
try:
    exon_length /= 3.0
except:
    exon_length = 'N/A'

# Check if exon is in frame or out-of-frame
try:
    if exon_length.is_integer():
        frame = 'In-frame'
    elif exon_length.isdecimal():
        frame = 'Out-of-frame'
except:
    frame = 'N/A'

# Get percentage of exon length compared to total protein length
try:
    percentage_length = round(exon_length / total_protein_length * 100, 2)
    exon_length = str(round(exon_length, 2))
except:
    percentage_length = 'N/A'
    exon_length = 'N/A'

# Get OMIM identifier
try:
    omim_id = data_variantvalidator[MANE_select_NM_variant]["gene_ids"]["omim_id"][0]  # Index of zero is necessary,
    # otherwise you get a list
except:
    omim_id = 'N/A'

# Below is about retrieving information about the exon skip
# This part needs to be revised and the updated VariantValidator API needs to be implemented (whenever
# VariantValidator is able to predict the consequence of skipping at protein level based on RNA-reference input
req_exon_variantvalidator = requests.get(
    f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/{NC_exon_NC_format}/ \ '
    f'mane_select?content-type=application%2Fjson')
data_exon_variantvalidator = json.loads(req_exon_variantvalidator.content)

# Get consequence of skipping
try:
    for key in data_exon_variantvalidator.keys():
        if key.startswith('NM'):
            MANE_select_NM_exon = key
    consequence_skipping = data_exon_variantvalidator[MANE_select_NM_exon] \
        ['hgvs_predicted_protein_consequence']['tlr']
except:
    MANE_select_NM_exon = 'N/A'
    consequence_skipping = 'N/A'
