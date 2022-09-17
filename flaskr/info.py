import requests
import json
import xmltodict
import urllib.parse
import urllib.request

def check_for_hgvs_format(uploaded_variant):
    """
    This function checks the syntax of the uploaded variant by making use of checkHGVS API from LOVD (credits: TO DO)
    Input: uploaded variant (protein coding RNA), reference can be the MANE select as any other transcript
    Output: When the original syntax is wrong:
        This function returns a string containing all warnings and error messages, the suggested syntax correction if
        any with its corresponding confidence value.
    When the original syntax is correct:
        This function returns an empty string
    """

    syntax_message = ''

    # Get syntax data from checkHGVS
    try:
        req = requests.get(f'https://api.lovd.nl/v1/checkHGVS/{uploaded_variant}')
        data = json.loads(req.content)
        warnings = data['data'][uploaded_variant]["warnings"]
        errors = data['data'][uploaded_variant]["errors"]
        suggested_corrections = data['data'][uploaded_variant]["data"]["suggested_correction"]
    except:
        warnings = []
        errors = []
        suggested_corrections = []

    # Merge all warnings in one syntax message if any warning exists
    if warnings != []:
        syntax_message += f'Warning(s) = '
        for warning in warnings.values():
            syntax_message += warning + ' '

    # Add error messages to the syntax message if any error exists
    if errors != []:
        syntax_message += f'Error(s) = '
        for error in errors.values():
            syntax_message += error + ' '

    # Add suggested correction to the syntex message if any suggestion exists
    # A new try-except statement is used since suggested correction are not provided for all wrong syntaxes
    if suggested_corrections != []:
        try:
            syntax_message += f'Suggested correction = {suggested_corrections["value"]} (confidence: ' \
                              f'{suggested_corrections["confidence"]})'
        except:
            syntax_message += 'No suggested corrections available'

    return syntax_message

def check_for_match_variant_and_transcript(uploaded_variant):
    """
    This function checks if the uploaded variant does make sense on a biological level by employing VariantValidator
    (credits: TO DO)
    Input: uploaded variant (protein coding RNA), reference can be the MANE select as any other transcript
    Output: When the variant does not make sense on a biological level:
        This function returns a string with the error message
    When the variant does make sense on a biological level:
        This function returns an empty string
    """

    try:
        req = requests.get(f'https://reg.genome.network/allele?hgvs={uploaded_variant}')
        data = json.loads(req.content)
        match_message = data['message']
    except:
        match_message = ''

    return match_message

def get_MANE_select_identifiers(uploaded_variant):
    """
    This function selects the MANE select transcript of the uploaded variant by employing VariantValidator
    (credits: TO DO)
    Input: uploaded variant (protein coding RNA), reference can be the MANE select as any other transcript
    Output: the variant based on the MANE select transcript in NM- and ENST-based format
    """

    try:
        req = requests.get(f'https://reg.genome.network/allele?hgvs={uploaded_variant}')
        data = json.loads(req.content)
        for transcript in data['transcriptAlleles']:
            if 'MANE' in transcript:
                MANE_select_NM_variant = transcript['MANE']['nucleotide']['RefSeq']['hgvs']
                MANE_select_ENST_variant = transcript['MANE']['nucleotide']['Ensembl']['hgvs']
                break  # ClinGen provides the (same) MANE select twice, only one is needed
    except:
        MANE_select_NM_variant = 'N/A'
        MANE_select_ENST_variant = 'N/A'

    return MANE_select_NM_variant, MANE_select_ENST_variant

def get_strand(ENSG_gene_id):
    try:
        req = requests.get(f'https://rest.ensembl.org/lookup/id/{ENSG_gene_id}?content-type=application/json')
        data = json.loads(req.content)
        if data['strand'] == -1:
            return 'Reverse'
        elif data['strand'] == 1:
            return 'Forward'
    except:
        return 'N/A'

def exploit_variant_validator(MANE_select_NM_variant):
    NC_variant = 'N/A'
    hg38_variant = 'N/A'
    ENSG_gene = 'N/A'
    omim_id = 'N/A'
    gene_symbol = 'N/A'
    consequence_variant = 'N/A'
    exon_number = 'N/A'
    total_exons = 'N/A'
    NC_exon_NC_format = 'N/A'
    exon_length = 'N/A'
    total_protein_length = 'N/A'
    percentage_length = 'N/A'
    frame = 'N/A'
    consequence_skipping = 'N/A'
    MANE_select_NM_exon = 'N/A'

    NM_id = MANE_select_NM_variant.split(':')[0]

    req = requests.get(f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/{MANE_select_NM_variant}/{NM_id}?content-type=application%2Fjson')
    data = json.loads(req.content)

    req2 = requests.get(f'https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts_v2/{MANE_select_NM_variant}/{NM_id}?content-type=application%2Fjson')
    data2 = json.loads(req2.content)

    # Get hg38 positions
    try:
        NC_variant = data[MANE_select_NM_variant]["primary_assembly_loci"]["hg19"]["hgvs_genomic_description"]
    except:
        NC_variant = 'NA'

    # Get hg38_variant
    try:
        hg38_chr = data[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["vcf"]["chr"][3:]
        hg38_pos = data[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["vcf"]["pos"]
        hg38_ref = data[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["vcf"]["ref"]
        hg38_alt = data[MANE_select_NM_variant]["primary_assembly_loci"]["hg38"]["vcf"]["alt"]
        hg38_variant = hg38_chr + '-' + hg38_pos + '-' + hg38_ref + '-' + hg38_alt
    except:
        hg38_variant = 'N/A'

    # Get ENSG identifier
    try:
        ENSG_gene = data[MANE_select_NM_variant]["gene_ids"]["ensembl_gene_id"]
    except:
        ENSG_gene = 'N/A'

    # Get OMIM identifier
    try:
        omim_id = data[MANE_select_NM_variant]["gene_ids"]["omim_id"][0] # Index of zero is necessary, otherwise you get a list
    except:
        omim_id = 'N/A'

    # Get gene symbol
    try:
        gene_symbol = data[MANE_select_NM_variant]["gene_symbol"]
    except:
        gene_symbol = 'N/A'

    # Get mutation type
    try:
        consequence_variant = data[MANE_select_NM_variant]["hgvs_predicted_protein_consequence"]["tlr"]
        consequence_variant = consequence_variant.split('.')[-1][1:-1] # Extract only mutation type part and remove the brackets
    except:
        consequence_variant = 'N/A'

    # Get the latest reference sequence
    try:
        reference_sequences = data[MANE_select_NM_variant]["variant_exonic_positions"].keys()
        latest_reference_sequence = ''
        for reference_sequence in reference_sequences:
            if reference_sequence.startswith('NC'):
                if reference_sequence > latest_reference_sequence:
                    latest_reference_sequence = reference_sequence
    except:
        latest_reference_sequence ='N/A'


    # Get exon number
    try:
        start_exon_number = data[MANE_select_NM_variant]["variant_exonic_positions"][latest_reference_sequence]["start_exon"]
        end_exon_number = data[MANE_select_NM_variant]["variant_exonic_positions"][latest_reference_sequence]["end_exon"]

        total_exons = str(data2["transcripts"][0]["genomic_spans"][latest_reference_sequence]["total_exons"])

        if start_exon_number == end_exon_number:
            exon_number = start_exon_number
        else:
            exon_number = start_exon_number + '+' + end_exon_number
    except:
        exon_number = 'N/A'
        total_exons = 'N/A'

    # Exon to be skipped and get total protein length
    NC_exon_NC_format = 'unknown'
    exon_length = 0.0

    try:
        coding_end = data2["transcripts"][0]["coding_end"]
        coding_start = data2["transcripts"][0]["coding_start"]
        total_protein_length = round((abs(coding_end - coding_start) + 1) / 3)
    except:
        total_protein_length = 'N/A'

    exon_number_interpretation = ''

    try:
        for exon in data2["transcripts"][0]["genomic_spans"][latest_reference_sequence]["exon_structure"]:
            if str(exon["exon_number"]) == exon_number:
                genomic_end = str(exon["genomic_end"])
                genomic_start = str(exon["genomic_start"])
                exon_length = int(exon["cigar"][:-1])
                NC_exon_NC_format = latest_reference_sequence + ':g.' + genomic_start + '_' + genomic_end + 'del'
                if exon_number == '1' or exon_number == total_exons:
                    exon_length = exon_length - int(data2["transcripts"][0]["coding_start"]) + 1
                elif exon_number == total_exons:
                    exon_length = exon_length - int(data2["transcripts"][0]["coding_end"]) + 1

                if exon_number == '1':
                    exon_number_interpretation = "First exon can't be skipped"
                if exon_number == total_exons:
                    exon_number_interpretation = "Last exon can't be skipped"
    except:
        exon_length = 'N/A'
        NC_exon_NC_format = 'N/A'

    # looking at amino acids instead of nucleotides
    try:
        exon_length /= 3.0
    except:
        exon_length = 0

    # Check if exon is in frame or out-of-frame
    try:
        if exon_length.is_integer():
            frame = 'In-frame'
        else:
            frame = 'Out-of-frame'
    except:
        frame = 'N/A'

    # check percentage of protein length
    try:
        percentage_length = round(exon_length / total_protein_length * 100, 2)

        exon_length = str(round(exon_length, 2))
    except:
        percentage_length = 'N/A'
        exon_length = 'N/A'

    # Exon to be skipped
    # req_NC_exon = requests.get(f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/{NC_exon_NC_format}/mane_select?content-type=application%2Fjson')
    # data_NC_exon = json.loads(req_NC_exon.content)


    # consequence of skipping
    req3 = requests.get(
        f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/{NC_exon_NC_format}/mane_select?content-type=application%2Fjson')
    data3 = json.loads(req3.content)

    try:
        for key in data3.keys():
            if key.startswith('NM'):
                MANE_select_NM_exon = key

        consequence_skipping = data3[MANE_select_NM_exon]['hgvs_predicted_protein_consequence']['tlr']
    except:
        consequence_skipping = 'N/A'

    return \
        NC_variant, \
        hg38_variant, \
        ENSG_gene, \
        omim_id, \
        gene_symbol, \
        consequence_variant, exon_number, total_exons, exon_number_interpretation, NC_exon_NC_format, \
        exon_length, total_protein_length, percentage_length, frame, consequence_skipping, MANE_select_NM_exon

def get_positions_for_lovd(hg38_genomic_description):
    try:
        hg38_coordinates = hg38_genomic_description.split('.')[-1].split('_')

        if len(hg38_coordinates) == 1:
            hg38_coordinates_for_lovd = ''.join([i for i in hg38_coordinates[0] if i.isdigit()])
        else:
            coordinate1 = hg38_coordinates[0]
            coordinate2 = ''.join([i for i in hg38_coordinates[1] if i.isdigit()])
            hg38_coordinates_for_lovd = coordinate1 + '_' + coordinate2
    except:
        hg38_coordinates_for_lovd = 'N/A'

    return hg38_coordinates_for_lovd

def get_lovd_info(hg38_genomic_description, gene_symbol):
    hg38_coordinates_for_lovd = get_positions_for_lovd(hg38_genomic_description)
    no_partial_lovd_matches = 0
    no_exact_lovd_matches = 0

    # Check if gene symbol is available in the LOVD database
    try:
        req_gene = requests.get(
            f'https://databases.lovd.nl/shared/api/rest.php/genes/{gene_symbol}')
        lovd_gene_link = f'https://databases.lovd.nl/shared/variants/{gene_symbol}/unique'
    except:
        lovd_gene_link = 'Gene is not available in LOVD'

    # Check if variant position EXACTLY matches other variants
    req_exact = requests.get(
        f'https://databases.lovd.nl/shared/api/rest.php/variants/{gene_symbol}?search_position=g.{hg38_coordinates_for_lovd}&format=application/json')

    try:
        data_exact = json.loads(req_exact.content)

        exact_matching_lovd_variants = []
        exact_lovd_match_link = 'N/A'



        if data_exact:
            for variant in data_exact:
                exact_matching_lovd_variants.append(variant['id'])
                no_exact_lovd_matches += 1
                lovd_DBID = variant["Variant/DBID"]
                exact_lovd_match_link = f'https://databases.lovd.nl/shared/view/{gene_symbol}?search_VariantOnGenome%2FDBID=%22{lovd_DBID}%22' # Maybe move to outside of the loop
        else:
            exact_matching_lovd_variants.append('N/A')
    except:
        exact_lovd_match_link = 'N/A'
        exact_lovd_match_link =  'N/A'
        no_exact_lovd_matches = 0


    # Check if variant position PARTIALLY matches other variants
    partial_matching_lovd_variants = ['', '', '', '', '', '', '', '', '', '']

    try:
        req = requests.get(
            f'https://databases.lovd.nl/shared/api/rest.php/variants/{gene_symbol}?search_position=g.{hg38_coordinates_for_lovd}&position_match=partial&format=application/json')

        data = json.loads(req.content)
    except:
        data = False


    if data:
        for variant in data:
            # 'https://databases.lovd.nl/shared/variants/' + lovd_id
            try:
                partial_matching_lovd_variants[no_partial_lovd_matches] = variant['id']
                no_partial_lovd_matches += 1
            except:
                continue
    if not partial_matching_lovd_variants:
        partial_matching_lovd_variants.append('N/A')



    partial_lovd_match1, partial_lovd_match2, partial_lovd_match3, partial_lovd_match4, \
    partial_lovd_match5, partial_lovd_match6, partial_lovd_match7, partial_lovd_match8, \
    partial_lovd_match9, partial_lovd_match10 =  partial_matching_lovd_variants



    return lovd_gene_link, no_exact_lovd_matches, exact_lovd_match_link, no_partial_lovd_matches, partial_lovd_match1, partial_lovd_match2, partial_lovd_match3, partial_lovd_match4, \
    partial_lovd_match5, partial_lovd_match6, partial_lovd_match7, partial_lovd_match8, \
    partial_lovd_match9, partial_lovd_match10

def get_uniprot_info(ENSG_gene_id):
    # get uniprot id
    try:
        req = requests.get(f'https://mygene.info/v3/gene/{ENSG_gene_id}?fields=uniprot')
        data = json.loads(req.content)
        uniprot_id = data['uniprot']['Swiss-Prot']
        # get uniprot gene link
        uniprot_link = f'https://www.uniprot.org/uniprotkb/{uniprot_id}/entry'
    except:
        uniprot_id = 'N/A'
        uniprot_link = 'N/A'

    # get domain info
    try:
        file = urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.xml')
        data = file.read()
        file.close()
        uniprot_dict = xmltodict.parse(data)

        domain_info = uniprot_dict['uniprot']['entry']['comment'][7]['text']['#text']
    except:
        domain_info = 'N/A'

    return uniprot_link, domain_info

def get_gtexportal_json(ENSG_gene_id):
    url_gtexportal = f'https://gtexportal.org/rest/v1/expression/medianTranscriptExpression?datasetId=gtex_v8&gencodeId={ENSG_gene_id}&format=json'
    r_gtex = requests.get(url_gtexportal)
    gtex_data = json.loads(r_gtex.text)
    return gtex_data

def get_gene_expression(ENSG_gene_id, MANE_select_ENST_variant):
    # This function checks in which tissues the gene is expressed
    # It uses the GTExPortal 'Bulk tissue gene expression' data
    # (i.e. from https://gtexportal.org/rest/v1/expression/geneExpression?datasetId=gtex_v8&gencodeId=ENSG00000196998.17&format=json)
    # Need to be improved, wait for Marlen

    ENSG_version = 0
    ENST_without_version = MANE_select_ENST_variant.split('.')[0] + '.'

    while get_gtexportal_json(ENSG_gene_id + '.' + str(ENSG_version))[
        'medianTranscriptExpression'] == []:
        ENSG_version += 1

    # IMPLEMENT THIS LATER IN THE TOOL
    # if ENSG_gene_id != ENSG_gene_id + '.' + str(ENSG_version):
    #     print('\n***!WARNING!***\n' + ENSG_gene_id + '.' + str(
    #         ENSG_version) + ' is utilized instead of ' + ENSG_gene_id +
    #           ' to consult GTEx portal\n')

    gtex_data = get_gtexportal_json(ENSG_gene_id + '.' + str(ENSG_version))

    expression_eye = 'https://www.eye-transcriptome.com/search_latest.php'
    expression_brain = 'no'
    expression_fibroblasts = 'no'
    expression_tibial_nerve = 'no'
    expression_blood = 'no'
    expression_transformed_lymphocytes = 'no'

    for transcript_in_tissue in gtex_data['medianTranscriptExpression']:
        try:
            if ENST_without_version in transcript_in_tissue[
                'transcriptId'] and 'Brain' in transcript_in_tissue[
                'tissueSiteDetailId'] and transcript_in_tissue['median'] != 0:
                expression_brain = 'yes'

            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                transcript_in_tissue[
                    'tissueSiteDetailId'] == 'Cells_Cultured_fibroblasts' and \
                transcript_in_tissue['median'] != 0:
                expression_fibroblasts = 'yes'

            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                transcript_in_tissue[
                    'tissueSiteDetailId'] == 'Nerve_Tibial' and transcript_in_tissue[
                'median'] != 0:
                expression_tibial_nerve = 'yes'

            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                transcript_in_tissue[
                    'tissueSiteDetailId'] == 'Whole_Blood' and transcript_in_tissue[
                'median'] != 0:
                expression_blood = 'yes'

            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                transcript_in_tissue[
                    'tissueSiteDetailId'] == 'Cells_EBV-transformed_lymphocytes' and \
                transcript_in_tissue['median'] != 0:
                expression_transformed_lymphocytes = 'yes'
        except:
            expression_brain = 'N/A'
            expression_fibroblasts = 'N/A'
            expression_tibial_nerve = 'N/A'
            expression_blood = 'N/A'
            expression_transformed_lymphocytes = 'N/A'

    return expression_eye, expression_brain, expression_fibroblasts, expression_tibial_nerve, expression_blood, expression_transformed_lymphocytes
