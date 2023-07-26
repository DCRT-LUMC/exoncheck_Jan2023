import requests
import json
import xmltodict
import urllib.parse
import urllib.request
import pandas as pd
import re

def gene_to_exons(input_gene):
    req = requests.get(f'https://reg.test.genome.network/gene?HGNC.symbol={input_gene}')
    data = json.loads(req.content)

    NM_id = (data['externalRecords']['MANE'][0]['nucleotide']['RefSeq']['id'])

    data_gene2transcripts = fetch_gene2transcript(NM_id)

    NC_id = list(data_gene2transcripts["transcripts"][0]["genomic_spans"].keys())[1]

#    in_frame = {}
#    out_of_frame = {}

    in_frame_exons = []
    out_of_frame_exons = []
    for exon in data_gene2transcripts["transcripts"][0]["genomic_spans"][NC_id]["exon_structure"]:
    #    coding_start = data_gene2transcripts["transcripts"][0]["coding_start"]
        exon_length = int(exon["cigar"][:-1])
        if (exon_length/3.0).is_integer():
#            in_frame[exon["exon_number"]] = [(transcript_start), (transcript_end)]
            in_frame_exons.append(exon["exon_number"])
        else:
#            out_of_frame[exon["exon_number"]] = [(transcript_start), (transcript_end)]
            out_of_frame_exons.append(exon["exon_number"])
    return NM_id, NC_id, in_frame_exons, out_of_frame_exons

def exon_skip(NM_id, latest_reference_sequence, exon_number):
    data_gene2transcripts = fetch_gene2transcript(NM_id)

    # Get total number of exons
    total_exons = str(data_gene2transcripts["transcripts"][0]["genomic_spans"] \
                  [latest_reference_sequence]["total_exons"])

    # Get total protein length
    try:
        coding_end = data_gene2transcripts["transcripts"][0]["coding_end"]
        coding_start = data_gene2transcripts["transcripts"][0]["coding_start"]
        # - 3 (stop codon does not produce an amino acid) + 1 (to account for mathematical subtraction)
        total_protein_length = (abs(coding_end - coding_start) - 2) / 3
    except:
        total_protein_length = 'N/A'

    # Get the coding exon positions
    try:
        for exon in data_gene2transcripts["transcripts"][0]["genomic_spans"][latest_reference_sequence]["exon_structure"]:
            transcript_start = exon["transcript_start"]
            transcript_end = exon["transcript_end"]
            if coding_start >= transcript_start and coding_start <= transcript_end:
                first_coding = exon["exon_number"]
            if coding_end >= transcript_start and coding_end <= transcript_end:
                last_coding = exon["exon_number"]

        coding_exons = str(first_coding) + " to " + str(last_coding)
    except:
        coding_exons = 'N/A'

    # Format the NC_exon_description
    try:
        for exon in data_gene2transcripts["transcripts"][0]["genomic_spans"][latest_reference_sequence] \
                ["exon_structure"]:
            if str(exon["exon_number"]) == exon_number:
                exon_end = exon["genomic_end"]
                exon_start = exon["genomic_start"]
                exon_length_nu = int(exon["cigar"][:-1])
                NC_exon_NC_format = latest_reference_sequence + ':g.' + str(exon_start) + '_' + str(exon_end) + 'del'
                # Include only the coding part for the exon length
                if exon_number == str(first_coding):
                    exon_length_nu = exon_length_nu - coding_start + 1
                    exon_number_interpretation = "First exon can't be skipped."
                elif exon_number == str(last_coding):
                    exon_length_nu = exon_length_nu - coding_end + 1
                    exon_number_interpretation = "Last exon can't be skipped."
                else:
                    exon_number_interpretation = ""

                # Get exon length and percentage in terms of nucleotides
                coding_exon_length = coding_end - coding_start
                percentage_length_nu = percentage_length = round(exon_length_nu / coding_exon_length * 100, 2)
    except:
        NC_exon_NC_format = 'N/A'
        exon_length_nu = 'N/A'
        coding_exon_length = 'N/A'
        percentage_length_nu = 'N/A'

    # Convert exon length from nucleotides to amino acids
    try:
        exon_length = exon_length_nu/3.0
    except:
        exon_length = 'N/A'

    # Check if exon is in frame or out-of-frame
    try:
        if exon_length.is_integer():
            frame = 'In-frame'
        elif exon_length != 'N/A':
            frame = 'Out-of-frame'
    except:
        frame = 'N/A'

    # Get interpretation of distance to nearest splice boundary
    try:
#        if frame == 'Out-of-frame':
        if distance_1 < distance_2: # closer to 5' end
            if distance_1/exon_length_nu <= 0.3:
                splice_dist_interpretation = 'Variant might be too close to the donor site.'
#            else:
#                splice_dist_interpretation = 'Variant in the last 30% of the coding exon. Possibly eligible.'
        elif distance_2 <= distance_1:
            if distance_2/exon_length_nu <= 0.3:
                splice_dist_interpretation = 'Variant might be too close to the acceptor site.'
#            else:
#                splice_dist_interpretation = 'Variant in the first 70% of the coding exon. Not eligible.'
    except:
        splice_dist_interpretation = ''

    # Get percentage of exon length compared to total protein length
    try:
        percentage_length = round(exon_length / total_protein_length * 100, 2)
        exon_length = str(round(exon_length, 2))
    except:
        percentage_length = 'N/A'
        exon_length = 'N/A'

    # Check exon length conditions
    if percentage_length <= 10:
        length_condition = 'Exon length might be small enough to be skipped.'
    elif percentage_length > 10 and percentage_length < 30:
        length_condition = 'Subject to your interpretation.'
    else:
        length_condition = 'Exon length might be too large to be skipped.'

    req_exon_variantvalidator = requests.get(
        f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/{NC_exon_NC_format}/ \ '
        f'mane_select?content-type=application%2Fjson')
    data_exon_variantvalidator = json.loads(req_exon_variantvalidator.content)

    # Get consequence of skipping in protein (issue solved by VV)
    try:
        for key in data_exon_variantvalidator.keys():
            if key.startswith(NM_id):
                MANE_select_NM_exon = key
        try:
            if '+' in MANE_select_NM_exon:
                # Format the r. exon skip id (need to fix this as per Ivo's instructions)
                MANE_select_exon_split = re.split('[_:.+]', MANE_select_NM_exon)
                NM_r_exon = (NM_id + ":r." + str(int(MANE_select_exon_split[4]) - int(MANE_select_exon_split[-1][-4])) + "_" + MANE_select_exon_split[5] + "del")
            else:
                NM_r_exon = MANE_select_NM_exon.replace("c", "r")

            # Query the API to get the consequence of skipping
            req_rnavariant = requests.get(f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/'
                                        f'{NM_r_exon}/{NM_id}?content-type=application%2Fjson') 
            data_rnavariant = json.loads(req_rnavariant.content)

            consequence_skipping = data_rnavariant[MANE_select_NM_exon]["rna_variant_descriptions"]["translation"]
            r_exon_skip = data_rnavariant[MANE_select_NM_exon]["rna_variant_descriptions"]["rna_variant"]
        except:
            consequence_skipping = 'N/A'
            r_exon_skip = 'N/A'
    except:
        MANE_select_NM_exon = 'N/A'

    return total_exons,\
    total_protein_length,\
    coding_exons,\
    NC_exon_NC_format,\
    coding_exon_length,\
    exon_number_interpretation,\
    exon_length,\
    frame,\
    percentage_length,\
    length_condition,\
    MANE_select_NM_exon,\
    consequence_skipping,\
    r_exon_skip


def get_NM_for_NP(uploaded_variant):
    req = requests.get(f'https://reg.genome.network/allele?hgvs={uploaded_variant}')
    data = json.loads(req.content)
    message = ""
    try:
        NM_variants = data["aminoAcidAlleles"][0]['hgvsMatchingTranscriptVariant']
#        NM_variants = json.dumps(NM_variants)
    except:
        message = data["message"]
        NM_variants = ""

    return message, NM_variants

def check_for_hgvs_format(uploaded_variant):
    """
    This function checks whether the syntax of the uploaded variant conforms to HGVS nomenclature by making use of checkHGVS API from LOVD (credits: TO DO)
    Input: uploaded variant (NM and NC), reference can be the MANE select or any other transcript
    Output: When the original syntax is wrong:
        This function returns a string containing all warnings and error messages, the suggested syntax correction if
        any with its corresponding confidence value.
    When the original syntax is correct:
        This function returns an empty string

    Works for both hg19 and hg38.
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

    # Add suggested correction to the syntax message if any suggestion exists
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
    Input: uploaded variant (genomic and transcript alleles in HGVS format), reference can be the MANE select or any other transcript
    Output: When the variant does not make sense on a biological level:
        This function returns a string with the error message
    When the variant does make sense on a biological level:
        This function returns an empty string

    Works for both hg19 and hg38.
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
    Input: uploaded variant (genomic and transcript alleles in HGVS format), reference can be the MANE select as any other transcript
    Output: the variant based on the MANE select transcript in NM- and ENST-based format

    Works for both hg19 and hg38.
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

# Can cache this (?)
def get_strand(ENSG_gene_id):
    """
    This function retrieves the strand of the gene from Ensembl
    Input: ENSG gene identifier
    Output: Forward/Reverse/N/A

    Works for both hg19 and hg38.
    """

    try:
        req = requests.get(f'https://rest.ensembl.org/lookup/id/{ENSG_gene_id}?content-type=application/json')
        data = json.loads(req.content)
        if data['strand'] == -1:
            return 'Reverse'
        elif data['strand'] == 1:
            return 'Forward'
    except:
        return 'N/A'

# hg38
def reformat_hg38_positions(NC_variant):
    """
    This function converts the hg38 genomic description (i.e. NC_000023.11:g.49075445_49075447del) to a format
    that is accepted by LOVD (i.e. g.49075445_49075447). Note that chromosome information is not needed since
    LOVD already known to which chromosome the gene belongs
    Input: hg38 genomic description
    Output: LOVD conform hg38 description
    """
    try:
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


def fetch_variantvalidator(transcript):
    # VariantValidation API endpoint helper strings
    vv = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"
    content_type = "?content-type=application%2Fjson"

    nm_id = transcript.split(':')[0]

    url = f"{vv}/hg38/{transcript}/{nm_id}{content_type}"
    req_variantvalidator = requests.get(url)

    data = json.loads(req_variantvalidator.content)

    with open('variantvalidator.json', 'wt') as fout:
        print(json.dumps(data, indent=True), file=fout)

    return data


def fetch_gene2transcript(transcript):
    vv = "https://rest.variantvalidator.org/VariantValidator/tools"
    content_type = "?content-type=application%2Fjson"

    nm_id = transcript.split(':')[0]

    #url = f"{vv}/gene2transcripts_v2/{nm_id}/{nm_id}/{content_type}"
    url = f"{vv}/gene2transcripts/{nm_id}{content_type}"
    req_gene2transcripts = requests.get(url)

    data = json.loads(req_gene2transcripts.content)

    # Throw out all non-MANE transcripts, IF there is at least one MANE transcript
    transcripts = data.get("transcripts", list())
    transcripts = [ts for ts in transcripts if ts["annotations"]["mane_select"]]
    if transcripts:
        data["transcripts"] = transcripts

    with open('gene2transcripts.json', 'wt') as fout:
        print(json.dumps(data, indent=True), file=fout)

    return data

def exploit_variant_validator(variant, build):
    """
    This function retrieves all VariantValidator information
    Input: The variant with the MANE select as reference
    Output: TO DO
    """

    print(f"Genome build: {build}")
    data_variantvalidator = fetch_variantvalidator(variant)

    data_gene2transcripts = fetch_gene2transcript(variant)


    # Get ENSG identifier
    try:
        ENSG_gene = data_variantvalidator[variant]["gene_ids"]["ensembl_gene_id"]
    except:
        ENSG_gene = 'N/A'

    # Get gene symbol
    try:
        gene_symbol = data_variantvalidator[variant]["gene_symbol"]
    except:
        gene_symbol = 'N/A'

    # Get variant
    try:
        NC_variant = data_variantvalidator[variant]["primary_assembly_loci"][build] \
            ["hgvs_genomic_description"]
        latest_reference_sequence = NC_variant.split(':')[0]
    except:
        NC_variant = 'N/A'
        latest_reference_sequence = 'N/A'

    # Get variant position information
    try:
        chr = data_variantvalidator[variant]["primary_assembly_loci"][build]["vcf"]["chr"][3:]
        pos = data_variantvalidator[variant]["primary_assembly_loci"][build]["vcf"]["pos"]
        ref = data_variantvalidator[variant]["primary_assembly_loci"][build]["vcf"]["ref"]
        alt = data_variantvalidator[variant]["primary_assembly_loci"][build]["vcf"]["alt"]
        hg_variant = chr + '-' + pos + '-' + ref + '-' + alt
    except:
        hg_variant = 'N/A'

    # Get consequence of variant at protein level
    try:
        consequence_variant = data_variantvalidator[variant]["hgvs_predicted_protein_consequence"]["tlr"]
        consequence_variant = consequence_variant.split('.')[-1][1:-1]  # Extract only mutation type part and
        # remove the brackets
    except:
        consequence_variant = 'N/A'

    # Get the latest NC reference sequence, which is later used for getting exon information
    try:
        latest_reference_sequence = data_variantvalidator[variant]["primary_assembly_loci"][build]["hgvs_genomic_description"].split(':')[0]
    except:
        latest_reference_sequence = 'N/A'

    # Get exon number
    try:
        start_exon_number = data_variantvalidator[variant]["variant_exonic_positions"] \
            [latest_reference_sequence]["start_exon"]
        end_exon_number = data_variantvalidator[variant]["variant_exonic_positions"] \
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

    print(f"exons: {exon_number}/{total_exons}")
    # Get total protein length
    try:
        coding_end = data_gene2transcripts["transcripts"][0]["coding_end"]
        coding_start = data_gene2transcripts["transcripts"][0]["coding_start"]
        # - 3 (stop codon does not produce an amino acid) + 1 (to account for mathematical subtraction)
        total_protein_length = (abs(coding_end - coding_start) - 2) / 3
    except:
        total_protein_length = 'N/A'

    print(f"protein length: {total_protein_length}")
    # Get the exon skip in NC format and save the exon length
    # Besides, calculate distance to nearest splice site
    if build == "hg19":
        coordinates = reformat_hg19_positions(NC_variant)
    elif build == "hg38":
        coordinates = reformat_hg38_positions(NC_variant)
    else:
        raise RuntimeError(f"Unknown reference build {build}")

    lower_limit_variant = coordinates.split('.')[-1].split('_')[0]
    try:
        upper_limit_variant = coordinates.split('.')[-1].split('_')[1]
    except:
        upper_limit_variant = lower_limit_variant

    print(f"Upper limit: {upper_limit_variant}")
    # Get the coding exon positions
    try:
        for exon in data_gene2transcripts["transcripts"][0]["genomic_spans"][latest_reference_sequence]["exon_structure"]:
            transcript_start = exon["transcript_start"]
            transcript_end = exon["transcript_end"]
            if coding_start >= transcript_start and coding_start <= transcript_end:
                first_coding = exon["exon_number"]
            if coding_end >= transcript_start and coding_end <= transcript_end:
                last_coding = exon["exon_number"]

        coding_exons = str(first_coding) + " to " + str(last_coding)
    except:
        coding_exons = 'N/A'

    print(f"coding_exons: {coding_exons}")
    # Format the NC_exon_description
    try:
        for exon in data_gene2transcripts["transcripts"][0]["genomic_spans"][latest_reference_sequence] \
                ["exon_structure"]:
            if str(exon["exon_number"]) == exon_number:
                exon_end = exon["genomic_end"]
                exon_start = exon["genomic_start"]
                exon_length_nu = int(exon["cigar"][:-1])
                NC_exon_NC_format = latest_reference_sequence + ':g.' + str(exon_start) + '_' + str(exon_end) + 'del'

                # Include only the coding part for the exon length
                coding_start = int(data_gene2transcripts["transcripts"][0]["coding_start"])
                coding_end = int(data_gene2transcripts["transcripts"][0]["coding_end"])
                if exon_number == str(first_coding):
                    exon_length_nu = exon_length_nu - coding_start + 1
                    exon_number_interpretation = "First exon can't be skipped."
                elif exon_number == str(last_coding):
                    exon_length_nu = exon_length_nu - coding_end + 1
                    exon_number_interpretation = "Last exon can't be skipped."
                else:
                    exon_number_interpretation = ""
                # Get exon length and percentage in terms of nucleotides
                # Get exon length and percentage in terms of nucleotides
                coding_exon_length = coding_end - coding_start
                percentage_length_nu = percentage_length = round(exon_length_nu / coding_exon_length * 100, 2)

                # Get nearest splice site and which side (5'/3')
                distance_1 = int(lower_limit_variant) - int(exon_start)
                distance_2 = int(exon_end) - int(upper_limit_variant)
                if distance_1 < distance_2:
                    nearest_splice_distant = distance_1
                    nearest_end = "5'"
                else:
                    nearest_splice_distant = distance_2
                    nearest_end = "3'"
    except:
        NC_exon_NC_format = 'N/A'
        exon_length_nu = 'N/A'
        nearest_splice_distant = 'N/A'
        nearest_end = 'N/A'
        coding_exon_length = 'N/A'
        percentage_length_nu = 'N/A'

    print(f"NC format: {NC_exon_NC_format}")
    # Convert exon length from nucleotides to amino acids
    try:
        exon_length = exon_length_nu/3.0
    except:
        exon_length = 'N/A'
    print("exon_length: " + str(exon_length))

    # Check if exon is in frame or out-of-frame
    try:
        if exon_length.is_integer():
            frame = 'In-frame'
        elif exon_length != 'N/A':
            frame = 'Out-of-frame'
    except:
        frame = 'N/A'

    print(f"frame: {frame}")
    # Get interpretation of distance to nearest splice boundary
    splice_dist_interpretation = ''
    try:
#        if frame == 'Out-of-frame':
        if distance_1 < distance_2: # closer to 5' end
            if distance_1 <= 15:
                splice_dist_interpretation = 'Variant might be too close to the donor site.'
            else:
                splice_dist_interpretation = ''
        else:
            if distance_2 <= 15:
                splice_dist_interpretation = 'Variant might be too close to the acceptor site.'
            else:
                splice_dist_interpretation = ''
    except:
        splice_dist_interpretation = ''

    print(f"splice_dist: {splice_dist_interpretation}")
    # Get percentage of exon length compared to total protein length
    try:
        percentage_length = round(exon_length / total_protein_length * 100, 2)
        exon_length = str(round(exon_length, 2))
    except:
        percentage_length = 'N/A'
        exon_length = 'N/A'

    # Check exon length conditions
    if percentage_length <= 10:
        length_condition = 'Exon length might be small enough to be skipped.'
    elif percentage_length > 10 and percentage_length < 30:
        length_condition = 'Subject to your interpretation.'
    else:
        length_condition = 'Exon length might be too large to be skipped.'

    # Get OMIM identifier
    try:
        omim_id = data_variantvalidator[variant]["gene_ids"]["omim_id"][0]  # Index of zero is necessary,
        # otherwise you get a list
    except:
        omim_id = 'N/A'

    print(f"OMIM: {omim_id}")
    # Below is about retrieving information about the exon skip
    # This part needs to be revised and the updated VariantValidator API needs to be implemented (whenever
    # VariantValidator is able to predict the consequence of skipping at protein level based on RNA-reference input)
    url = (
        f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/{build}/{NC_exon_NC_format}/'
        f'mane_select?content-type=application%2Fjson'
    )
    req_exon_variantvalidator = requests.get(url)
    data_exon_variantvalidator = json.loads(req_exon_variantvalidator.content)

    print(NC_exon_NC_format)
    print(url)
    NM_id = variant.split(':')[0]

    # Get consequence of skipping in protein (issue solved by VV)
    try:
        for key in data_exon_variantvalidator.keys():
            if key.startswith(NM_id):
                MANE_select_NM_exon = key
        try:
            if '+' in MANE_select_NM_exon:
                # Format the r. exon skip id (need to fix this as per Ivo's instructions)
                MANE_select_exon_split = re.split('[_:.+]', MANE_select_NM_exon)
                NM_r_exon = (NM_id + ":r." + str(int(MANE_select_exon_split[4]) - int(MANE_select_exon_split[-1][-4])) + "_" + MANE_select_exon_split[5] + "del")
            else:
                NM_r_exon = MANE_select_NM_exon.replace("c", "r")

            # Query the API to get the consequence of skipping
            req_rnavariant = requests.get(f'https://rest.variantvalidator.org/VariantValidator/variantvalidator/{build}/'
                                        f'{NM_r_exon}/{NM_id}?content-type=application%2Fjson')
            data_rnavariant = json.loads(req_rnavariant.content)

            consequence_skipping = data_rnavariant[MANE_select_NM_exon]["rna_variant_descriptions"]["translation"]
            r_exon_skip = data_rnavariant[MANE_select_NM_exon]["rna_variant_descriptions"]["rna_variant"]
        except:
            consequence_skipping = 'N/A'
            r_exon_skip = 'N/A'
            MANE_select_NM_exon = 'N/A'
        print(f"consequence: {consequence_skipping}, r_exon_skip: {r_exon_skip}")
    except:
        MANE_select_NM_exon = 'N/A'

    exon_number_interpretation = 'N/A'
    length_condition = 'N/A'
    return \
        NC_variant, \
        hg_variant, \
        ENSG_gene, \
        omim_id, \
        gene_symbol, \
        consequence_variant, \
        exon_number, \
        total_exons, \
        exon_number_interpretation, \
        coding_exons, \
        NC_exon_NC_format, \
        exon_length, \
        nearest_splice_distant, \
        nearest_end, \
        total_protein_length, \
        percentage_length, \
        length_condition, \
        frame, \
        splice_dist_interpretation, \
        consequence_skipping, \
        r_exon_skip, \
        MANE_select_NM_exon

def check_gene_eligibility(gene_symbol):
    """
    This function checks whether the gene is contained in the DCRT_gene_list and prints out its eligibility.
    Input: HGNC Gene Symbol
    Output: Eligibility 1/2/3/N/A as per gene list

    Gene list stored in data file, need to reupload in case of new additions or updates.
    """
    # Read excel file
    df = pd.read_excel("flaskr/data/DCRT_gene_list.xlsx", sheet_name='Genes', engine='openpyxl')
    try:
        elig = str(df.loc[int(df[df['Entity Name'] == gene_symbol].index[0]), 'Eligibility'])
    except:
        elig = 'N\A'

    return elig

def get_uniprot_info(ENSG_gene_id, r_exon_skip):
    """
    This function provides the UniProt link (TO DO: add credits) and exon-specific domain information
    Input: ENSG gene identifier (without version number)
    Output: UniProt protein link, exon-specific domain information, protein name and symbol
    """
    try:
        req = requests.get(f'https://mygene.info/v3/gene/{ENSG_gene_id}?fields=uniprot')
        data = json.loads(req.content)
        uniprot_id = data['uniprot']['Swiss-Prot']
        # Get UniProt gene link
        uniprot_link = f'https://www.uniprot.org/uniprotkb/{uniprot_id}/entry'
    except:
        uniprot_id = 'N/A'
        uniprot_link = 'N/A'

    print(f"r_exon_skip: {r_exon_skip}")
    # Get the exon coordinates in the protein sequence
    var1, var2 = re.findall(r'\d+', r_exon_skip.split(':')[1])
    begin = round(int(var1)/3)
    end = round(int(var2)/3)

    # Get general domain info
    domain_info = ''
    file = urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.xml')
    data = file.read()
    file.close()
    # Convert xml to dict
    uniprot_dict = xmltodict.parse(data)

    # To print the protein name
    try:
        prot_name = (uniprot_dict['uniprot']['entry']['protein']['recommendedName']['fullName'])
        if isinstance(prot_name, str) == False:
            prot_name = (uniprot_dict['uniprot']['entry']['protein']['recommendedName']['fullName']['#text'])
    except:
        prot_name = 'N/A'

    try:
        short_name = (uniprot_dict['uniprot']['entry']['protein']['recommendedName']['shortName'])
        if isinstance(short_name, str) == False:
            short_name = (uniprot_dict['uniprot']['entry']['protein']['recommendedName']['shortName']['#text'])
    except:
        short_name = 'N/A'

    # Retrieve exon-specific domain information (need to figure out how to access evidence links)
    for entry in uniprot_dict['uniprot']['entry']['feature']:
        try:
            start_exon_p = entry['location']['begin']['@position']
            end_exon_p = entry['location']['end']['@position']
            if not(int(start_exon_p) < int(begin) and int(end_exon_p) < int(begin)):
                if not(int(start_exon_p) > int(end) and int(end_exon_p) > int(end)):
                    domain_info += entry['@description'] + ': from ' + start_exon_p + ' to ' + end_exon_p + '\n\t\t\t\t\t\t'
        except:
            continue

    return uniprot_link, domain_info, prot_name, short_name


def get_gtexportal_json(ENSG_gene_id):
    """
    This function facilitates checking for which ENSG gene id version GTEx Portal data is available
    Input: ENSG gene id with version number
    Output: GTEx Portal data
    """
    url_gtexportal = f'https://gtexportal.org/rest/v1/expression/medianTranscriptExpression?datasetId=gtex_v8 \
                        &gencodeId={ENSG_gene_id}&format=json'
    r_gtex = requests.get(url_gtexportal)
    gtex_data = json.loads(r_gtex.text)
    return gtex_data


def get_gene_expression(ENSG_gene_id, MANE_select_ENST_variant):
    """
    This function checks in which (non-eye) tissues the MANE select transcript is expressed
    It uses the GTExPortal 'Bulk tissue gene expression' data (TO DO: add credits)
    Input: ENSG gene id with version number and the MANE select ENST variant
    Output: expression levels (yes/no) of non-eye tissues

    GTEx Portal has not always included the latest ENSG version. Therefore, select data from the
    ENSG version that is present
    """
    ENSG_version = 0
    ENST_without_version = MANE_select_ENST_variant.split('.')[0] + '.'

    while get_gtexportal_json(ENSG_gene_id + '.' + str(ENSG_version))[
        'medianTranscriptExpression'] == []:
        ENSG_version += 1

    # IMPLEMENT THIS LATER IN THE TOOL
    # TO DO: find a way to inform the user that the latest ENSG version is not used to employ GTEx Portal
    # if ENSG_gene_id != ENSG_gene_id + '.' + str(ENSG_version):
    #     print('\n***!WARNING!***\n' + ENSG_gene_id + '.' + str(
    #         ENSG_version) + ' is utilized instead of ' + ENSG_gene_id +
    #           ' to consult GTEx portal\n')

    gtex_data = get_gtexportal_json(ENSG_gene_id + '.' + str(ENSG_version))

    # Set expression level to 'no' unless expression is shown
    expression_brain = 'no'
    expression_fibroblasts = 'no'
    expression_tibial_nerve = 'no'
    expression_blood = 'no'
    expression_transformed_lymphocytes = 'no'

    # Check if the MANE select transcript is expressed in the following tissues
    for transcript_in_tissue in gtex_data['medianTranscriptExpression']:
        try:
            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                    'Brain' in transcript_in_tissue['tissueSiteDetailId'] and \
                    transcript_in_tissue['median'] != 0:
                expression_brain = 'yes'

            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                    transcript_in_tissue['tissueSiteDetailId'] == 'Cells_Cultured_fibroblasts' and \
                    transcript_in_tissue['median'] != 0:
                expression_fibroblasts = 'yes'

            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                    transcript_in_tissue['tissueSiteDetailId'] == 'Nerve_Tibial' and \
                    transcript_in_tissue['median'] != 0:
                expression_tibial_nerve = 'yes'

            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                    transcript_in_tissue['tissueSiteDetailId'] == 'Whole_Blood' and \
                    transcript_in_tissue['median'] != 0:
                expression_blood = 'yes'

            if ENST_without_version in transcript_in_tissue['transcriptId'] and \
                    transcript_in_tissue['tissueSiteDetailId'] == 'Cells_EBV-transformed_lymphocytes' and \
                    transcript_in_tissue['median'] != 0:
                expression_transformed_lymphocytes = 'yes'

        # If data is not available, set expression levels to N/A
        except:
            expression_brain = 'N/A'
            expression_fibroblasts = 'N/A'
            expression_tibial_nerve = 'N/A'
            expression_blood = 'N/A'
            expression_transformed_lymphocytes = 'N/A'

    return expression_brain, \
           expression_fibroblasts, \
           expression_tibial_nerve, \
           expression_blood, \
           expression_transformed_lymphocytes

def get_eye_expression(ENSG_gene_id):
    """
    This function checks whether the gene is expressed in the periphery retina and the center retina
    Normalized eye expression data is retrieved from the Human Eye Transcriptome Atlas (TO DO: add credits)
    When two out of three samples of either periphery or center data show expression, expression is assumed to be true
    Input: ENSG gene identifier
    Output: periphery and center retina expressions (Yes/No/N/A)
    """

    query_gene = ENSG_gene_id.split('.')[0]
    eye_data = open('flaskr/data/retina_data_threshold_implemented.csv').readlines()

    periphery_retina_expression = 'N/A'
    center_retina_expression = 'N/A'

    for line in eye_data:
        gene = line.split(',')[0]
        if gene == query_gene:
            periphery = line.split(',')[1]
            center = line.split(',')[2]
            if periphery == 'True':
                periphery_retina_expression = 'yes'
            else:
                periphery_retina_expression = 'no'

            if center == 'True':
                center_retina_expression = 'yes'
            else:
                center_retina_expression = 'no'
    return periphery_retina_expression, center_retina_expression


def get_lovd_info(hg38_variant, NC_variant):
    """
    This function checks if there are exact matches available in the LOVD database (TO DO: add credits)
    Currently it's only checking the hg19 based database, this should be elaborated with hg18 and hg17.
    Input: the variant in hg38 format, the variant with NC as reference
    Output: A string containing the number of hits per found gene and the corresponding link
    TO DO: IMPROVE THIS OUTPUT FORMAT, THINK OF A WAY TO CONVENIENTLY STORE THIS IN A SQL FORMAT
    """

    exact_lovd_match_link = 'N/A'

    # First get the coordinates in the right format
    chromosome_of_variant = hg38_variant.split('-')[0]
    hg38_coordinates_for_gene_lovd = reformat_hg38_positions(NC_variant)
    hg38_coordinates_for_general_lovd = 'chr' + chromosome_of_variant + ':' + hg38_coordinates_for_gene_lovd[2:]

    # Find in which genes exact matches are found
    try:
    # Changing the nm_accession ID edit
        genes_containing_exact_hits = {}
        url = f'http://lovd.nl/search.php?build=hg38&position={hg38_coordinates_for_general_lovd}'
        url_data = pd.read_table(url, sep='\t')

        for i, gene in enumerate(url_data['gene_id']):
            genes_containing_exact_hits[gene] = url_data["DNA"][i]

        # Remove duplicates if any
        # genes_containing_exact_hits = list(dict.fromkeys(genes_containing_exact_hits))

    except:
        genes_containing_exact_hits = 'N/A'

    # Loop over the genes containing hits
    if genes_containing_exact_hits != 'N/A':
        output_exact_hits = ''

        for gene in genes_containing_exact_hits:
            # Start counting from zero again when querying a new gene
            number_exact_lovd_matches = 0
            # Get nm_accession id from dictionary for corresponding gene
            dna_id = genes_containing_exact_hits.get(gene)
            # Check if variant position EXACTLY matches other variants
            try:
                req_exact = requests.get(f'http://databases.lovd.nl/shared/api/rest.php/variants/{gene}?search_position={dna_id}&format=application/json')
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

    else:
        output_exact_hits = 'N/A'

    return output_exact_hits

# hg19
def reformat_hg19_positions(NC_variant):
    """
    This function converts the hg19 genomic description (i.e. NC_000003.11:g.48612270_48612271del) to a format
    that is accepted by LOVD (i.e. g.48612270_48612271). Note that chromosome information is not needed since
    LOVD already known to which chromosome the gene belongs
    Input: hg19 genomic description
    Output: LOVD conform hg19 description
    """
    try:
        hg19_coordinates = NC_variant.split('.')[-1].split('_')

        if len(hg19_coordinates) == 1:
            coordinate1 = ''.join([i for i in hg19_coordinates[0] if i.isdigit()])
            hg19_coordinates_for_gene_based_lovd = 'g.' + coordinate1
        else:
            coordinate1 = hg19_coordinates[0]
            coordinate2 = ''.join([i for i in hg19_coordinates[1] if i.isdigit()])
            hg19_coordinates_for_gene_based_lovd = 'g.' + coordinate1 + '_' + coordinate2
    except:
        hg19_coordinates_for_gene_based_lovd = 'N/A'

    return hg19_coordinates_for_gene_based_lovd

def get_lovd_info_hg19(hg19_variant, NC_variant):
    """
    This function checks if there are exact matches available in the LOVD database (TO DO: add credits)
    Currently it's only checking the hg19 based database, this should be elaborated with hg18 and hg17.
    Input: the variant in hg38 format, the variant with NC as reference
    Output: A string containing the number of hits per found gene and the corresponding link
    TO DO: IMPROVE THIS OUTPUT FORMAT, THINK OF A WAY TO CONVENIENTLY STORE THIS IN A SQL FORMAT
    """

    exact_lovd_match_link = 'N/A'

    # First get the coordinates in the right format
    chromosome_of_variant = hg19_variant.split('-')[0]
    hg19_coordinates_for_gene_lovd = reformat_hg19_positions(NC_variant)
    hg19_coordinates_for_general_lovd = 'chr' + chromosome_of_variant + ':' + hg19_coordinates_for_gene_lovd[2:]

    # Find in which genes exact matches are found
    try:
    # Changing the nm_accession ID edit
        genes_containing_exact_hits = {}
        url = f'http://lovd.nl/search.php?build=hg19&position={hg19_coordinates_for_general_lovd}'
        url_data = pd.read_table(url, sep='\t')

        for i, gene in enumerate(url_data['gene_id']):
            genes_containing_exact_hits[gene] = url_data["DNA"][i]

        # Remove duplicates if any
        # genes_containing_exact_hits = list(dict.fromkeys(genes_containing_exact_hits))

    except:
        genes_containing_exact_hits = 'N/A'

    # Loop over the genes containing hits
    if genes_containing_exact_hits != 'N/A':
        output_exact_hits = ''

        for gene in genes_containing_exact_hits:
            # Start counting from zero again when querying a new gene
            number_exact_lovd_matches = 0
            # Get nm_accession id from dictionary for corresponding gene
            dna_id = genes_containing_exact_hits.get(gene)
            # Check if variant position EXACTLY matches other variants
            try:
                req_exact = requests.get(f'http://databases.lovd.nl/shared/api/rest.php/variants/{gene}?search_position={dna_id}&format=application/json')
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

    else:
        output_exact_hits = 'N/A'

    return output_exact_hits
