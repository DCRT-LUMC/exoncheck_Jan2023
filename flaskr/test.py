import json

import requests
MANE_select_NM_variant = 'NM_004006.2:c.4634dup'
NM_id = 'NM_004006.2'
latest_reference_sequence = 'NC_000023.11'
exon_number = '33'
total_exons = '79'

req2 = requests.get(f'https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts_v2/{MANE_select_NM_variant}/{NM_id}?content-type=application%2Fjson')
data2 = json.loads(req2.content)

for exon in data2["transcripts"][0]["genomic_spans"][latest_reference_sequence]["exon_structure"]:
    if str(exon["exon_number"]) == exon_number:
        genomic_end = str(exon["genomic_end"])
        genomic_start = str(exon["genomic_start"])
        exon_length = int(exon["cigar"][:-1])
        NC_exon_NC_format = latest_reference_sequence + ':g.' + genomic_start + '_' + genomic_end + 'del'
        if exon_number == '1' or exon_number == total_exons:
            exon_length = int(exon_length) - int(data2["transcripts"][0]["coding_start"]) + 1
        elif exon_number == total_exons:
            exon_length = int(exon_length) - int(data2["transcripts"][0]["coding_end"]) + 1
print(exon_length)

try:
    exon_length /= 3.0
except:
    exon_length = 0

print(exon_length)


if exon_length.is_integer():
    frame = 'In-frame'
else:
    frame = 'Out-of-frame'

print(frame)