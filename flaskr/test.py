import json
import pandas as pd
import numpy as np

import requests

# This function checks whether the gene is expressed in the periphery retina and the center retina
# Normalized eye expression data is retrieved from the Human Eye Transcriptome Atlas (TO DO: add credits)
# When two out of three samples of either periphery or center data show expression, expression is assumed to be true
# Input: ENSG gene identifier
# Output: periphery and center retina expressions (Yes/No/N/A)

eye_data = open('data/retina_data_threshold_implemented.csv').readlines()
query_gene = "ENSG0000027x9602"

periphery_retina_expression = 'N/A'
center_retina_expression = 'N/A'

for line in eye_data:
    gene = line.split(',')[0]
    if gene == query_gene:
        periphery = line.split(',')[1]
        center = line.split(',')[2]
        if periphery == 'True':
            periphery_retina_expression = 'Yes'
        else:
            periphery_retina_expression = 'No'

        if center == 'True':
            center_retina_expression = 'Yes'
        else:
            center_retina_expression = 'No'

print(center_retina_expression, periphery_retina_expression)


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
