import requests
import json
import csv
import argparse
import pandas as pd
url = 'https://maayanlab.cloud/L1000CDS2/query'

def upperGenes(genes):
    return [gene.upper() for gene in genes]

parser = argparse.ArgumentParser()
parser.add_argument('--TF_path1')
parser.add_argument('--TF_path2')
args = parser.parse_args()

g1 = pd.read_csv(args.TF_path1, sep='\t')
tf1 = g1.TF1.to_list()
val1 = g1.tfs_score.to_list()

g2 = pd.read_csv(args.TF_path2, sep='\t')
tf2 = g2.TF1.to_list()
val2 = g2.tfs_score.to_list()

data = {"upGenes":tf1,
"dnGenes":tf2}
data['upGenes'] = upperGenes(data['upGenes'])
data['dnGenes'] = upperGenes(data['dnGenes'])
config = {"aggravate":True,"searchMethod":"geneSet","share":True,"combination":True,"db-version":"latest"}
#metadata = [{"key":"Tag","value":"gene-set python example"},{"key":"Cell","value":"MCF7"}]
metadata = []
payload = {"data":data,"config":config,"meta":metadata}
headers = {'content-type':'application/json'}
r = requests.post(url,data=json.dumps(payload),headers=headers)
resGeneSet = r.json()

for Meta in resGeneSet['topMeta']:
    for name in Meta['pert_desc']:
        if Meta['pert_desc'] == "-666":
            Meta['pert_desc'] = "Unknown"

print('Drug_name , Pubchem_id , Score, Cell_line')
for Meta in resGeneSet['topMeta']:
    print(Meta['pert_desc'], ',', Meta['pubchem_id'], ',', Meta['score'], ',', Meta['cell_id']) 
