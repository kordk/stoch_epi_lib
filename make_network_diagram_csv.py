##
##
##  Given gene and site files (i.e. parameters for go program), creates 2 .csv files, "network_picture_edges.csv" and one with "network_picture_nodes.csv" (or "optional_name_edges.csv",
##      "optional_name_nodes.csv"). Import these into cytoscape using network from file -> edges.csv, then table from file -> nodes.csv.
##      Network can be styled in cytoscape.
##





import sys
import pandas as pd
import numpy as  np
import json

gnfl = sys.argv[1]
sifl = sys.argv[2]


try:
    fl_save = sys.argv[3]
except:
    fl_save = "network_picture"

with open(gnfl) as fl:
    genes = json.load(fl)

with open(sifl) as fl:
    sites = json.load(fl)

edge_table = pd.DataFrame(columns = ["Source","Target","Type","Val"])
for site in sites:
    src = genes[site["Kaploc"][0]]['GeneID']
    targ = src + site["Aid"]
    siedge = pd.DataFrame([[src,targ,'Kappa',1]],columns = ["Source","Target","Type","Val"])
    edge_table = edge_table.append(siedge,ignore_index=True)

for gn in genes:
    for jj in range(len(gn["Philoc"])):
        j = gn["Philoc"][jj]
        st = sites[j]
        src0 = genes[st["Kaploc"][0]]['GeneID']
        src = src0 + st["Aid"]
        gnedge =  pd.DataFrame([[src,gn["GeneID"],'Phi',gn["Phival"][jj]]],columns = ["Source","Target","Type","Val"])
        edge_table = edge_table.append(gnedge,ignore_index=True)

edge_table.to_csv(fl_save+'_edges.csv',index = False)

genenames = [gn["GeneID"] for gn in genes]
nodenms = np.unique(np.concatenate([edge_table.loc[:,'Source'].values,edge_table.loc[:,'Target'].values]))
picture_nodes = pd.DataFrame(columns = ['Type'], index = nodenms)

for nd in nodenms:
    if nd in genenames:
        picture_nodes.loc[nd,'Type'] = 'Gene'
    else:
        picture_nodes.loc[nd,'Type'] = 'Site'

picture_nodes.to_csv(fl_save+'_nodes.csv')
