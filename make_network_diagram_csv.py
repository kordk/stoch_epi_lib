##
##
##  Given same inputs as "build model", creates 2 .csv files, "network_picture_edges.csv" and one with "network_picture_nodes.csv" (or "optional_name_edges.csv",
##      "optional_name_nodes.csv"). Import these into cytoscape using network from file -> edges.csv, then table from file -> nodes.csv.
##      Network can be styled in cytoscape.
##





import sys
import pandas as pd
import numpy as  np

gene_gene_file = sys.argv[1]
epigen_file = sys.argv[2]#'interaction_parameters/toy_site_agg.csv'


if len(sys.argv) >3:
    fl_save = sys.argv[3]
else:
    fl_save = "network_picture"

gene_gene = pd.read_csv(gene_gene_file)
epigen = pd.read_csv(epigen_file, sep = '\t')


####Build model. First, get the list of genes
gene_ids = pd.concat([gene_gene.Target, gene_gene.Regulator]).unique()

####Every arrow of GRN is duplicated to number of binding sites associated with target.
#### So we need to assign sites to all the interactions they are in the middle of.
sites = pd.Series(index = gene_gene.index)
for i in gene_gene.index:
    sites[i] = epigen[epigen.loc[:,'gene']==gene_gene.loc[i,'Target']].loc[:,'region'].values

gene_gene.loc[:,'Sites'] = sites

site_inds = pd.Series(index = gene_gene.index,dtype = 'object')
site_score = pd.Series(index = gene_gene.index,dtype = 'object')

site_to_ind = {}
ct = 0
for i in gene_gene.index:
    sites_here = []
    for si in gene_gene.loc[i,'Sites']:
        if gene_gene.loc[i,'Regulator'] + si in site_to_ind.keys():
            sites_here += [site_to_ind[gene_gene.loc[i,'Regulator'] + si]]
        else:
            sites_here += [ct]
            site_to_ind[gene_gene.loc[i,'Regulator'] + si] = ct
            ct += 1
    site_inds[i] = sites_here
    site_score[i] = [gene_gene.loc[i,'IntScore']]*len(site_inds[i])


gene_gene.loc[:,'site_inds'] = site_inds
gene_gene.loc[:,'site_score'] = site_score




picture_edges = pd.DataFrame(columns = ['Source','Target','Type','Val'])
for i in gene_gene.index:
    nsts = len(gene_gene.loc[i,'site_inds'])
    df1 = pd.DataFrame({'Source' : [gene_gene.loc[i,'Regulator']]*nsts,'Target' : [gene_gene.loc[i,'Regulator'] + Si for Si in gene_gene.loc[i,'Sites']],'Type' : ['Kappa']*nsts,'Val' : [1]*nsts})
    df2 = pd.DataFrame({'Source' : [gene_gene.loc[i,'Regulator'] + Si for Si in gene_gene.loc[i,'Sites']],'Target' : [gene_gene.loc[i,'Target']]*nsts,'Type' : ['Phi']*nsts,'Val' : gene_gene.loc[i,'site_score']})
    picture_edges = picture_edges.append(df1,ignore_index = True)
    picture_edges = picture_edges.append(df2,ignore_index = True)

picture_edges.to_csv(fl_save + '_edges.csv', index = False)

genenames = np.unique(np.concatenate([gene_gene.loc[:,'Regulator'].values,gene_gene.loc[:,'Target'].values]))
nodenms = np.unique(np.concatenate([picture_edges.loc[:,'Source'].values,picture_edges.loc[:,'Target'].values]))
picture_nodes = pd.DataFrame(columns = ['Type'], index = nodenms)

for nd in nodenms:
    if nd in genenames:
        picture_nodes.loc[nd,'Type'] = 'Gene'
    else:
        picture_nodes.loc[nd,'Type'] = 'Site'

picture_nodes.to_csv(fl_save+'_nodes.csv')
