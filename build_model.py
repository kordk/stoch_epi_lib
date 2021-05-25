
####
####        Script to generate model using ``site splitting"
####

###     Builds structure of model as gene_gene dataframe
###     site_inds is index of site associated with the gene-gene interaction (based on target)
###        site_score is the IntScore but expanded for each site.
###     Compiles this into
####        phis_ids & phis_sc
###
###         phis_ids is a list of lists - list j contains the index of sites that regulate gene j
###         phis_sc is a list of lists - list j contains the ``scores" (+/- 1) of the regulation
####        of gene j by sites in list phis_ids[j]
###                     so phi_j \cdot B becomes dot(B[phis_ids[j]],phis_sc[j])

###
###         binder_inds:
###                 This is a list of lists - list i contains the genes that bind to site i
###                             (by index). In this model construction, each site is bound by
###                             a single gene. Epigenetic data is repeated for sites that are
###                             duplicated.


####        arrow_sites - matches model sites with real sites - in case sites are split.


#### Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import fsolve
from numpy.random import rand
import seaborn as sb
import pandas as pd
import time
from joblib import Parallel, cpu_count, delayed

import json

from cycler import cycler

import sys

from grn_with_epi_functions import *



###### READ IN AND CONSTRUCT THE MODEL
print("Reading Network Topology and Region Info")


#Our underlying (full) GRN should be available as a csv file with columns:
# ["Regulator","Target","IntScore"] (others will be ignored) populated with [geneID,geneID,strength and direction of interaction]

#The binding site information should be available in .tsv (separated by tabs unless specified) with columns
# ["gene","region"] (others will be ignored) populated with [downstream geneID, region ID]


try:
    #First look for paths to the gene-gene network file and binding site information file in a
    #txt file. Format should be:
    #region_info=PATH
    #gene_gene_network=PATH
    #on separate lines (order does not matter)
    with open("data_file_paths.txt") as fl:
        lines = fl.readlines()

    filelocs = dict([ln.replace('\n','').split("=") for ln in lines])

    try:
        region_info_df = pd.read_csv(filelocs["region_info"], sep = sys.argv[1])
    except:
        region_info_df = pd.read_csv(filelocs["region_info"], sep = '\t')


    gene_gene_network_df = pd.read_csv(filelocs["gene_gene_network"])


except:
    #if the file path information is not in a txt file, it can be given as ordered options
    gene_gene_file = sys.argv[1]#Path to GRN information
    epigen_file = sys.argv[2]#Path to region info

    gene_gene_network_df = pd.read_csv(gene_gene_file)

    try:
        region_info_df = pd.read_csv(epigen_file, sep = sys.argv[3])
    except:
        region_info_df = pd.read_csv(epigen_file, sep = '\t')


### look for a list of genes to exclude - given as txt file with one gene ID per line
try:
    with open("leave_out_genes.txt") as fl:
        lv_out_genes = [ln.replace("\n",'') for ln in fl.readlines()]
except:
    lv_out_genes = []

### look for a list of sites to exclude - given as txt file with one region ID per line
try:
    with open("leave_out_sites.txt") as fl:
        lv_out_sites = [ln.replace("\n",'') for ln in fl.readlines()]
except:
    lv_out_sites = []




print("Creating network structure")

# ####Every arrow of GRN is duplicated to number of binding sites associated with target.
# #### So we need to assign sites to all the interactions they are in the middle of.

gene_gene_network_df = gene_gene_network_df[[(gn not in lv_out_genes) for gn in gene_gene_network_df["Regulator"].values]]
gene_gene_network_df = gene_gene_network_df[[(gn not in lv_out_genes) for gn in gene_gene_network_df["Target"].values]]

#
sites = pd.Series(index = gene_gene_network_df.index,dtype = 'object')
for i in gene_gene_network_df.index:
    rawSites = region_info_df[region_info_df.loc[:,'gene']==gene_gene_network_df.loc[i,'Target']].loc[:,'region'].values
    sites[i] = [si for si in rawSites if si not in lv_out_sites]

gene_gene_network_df.loc[:,'Sites'] = sites

initial_len  = len(gene_gene_network_df)

print("Removing interactions with no site information")

gene_gene_network_df = gene_gene_network_df[gene_gene_network_df["Sites"].map(len)>0]

removed = initial_len - len(gene_gene_network_df)
print("Removed {0} interactions due to lack of  site information".format(removed))

####Build model. First, get the list of genes
gene_ids = pd.concat([gene_gene_network_df.Target, gene_gene_network_df.Regulator]).unique()




# ##In our model we have a vector B for ``sites" that control a regulation (binding sites). These are
# ## artificial because we duplicate them.
# ####
# ###     This script assumes that model is built from the gene-gene interactions, and sites
# ###     duplicated for each arrow of the GRN. That is, if we have the GRN A->B,A->C,B->C and
# ###     sites x,y,z associated with B, and x,w associated with C, we have 7 "arrows":
# ###         x_(AB),y_(AB),z_(AB),x_(AC),w_(AC),x_(BC),w_(BC)
# ###
# ####      NOTE: The mdoel fitting and run can handle more complex bipartite graph topology.
# ##
# ###           Gene parameters are a list of dicts corresponding to each gene. Dict entries must include dynamic info:
# ##                - GeneID
# ##                - Dil (dilution rate)
# ##                - Gamma (intrinsic preduction rate.)
# ##                And STRUCTURAL INFO
# ##                - Philoc (list of site indices such that we have Site -> gene arrow in the model)
# ##                - Phival (list of Corresponding strength of interactions, signed to allow repression)
# ##
# ##            Site parameters are a list of dicts corresponding to each site in the network. In this
# ##            script, each actual site is copied for each gene(regulator)-gene(target) interaction
# ##            such that the actual site is associated with the target. Dict entries must incliude dynamic info:
# ##                - Lambda (Intrinsic activation propensity)
# ##                - Lambdahat (Intrinsic deactivation propensity)
# ##                - Mu (k0 parameter for epigenetic effect on activation)
# ##                - Nu (hill power for epigenetic effect  on activation)
# ##                - Alpha  (epigenetic parameter)
# ##                - Aid (actual Region ID)
# ##                And STRUCTURAL INFO
# ##                - Kaploc (index of genes such that we have Gene -> site arrow in the model)
# ##                - Kapval (corresponding strength of interaction)
# ##
# ##
# ##           THUS: Network topology is determined by the set of Philoc/Phival pairs and Kaploc/Kapval pairs

#### To construct our model files, we will simply interate throught the set  of genes.


expanded_gene_gene_network_df = pd.DataFrame(columns = gene_gene_network_df.columns)
for gngn in gene_gene_network_df.index:
    for siteID in gene_gene_network_df.loc[gngn,'Sites']:
        row = pd.DataFrame([list(gene_gene_network_df.loc[gngn].values[:-1])  + [siteID]],columns = gene_gene_network_df.columns)
        expanded_gene_gene_network_df = expanded_gene_gene_network_df.append(row, ignore_index = True)


gene_params = []
gene_indices = {}

print("Creating gene information")
for gn in gene_ids:
    gene_indices[gn] = len(gene_params)
    gdict = {"GeneID":gn,"Dil":np.random.rand()}
    phis = list(expanded_gene_gene_network_df.index[expanded_gene_gene_network_df['Target']==gn])
    gdict["Philoc"] = phis
    gdict["Phival"] = list(expanded_gene_gene_network_df.loc[phis,"IntScore"].values)
    gdict["Gamma"] = -sum([min(0,x) for  x in gdict["Phival"]]) + 1
    # if len(phis) == 0:
    #     gdict["Gamma"] = 1
    gene_params +=  [gdict]

print("Model built with {0} genes".format(len(gene_params)))

print("Creating site information")
site_params = []
for edge in expanded_gene_gene_network_df.index:
    sdict = {"Aid":expanded_gene_gene_network_df.loc[edge,"Sites"],"Lambda":np.random.rand(),"Mu":np.random.rand(),"Alpha":np.random.rand(),"Lambdahat":np.random.rand(),"Nu":np.random.rand()-0.5}
    kaps = [gene_indices[expanded_gene_gene_network_df.loc[edge,"Regulator"]]]
    sdict["Kaploc"] = kaps
    sdict["Kapval"] = [1]
    site_params  += [sdict]

print("Model built with {0} sites".format(len(site_params)))

print("Saving as .json")
with open('genes.json','w') as handle:
    json.dump(gene_params,handle)
#
with open('sites.json', 'w') as handle:
    json.dump(site_params, handle)

#

######################
