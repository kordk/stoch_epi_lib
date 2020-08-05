
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

gene_gene_file = sys.argv[1]#'interaction_parameters/toy_network.csv'###Next thing to do is make this generic
epigen_file = sys.argv[2]#'interaction_parameters/toy_site_agg.csv'

####Read in data
gene_gene = pd.read_csv(gene_gene_file)
epigen = pd.read_csv(epigen_file, sep = '\t')


####Build model. First, get the list of genes
gene_ids = pd.concat([gene_gene.Target, gene_gene.Regulator]).unique()

####Every arrow of GRN is duplicated to number of binding sites associated with target.
#### So we need to assign sites to all the interactions they are in the middle of.
site_names = epigen.loc[:,'region'].unique()

sites = pd.Series(index = gene_gene.index)
for i in gene_gene.index:
    sites[i] = epigen[epigen.loc[:,'gene']==gene_gene.loc[i,'Target']].loc[:,'region'].values

gene_gene.loc[:,'Sites'] = sites

site_inds = pd.Series(index = gene_gene.index,dtype = 'object')
site_score = pd.Series(index = gene_gene.index,dtype = 'object')

site_to_ind = {}
ind_to_site = {}
ct = 0
for i in gene_gene.index:
    sites_here = []
    for si in gene_gene.loc[i,'Sites']:
        if gene_gene.loc[i,'Regulator'] + si in site_to_ind.keys():
            sites_here += [site_to_ind[gene_gene.loc[i,'Regulator'] + si]]
        else:
            sites_here += [ct]
            site_to_ind[gene_gene.loc[i,'Regulator'] + si] = ct
            ind_to_site[ct] = si
            ct += 1
    site_inds[i] = sites_here
    site_score[i] = [gene_gene.loc[i,'IntScore']]*len(site_inds[i])


gene_gene.loc[:,'site_inds'] = site_inds
gene_gene.loc[:,'site_score'] = site_score


##In our model we have a vector B for ``sites" that control a regulation (binding sites). These are
## artificial because we duplicate them.
####
###     This script assumes that model is built from the gene-gene interactions, and sites
###     duplicated for each arrow of the GRN. That is, if we have the GRN A->B,A->C,B->C and
###     sites x,y,z associated with B, and x,w associated with C, we have 6 "arrows":
###         x_(A),y_(A),z_(A),w_(A),x_(B),w_(B)
###
###
##
###
###         Also we aren't considering the effect of competition for binding site.
###
##
## Here, we have a vector of the ``real" site associated with each entry in that vector.
###

phis_ids = [np.concatenate(gene_gene[gene_gene.Target == trg].loc[:,'site_inds'].values) if trg in gene_gene.Target.unique() else np.array([],dtype = 'int') for trg in gene_ids]


print(max(np.concatenate(phis_ids)))


phis_sc = [np.concatenate(gene_gene[gene_gene.Target == trg].loc[:,'site_score'].values) if trg in gene_gene.Target.unique() else np.array([],dtype = 'int') for trg in gene_ids]


# arrow_sites = np.concatenate([li for li in gene_gene.loc[:,'Sites'].values])
arrow_sites = [ind_to_site[i] for i in range(len(ind_to_site))]

num_arrows = len(arrow_sites)

print(num_arrows)

print(all(np.concatenate([li for li in gene_gene.loc[:,'Sites'].values]) == arrow_sites))



num_genes = len(gene_ids)

#####



min_grth = [-float(sum(np.minimum(ps,0))) for ps in phis_sc]




####For each site variable, we need to know which gene binds it. We leave this as a list
### To accomodate the more general case in which we have a site variable bound by more than
#### one regulator.

binder = np.concatenate([[gene_gene.loc[i,'Regulator']]*len(gene_gene.loc[i,'Sites']) for i in gene_gene.index])
binder_inds = [np.argwhere(gene_ids.astype('str') == gen)[0].tolist() for gen in binder]

#### Choose some parameters randomly to generate some fake data. Or load from file (pickle file)

NewParams = True

if NewParams:
    site_params_rand = pd.DataFrame(rand(num_arrows,2), columns = ['lambda','lambda_hat'])
    epigen_params_rand = pd.DataFrame(rand(len(site_names),3) + np.array(len(site_names)*[[0,-0.5,0]]), columns = ['mu','nu','alpha'], index = site_names)
    gene_params_rand = pd.DataFrame(rand(num_genes,2), columns = ['gamma','decay'])



    ##need to make sure we can't have a negative amount of transcript. Down regulate reduces production rather
    ## than increasing degredation.
    ##
    ##Comment out if loading parameters from pickle
    ##

    gene_params_rand.gamma = gene_params_rand.gamma + min_grth ###
    gene_params_rand.decay = gene_params_rand.decay + np.ones(len(gene_params_rand.decay))#just so computation is so slow

    ##Comment out if loading parameters from pickle
    site_params_rand.loc[:,'mu'] = [epigen_params_rand.loc[arrow_sites[i],'mu'] for i in site_params_rand.index]
    site_params_rand.loc[:,'nu'] = [epigen_params_rand.loc[arrow_sites[i],'nu'] for i in site_params_rand.index]
    site_params_rand.loc[:,'alpha'] = [epigen_params_rand.loc[arrow_sites[i],'alpha'] for i in site_params_rand.index]


else:
    #### Load parameters from pickle so they're the same
    site_params_rand = pd.read_pickle('saved_parameters/random_site_parameters.pkl')
    epigen_params_rand = pd.read_pickle('saved_parameters/random_epigenetic_parameters.pkl')
    gene_params_rand = pd.read_pickle('saved_parameters/random_gene_parameters.pkl')

SaveParams = False
if SaveParams:
    ###Comment out if loading parameters from pickle
    site_params_rand.to_pickle('random_site_parameters.pkl')
    epigen_params_rand.to_pickle('random_epigenetic_parameters.pkl')
    gene_params_rand.to_pickle('saved_parameters/random_gene_parameters.pkl')




#grab the ``alphas" - these are the epigenetic modification scores from data.
alphas = site_params_rand.loc[:,'alpha'].values

# gene_params_rand['gamma'] = min_grth

####### Save parameters as JSON files for go program to run the model (if new paramerters.)

sites_toy = site_params_rand.copy()
sites_toy.columns = ['Lambda','Lambdahat','Mu','Nu','Alpha']
stdict = sites_toy.to_dict(orient='list')
stdict['Kaploc'] = binder_inds
stdict['Kapval'] = [[1]*len(binder_inds[i]) for i in range(len(binder_inds))]

dlist = [sites_toy.loc[i].to_dict() for i in sites_toy.index]
for i in range(len(dlist)):
    dlist[i]['Kaploc'] = binder_inds[i]
    dlist[i]['Kapval'] = [1]*len(binder_inds[i])
    dlist[i]['Aid'] = arrow_sites[i]


with open('sites.json', 'w') as handle:
    json.dump(dlist, handle)


toy_genes = gene_params_rand.copy()
toy_genes.columns = ['Gamma','Dil']
gene_dicts = [toy_genes.loc[i].to_dict() for i in toy_genes.index]

for i in range(len(gene_dicts)):
    gene_dicts[i]['Philoc'] = phis_ids[i].tolist()
    gene_dicts[i]['Phival'] = phis_sc[i].tolist()
    gene_dicts[i]['GeneID'] = gene_ids[i]

with open('genes.json','w') as handle:
    json.dump(gene_dicts,handle)
######################
