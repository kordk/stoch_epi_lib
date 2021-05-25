import numpy as np
import json


with open("genes.json") as fl:
    intitial_genes = json.load(fl)

with open("sites.json") as fl:
    initial_sites = json.load(fl)

with open("fitted_genes.json") as fl:
    fitted_genes = json.load(fl)

with open("fitted_sites.json") as fl:
    fittted_sites = json.load(fl)

diff_genes = []
for i in range(len(intitial_genes)):
    diffs = {}
    for ky in ["Dil","Gamma"]:
        diffs[ky] = (intitial_genes[i][ky]-fitted_genes[i][ky])**2
    diff_genes += [diffs]

diff_sites = []
for i in range(len(initial_sites)):
    diffs = {}
    for ky in ["Lambda","Lambdahat","Mu","Nu"]:
        diffs[ky] = (initial_sites[i][ky]-fittted_sites[i][ky])**2
    diff_sites += [diffs]


sum1 = 0
sum2 = 0
for i in range(len(diff_genes)):
    sum1 += diff_genes[i]["Dil"]
    sum2 += diff_genes[i]["Gamma"]
mean_diff_genes = {"Dil":sum1**(1/2)/len(diff_genes), "Gamma" : sum2**(1/2)/len(diff_genes)}

sum1 = 0
sum2 = 0
sum3 = 0
sum4 = 0
for i in range(len(diff_sites)):
    sum1 += diff_sites[i]["Lambda"]
    sum2 += diff_sites[i]["Lambdahat"]
    sum3 += diff_sites[i]["Mu"]
    sum4 += diff_sites[i]["Nu"]
mean_diff_sites = {"Lambda":sum1**(1/2)/len(diff_sites), "Lambdahat" : sum2**(1/2)/len(diff_sites), "Mu" : sum3**(1/2)/len(diff_sites), "Nu" : sum3**(1/2)/len(diff_sites)}

print(mean_diff_genes)
print(mean_diff_sites)
