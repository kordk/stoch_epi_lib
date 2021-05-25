Store the path to the network and site list in data_file_paths and run:

$ python ../build_model.py

OUTPUT:
/Users/m197894/anaconda3/lib/python3.7/site-packages/sklearn/utils/deprecation.py:143: FutureWarning: The sklearn.neighbors.kde module is  deprecated in version 0.22 and will be removed in version 0.24. The corresponding classes / functions should instead be imported from sklearn.neighbors. Anything that cannot be imported from sklearn.neighbors is now part of the private API.
  warnings.warn(message, FutureWarning)
Reading Network Topology and Region Info
Creating network structure
Removing interactions with no site information
Removed 0 interactions due to lack of  site information
Creating gene information
Model built with 3 genes
Creating site information
Model built with 4 sites
Saving as .json

This will create genes.json and sites.json, which contain (with random-parameters) a representation of the model readable by the executable.

We may then run:

$ python ../make_network_diagram_csv.py genes.json sites.json

to generate a network diagram in format readable by cytoscape. This can be imported into cytoscape using the "import network from file" button and selecting network_picture_edges.csv. Next, import node types using the "import table from file" button and selecting network_picture_nodes.csv.

To demonstrate parameter fitting, we can generate some fake data:

$ cd data

$ ../../stoch_epi_lib -mode=Fakeit -sitefl=../sites.json -genefl=../genes.json -svfl=generated_data -numfake=50

$ cd ..

Now, we can run a parameter fitting with:

$ ../stoch_epi_lib -mode=ParamFit -DataScale=none -datafl=data/generated_data.json -DistBand=10 -EqLength=2000 -NumberCores=-1 -sitefl=sites.json -genefl=genes.json -svfl=fitted -DistBand=1 -StoppingCondition=30

And then generate an approximate equilibrium distribution:

$ ../stoch_epi_lib -mode=Eq -DataScale=none -datafl=data/generated_data.json -EqLength=2000 -NumberCores=-1 -sitefl=fitted_sites.json -genefl=fitted_genes.json -DistBand=1 -svfl=EqOut/example 

Approximate equilibrium distribution information for sample i is given in a json file named EqOut/example_transcript_i.json. This can be imported into python:

import json
with open("EqOut/example_transcript_1.json") as fl:
    sample1 = json.load(fl)

Each element of the list corresponds to one of the gene and contains fields:
dict_keys(['GeneID', 'DistEst', 'XDomain', 'RealMeans', 'RealVars', 'TimePts', 'MeanScaled', 'StdScaled'])

We can then, for example, generate a plot of the approximate equilibrium distribution for gene g2 by:

import matplotlib.pyplot as plt
plt.plot(sample1[0]["XDomain"],sample1[0]["DistEst"])

Note that gene 1 is not included, as it is not regulated by the network (appearing only as a source node, never a target).

To see how the estimate of the mean changes as we increase the length of the simulation (to get an idea about convergence to the equilibrium distribution) we can do:

plt.plot(sample1[0]["TimePts"],sample1[0]["RealMeans"])

In order to generate trajectories of the model, we can run:

$ ../stoch_epi_lib -mode=Realization -DataScale=none -datafl=data/generated_data.json -EqLength=2000 -NumberCores=-1 -sitefl=fitted_sites.json -genefl=fitted_genes.json -DistBand=1 -svfl=realizations/example

This will generate realization file for sample i named realizations/example_even_i.json (with data at evenly spaced timepoints) and realizations/example_jumps_i.json (with data at jump times). We can inspect these realizations in python with:

import json
import numpy as np
import matplotlib.pyplot as plt

with open("realizations/example_even_0.json") as fl:
    realization1 = json.load(fl)

fig,ax = plt.subplots(figsize = (15,7))
ax.set_facecolor('#efeff0')
ax.plot(realization1['T'],np.array(realization1["G"])[:,0],label = realization1["GeneIDs"][0],color = '#1f34d4',linewidth = 4)
ax.plot(realization1['T'],np.array(realization1["G"])[:,1],label = realization1["GeneIDs"][1],color = '#ee7a05',linewidth = 4)
ax.step(realization1['T'],np.array(realization1["SG"])[:,0],color = '#1f34d4',linewidth = 1)
ax.step(realization1['T'],np.array(realization1["SG"])[:,1],color = '#ee7a05',linewidth = 1)
ax.legend(fontsize = 20)
ax.set_xlabel("Time",fontsize = 18)
ax.set_ylabel("Transcript",fontsize = 18)
ax.tick_params(labelsize=14)

