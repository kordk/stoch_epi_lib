# Frequently Asked Quesitons (FAQ)

#### 1. "TypeError: 'bool' object is not iterable" when creating the gene regulatory network with build_model.py.  

<pre>
build_model.py:129: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  print(all(np.concatenate([li for li in gene_gene.loc[:,'Sites'].values]) == arrow_sites))
611.0
612
Traceback (most recent call last):
  File "build_model.py", line 129, in <module>
    print(all(np.concatenate([li for li in gene_gene.loc[:,'Sites'].values]) == arrow_sites))
TypeError: 'bool' object is not iterable
</pre>

I think what’s happening is we have a regulator binding to a site and effecting two genes. That print was checking for that, but apparently doing a bad job of it. And it hadn’t found that before. You can remove that print statement or replace with:
<pre>
print([si for si in allsites if si not in arrow_sites])
print(len(allsites)==len(arrow_sites))
</pre>

#### 2. The ParamFit method is returning an infinite -log likelihood and stops. 
<pre>
[MasterFit] Computing initial likelihood.
[MasterFit] Initial likelihood: +Inf No jump on 973.65641/1000 draws.
[MasterFit] Generating sample parameter sets.
[MasterFit] Stopping condition is %d steps in a row with descrease of less than 10^%d % 10 -6
[MasterFit] 0/1000 Running. 2020.08.06 17:20:11
[MasterFit] l2 norm of Gradient is now NaN, Stopping.
</pre>

We’ll get infinite -log likelihood whenever the model never gets close to the data, which becomes much more likely as we increase the network size. 

First, generate an equilibrium distribution with random parameter for a few samples to see how far away it is. This should also help us catch anything weird going on with the data scaling. 

Then, try cranking up “EstimateLength” and “Distband”. EstimateLength will let the model run longer and so hopefully get near the data occasionally, Distband will widen what we mean by “near”. 

