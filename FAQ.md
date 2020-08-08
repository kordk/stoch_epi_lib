# Frequently Asked Quesitons (FAQ)

1. The ParamFit method is returning an infinite -log likelihood and stops. What's happenning?
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

