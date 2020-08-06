# stoch_epi_lib
stoch_epi_lib: An implementation of a dynamical systems method to estimate gene expression from epigenetic data and a gene regulatory network


This software and pipeline is to estimate gene expression using the stoch_epi_lib program. The method is described in a manuscript available in pre-print on bioRxiv Systems Biology ( https://biorxiv.org/cgi/content/short/2020.08.03.234740v1) and arXiv Molecular Networks (https://arxiv.org/abs/2008.01126) and submitted for publication at the <a href="https://psb.stanford.edu/">Pacific Symposium Biocomputing</a>. 

## Requirements

     - Go (https://golang.org/)
     - Go modules : gonum.org/v1/gonum/mat, 
                    github.com/montanaflynn/stats, 
                    gonum.org/v1/gonum/integrate/quad, 
                    github.com/aclements/go-gg/generic/slice
     - Python 2.7
     - Python modules : pandas, numpy, matplotlib, scipy, seaborn, joblib, cycler, smtplib, sklearn, simplejson
     
## Installation instructions

### Download a release

Periodically, we may provide a release of the tool as a stand-alone executable with an installation tool. Please see the <a href="https://github.com/kordk/stoch_epi_lib/releases">releases</a> to download the latest version.

### From Source Code

Once the the respository has been downloaded from GitHub,  change directory into the repository and compile the application:
<pre>
 cd stoch_epi_lib/
 go build stoch_epi_lib.go
</pre>

## Features

This distrubution was developed to be accessable, configurable, and scalable:
- The steps are modularized
- Parallel processing was implemented where available
- Support tools for constructing input data from public and standard file format
- Detailed control over parameters are provided via 
- Support tools for summarizing and reporting data (e.g., plots) 

## Usage

### Input

Data need to create the gene regulatory network (GRN):
- Transcription factor to target gene mappings.
- Binding site (i.e., methylation loci) to target gene mapping.

Data needed to perform the model parameter fitting:
- Samples with both methylation and gene expression data

Data needed to perform the gene expression estimates:
- GRN and fitted model parameters
- Methylation data for query samples

### Options

<pre>
./stoch_epi_lib.amd64.linux -h

Usage of ./stoch_epi_lib.amd64.linux:
  -Banded
        Whether or not to use a symmetric tri-diagonal matrix for UOBYQAFit or simply a diagonal. (default true)
  -BurnIn int
        Burn in for equilibrium and likelihood estimations
  -DataScale string
        How the transcript data is scaled relative to raw (Log2, Log, Log10, Linear, or none). Default none. Linear is treated as none, so model will be scaled the same way. Log scaled data is converted to unscaled for all calculation, returned equilibrium and trajectories are converted back to original scale. Parameter fitting is done to unscaled data. (default "none")
  -DistBand float
        Bandwidth for equilibrium distriubtion or likelihood estimation (default 5)
  -DistResolution int
        Resolution of Equilibrium Distriubtion Estimation. (default 50)
  -EqLength int
        Number of jumps in equilibrium distribution estimation. (default 1000)
  -EqWindow float
        Radius of equilibrium distribution to compute in units of standard deviation (i.e. number of standard deviations above & below mean to compute distribution). (default 4)
  -EstimateLength int
        Number of jumps in log-likelihood estimator realizations. (default 1000)
  -EvenRes float
        Resolution for process realizations between jumps. Set to 0 for no fill-in. (default 0.1)
  -JupPrint
        Puts a newline at the end of each line of paramfit estimating likelihoods count, so that it runs with carriage return correctly in Jupyter for debugging.
  -LoadRands string
        give an existing .json file name to the random parameters used in log-likelihood or parameter fitting calculation. (default "no")
  -MissingSite float
        Value to use for missing epigenetic data (alpha value in model) in sample. (default 1e-10)
  -MissingTranscript float
        Value to use for missing transcript in sample.
  -NumberCores int
        Extent to parallelize making quadratic in parameter fitting (number of go routines allowed). Defauts to not parallel, use -1 for max.
  -QuadPoints int
        Number of quadrature points for time averaging. (default 10)
  -SaveRands string
        give a valid .json file name to the random parameters used in log-likelihood or parameter fitting calculation. (default "no")
  -StoppingCondition int
        number of no-improvement steps before giving up (default 10)
  -datafl string
        path to and name of matched data file (default " ")
  -debug string
        Jumper 1 or Jumper 2 (default " ")
  -endtime float
        end time of realization (default 25)
  -genefl string
        file name gene parameters (default "go_in_jsons/genes.json")
  -maxsteps int
        total search steps (default 1000)
  -mode string
        mode: choose ParamFit, Fakeit, LogL, Eq, or none for sample path
  -numfake int
        number of fake samples to produce. (default 10)
  -numsamps int
        number of samples to draw for equilibrium distribution. (default 1000)
  -sitefl string
        file name of site parameters (default "go_in_jsons/sites.json")
  -svfl string
        name of save file (default "go_out_jsons/goout")
</pre>

## Demo

Please see the DEMO.md file for detailed instructions.


## Authors

James D. (Jim) Brunner: Brunner.James@mayo.edu

Jacob Kim: jk3709@cumc.columbia.edu

Kord M. Kober: Kord.Kober@ucsf.edu

## Acknowledgements
Support for this project was provided by the National Cancer Institute (CA233774). Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the National Institute of Health. This project was initially conceived as an interdisciplinary project as part of the "Short Course in Systems Biology - a foundation for interdisciplinary careers" at the Center for Complex Biological Systems at the University of California Irvine held Jan. 21 - Feb. 8, 2019 in Irvine, CA (NIH GM126365).

