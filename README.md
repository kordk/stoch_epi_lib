# stoch_epi_lib
stoch_epi_lib: An implementation of a dynamical systems method to estimate gene expression from epigenetic data and a gene regulatory network


This software and pipeline is to estimate gene expression using the stoch_epi_lib program. The method is described in a manuscript available in pre-print on BioRxiv (https://bioxriv.org/XXX) and submitted for publication at the Pacific Symposium Biocomputing. 

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

This demonstration assumes you are running a Linux distribution with a sh/bash/ksh shell.

Download the files used in this demo from Synapse: https://www.synapse.org/#!Synapse:syn22255244/files

Copy all of the files into the same directory as the code from this repository and compile the programn. E.g.:

<pre>
git clone https://github.com/kordk/stoch_epi_lib
cd stoch_epi_lib
go build stoch_epi_lib.go
synapse get XXX
mv syn22255244/* .
</pre>

There are four main steps to follow to generate the gene expression estimates. 

These steps can be performed on multiple shuffles of the data into training and testing sets to evaluate the performance of the model. 

### 1) Data preparation

#### Omics data

Match gene expression and methylation state data are from the GTP Study (PMID 21536970). 
Transcription factor binding site to target (gene) data are from the MESA Study (PMID 29914364).

| Filename | Description | Source |
| :----- | :---------- | :---------- |
| 29914364.MESA.eCpG.txt | Annotation files containing binding site to target (Phi) data from the expression-CpG pair predictions generated using the MESA study data. | https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4842-3#Sec29 | 
| GSE72680_beta_values.txt.gz | Methylation data from the GTP Study | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72680 |

Format the CpG annotation file from the downloaded manuscript data. This is already done for the demo.
 
<pre>
cut -f1,3,4,7,10 12864_2018_4842_MOESM2_ESM.txt >29914364.MESA.eCpG.txt
</pre>

The chrom and position fields are not used.

<pre>
aggregateAssayToRegion.py - aggregate assays in a region of a gene

usage: aggregateAssayToRegion [hD] -a <annotation file> -d <gzip'd metylation data file> -s <combination score method: [ average ]>"

annotation file format:
 assayName,chrom,position,gene,label

 where the columns are
    assayName   [string] assay identifieer (e.g., cg00000029) 
    chrom       [int] chrmosome
    position    [int] nucleotide position on the chromosome
    gene        [string] label of the associated gene (e.g., TRPV1A)
    label       [string] label with which the regions will be guided (e.g., Promoter)
<pre>

Many promoter regions contain multiple methylation loci. Aggregate the epigenetic loci together and summarize as the average methylation score. Merge the GTP methylation data into the MESA regions (Promoter) and create the Phi (binding site to target) data (MESA eCpG). The region (i.e., "Promoter") is hard-coded into the script but can be changed easily.

<pre>
aggregateAssayToRegion.py \
     -s average \
     -a 29914364.MESA.eCpG.txt \
     -d GSE72680_beta_values.txt.gz
</pre>


| Filename | Description |
| :----- | :---------- |
| demo/GSE58137_Raw_119_samplesremoved.csv | Gene expression data (V3) from the GTP Study (PMID 21536970) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58137 |
| demo/ht12_v3.txt | From GEO XXX and downloaded from GEO | XX |
| demo/sym2entrex.txt | From GEO XXX and downloaded from GEO | XX |
| demo/regionAggregateData.csv | aggregated methylation scores generated by aggregateAssayToRegion.py in the previous step | XX |
| demo/GSE72680.mid2pid.txt | Mapping between the IDs of the gene expression and methylation datasets | XX |



DoRothERA (PMID 29229604)

Merge the aggregated CpG regions and matching gene expression data by individual sample:
<pre>

    
mergeMatchingData.py \
     -g demo/GSE58137_Raw_119_samplesremoved.csv \
     -d 0.05 \
     -2 \
     -p demo/ht12_v3.txt \
     -s demo/sym2entrex.txt \
     -m demo/regionAggregateData.csv \
     -a demo/GSE72680.mid2pid.txt \
     -o merged.ht12_v3.json
</pre>

#### Gene regulatory network data

Transcription factor to target gene mappings.
Epigenetic loci to gene mappings.


| Filename | Description |
| :----- | :---------- |
| demo/A_confidenceregs.csv | XXX |
| demo/regionAggregateInfo.full.tsv  | XXX |
| k.csv | Mapping of the TF to the target genes for KEGG pathway XXX |

Collect the genes in the network. For this example we are using genes in the KEGG pathway XXX:
<pre>
demo/A_confidenceregs.csv
</pre>

You can build the GRN using our tool
<pre>
python3 build_model.py \
    k.csv demo/regionAggregateInfo.full.tsv \
    >build_model.py.log 2>&1
</pre>

Move these files (which describe the GRN) into a convenient location:
<pre>
mkdir go_in_jsons
mv genes.json go_in_jsons/
mv sites.json go_in_jsons/
</pre>

#### 2) Parameter fitting on training data

| Filename | Description |
| :----- | :---------- |
| demo/merged.ht12_v3.json | Generated from mergeMatchingData.py |

For our purposes, we split the data into training (80%) and testing (20%) datasets:
<pre>
splitArrayForTrainTest.py  \
    -j demo/merged.ht12_v3.json -f 0.80
</pre>

| Filename | Description |
| :----- | :---------- |
| rand.param.json | The random parameters generated for this run. Save for use in the future. |
| train.json | Generated by splitArrayForTrainTest.py |
| go_in_jsons/sites.json | TF to binding sites maps of the GRN. Generated by build_model.py |
| go_in_jsons/genes.json | TF to target gene maps of the GRN. Generated by build_model.py |

Then, run the parameter fitting on the training data:
<pre>
mkdir go_out_jsons
stoch_epi_lib \
    -SaveRands="rand.param.json" \
    -mode=ParamFit \
    -DataScale=Log2 \
    -datafl="train.json" \
    -sitefl="go_in_jsons/sites.json" \
    -genefl="go_in_jsons/genes.json" \
    -svfl="go_out_jsons/fitted" >stoch_epi_lib.ParamFit.train.log 2>&1
</pre>

#### 3) Generate log likelihoods 

| Filename | Description |
| :----- | :---------- |
| rand.param.json | The random parameters generated for this run. Save for use in the future. |
| train.json | Generated by splitArrayForTrainTest.py |
| go_out_jsons/fitted_sites.json | Fitted parameters for the TF to binding sites maps of the GRN. Generated by stoch_epi_lib with -mode=ParamFit |
| go_out_jsons/fitted_genes.json | Fitted parameters for the TF to target gene maps of the GRN. Generated by stoch_epi_lib with -mode=ParamFit|

Run for the training data:
<pre>
stoch_epi_lib \
    -mode=LogL \
    -DataScale Log2 \
    -LoadRands="rand.param.json" \
    -datafl="train.json" \
    -sitefl="go_out_jsons/fitted_sites.json" \
    -genefl="go_out_jsons/fitted_genes.json" >stoch_epi_lib.Log2.train.log 2>&1
</pre>

| Filename | Description |
| :----- | :---------- |
| rand.param.json | The random parameters generated for this run. Save for use in the future. |
| test.json | Generated by splitArrayForTrainTest.py |
| go_out_jsons/fitted_sites.json | Fitted parameters for the TF to binding sites maps of the GRN. Generated by stoch_epi_lib with -mode=ParamFit |
| go_out_jsons/fitted_genes.json | Fitted parameters for the TF to target gene maps of the GRN. Generated by stoch_epi_lib with -mode=ParamFit|

Run for the testing data:
<pre>
stoch_epi_lib \
    -mode=LogL \
    -DataScale Log2 \
    -LoadRands="rand.param.json" \
    -datafl="test.json" \
    -sitefl="go_out_jsons/fitted_sites.json" \
    -genefl="go_out_jsons/fitted_genes.json" >stoch_epi_lib.Log2.test.log 2>&1
</pre>


#### 4) Generate gene expression equilibrium distribution

| Filename | Description |
| :----- | :---------- |
| rand.param.json | The random parameters generated for this run. Save for use in the future. |
| test.json | Generated by splitArrayForTrainTest.py |
| go_out_jsons/fitted_sites.json | Fitted parameters for the TF to binding sites maps of the GRN. Generated by stoch_epi_lib with -mode=ParamFit |
| go_out_jsons/fitted_genes.json | Fitted parameters for the TF to target gene maps of the GRN. Generated by stoch_epi_lib with -mode=ParamFit|

Create the gene expression equilibrium distributions for the testing data:
<pre>
stoch_epi_lib \
    -mode=Eq \
    -EqLength=2000 \
    -DataScale Log2 \
    -datafl="test.json" \
    -sitefl="go_out_jsons/fitted_sites.json" \
    -genefl="go_out_jsons/fitted_genes.json" \
    -svfl="go_out_jsons/Eqs" >stoch_epi_lib.Eq.test.log 2>&1
</pre>

Generate expression estimates and plot the distributions for each individual for each gene
<pre>
ALL=$PWD/plotAllEq.sh
cp /dev/null $ALL
ONE=plotOneEq.sh
N=0
PIDLIST=pid.list
ls -1 go_out_jsons/ | grep "Eqs_" | cut -d'_' -f2 | sed 's/.json//g'>$PIDLIST
wc -l $PIDLIST
P=0
for PID in `cat $PIDLIST | xargs`; do
    mkdir $PID
    pushd .
    cd $PID
    echo "PID: `pwd`"
    echo "../go_out_jsons/Eqs_${PID}.json" >eqs.list

  cat << EOF > $ONE
#!/bin/sh -x
estGxFromEqDist.py -l eqs.list -D -m ../test.json -i $PID
EOF

  chmod 755 $ONE
  echo "cd $PWD; ./$ONE >$ONE.log 2>&1" >>$ALL
  popd
  let P+=1
done
chmod 755 $ALL
./$ALL >$ALL.log 2>&1
</pre>

Generate observed vs. predicted Plots for all genes
<pre>
ONE=plotOneObsPred.sh
N=0
find $PWD -name eqDistGxSummary.pickle > eqDistGxSummary.pickle.list
wc -l eqDistGxSummary.pickle.list
cat << EOF > $ONE
#!/bin/sh -x
estGxFromEqDistPlotByGene.py -D -l eqDistGxSummary.pickle.list -p plots
EOF 
chmod 755 $ONE
./$ONE >$ONE.log 2>&1
</pre>

If you have performed multiple shuffles (e.g., n=100) each in separate subdirectories (e.g., 0/ to 99/), the results can be plotted across all of them and the mean RSME calculated:
<pre>
OUT=eqDistGxSummary.pickle.fold.list
cp /dev/null $OUT
N=0
while [ $N -lt 100 ]; do
    echo "$N,$N/eqDistGxSummaryJoined.pickle" >>$OUT
    let N+=1
done

estGxFromEqDistPlotByGeneJoinFolds.py -D -s -l eqDistGxSummary.pickle.fold.list \
    >estGxFromEqDistPlotByGeneJoinFolds.small.py.log 2>&1
    
grep "Mean RSME" estGxFromEqDistPlotByGeneJoinFolds.small.py.log  | awk '{print $4, $10}' | column -t
# AHR      2.505461
# AK3      0.599878
# ALOX5    2.846304
# BAG3     1.573793
# BAK1     4.511250
# CCM2     1.201416
</pre>

## Authors

James D. (Jim) Brunner: Brunner.James@mayo.edu

Jacob Kim: jk3709@cumc.columbia.edu

Kord M. Kober: Kord.Kober@ucsf.edu

## Acknowledgements
Support for this project was provided by the National Cancer Institute (CA233774). Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the National Institute of Health. This project was initially conceived as an interdisciplinary project as part of the "Short Course in Systems Biology - a foundation for interdisciplinary careers" at the Center for Complex Biological Systems at the University of California Irvine held Jan. 21 - Feb. 8, 2019 in Irvine, CA (NIH GM126365).

