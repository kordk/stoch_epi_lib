#!/usr/bin/env python3

# kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string,gzip
import statistics as stats
import pandas as pd

DEBUG=0

## Only look at Promoter and TSS sites as TFBS for now
labelsToEvalH={
    "TSS":1, 
    "Promoter":1
}

def get_size(obj, seen=None):
    ### https://goshippo.com/blog/measure-real-size-any-python-object/
    "Recursively finds size of objects"
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def loadAssayData(afn):
    "Read assay data in manually"
    lineCount=0
    headerA=[]
    probesAnnotH={}
    probesMatchH={}
    with gzip.open(afn,'r') as fin:        
        for line in fin:        
            lineA=line.strip("\n").split("\t")
            lineCount+=1
            ## read in the header
            if len(headerA) == 0:
                headerA=lineA
                continue

            ## only include probes for which we have annotation data
            if lineA[0] not in probesAnnotH:
                continue

            probe=lineA.pop(0)
            if DEBUG: print( "[main] Matched probe:",probe)
            probesMatchH[probe]=1

            ## store the data for this probe
            arrayDataH[probe]=lineA
            
        print( "Processed",lineCount,"lines.")
        print( "Matched",len(probesMatchH),"probes.")
        s=get_size(arrayDataH)
        print("Cached array data for eCpG probes (",s/1024/1024,"Mbytes)." )
    return assayFileName,probesMatchH


def readAnnotatioFileToDict(annFileName):
    "process the assay file for regions"
    annA = []
    annFP = open(annFileName,"r")
    annH={}
    probesH={}

    ## expected format:
    # CpG.probe   exp.probe.chrm  exp.probe.start  annot.gene  status
    # cg00000236  2               231103633        SP140       TRANS
    # cg00000289  15              41809627         RPAP1       TRANS
    # cg00000363  4               146055861        OTUD4       TRANS

    for line in annFP:
        annA = line.strip('\n').split('\t')
        probe = annA[0]
        chrom = annA[1]
        start = annA[2]
        gene = annA[3]
        label = annA[4]

        #if not annH.has_key(probe):
        #    annH[probe]={}
        #    probesH[probe]=1
        #if not annH[probe].has_key(gene):
        #    annH[probe][gene]={}
        #if not annH[probe][gene].has_key(label):
        #    annH[probe][gene][label]={}
        #if not annH[probe][gene][label].has_key("pos"):
        #    annH[probe][gene][label]["pos"]=[]
        #annH[probe][gene][label]["pos"].append(str(chrom+":"+start))
        #if DEBUG: print "[readAnnotatioFileToDict]",probe,annH[probe])
        
        if gene not in annH:
            annH[gene]={}
        if label not in annH[gene]:
            annH[gene][label]={}
            probesH[probe]=1
        if probe not in annH[gene][label]:
            annH[gene][label][probe]={}
        if "pos" not in annH[gene][label][probe]:
            annH[gene][label][probe]["pos"]=[]
        annH[gene][label][probe]["pos"].append(str(chrom+":"+start))
        if DEBUG: print("[readAnnotatioFileToDict]",gene,annH[gene])
    return annH,probesH


#### Give some usage information for this script #######################################
def usage(errorNum):
    print("""
aggregateAssayToRegion.py - aggregate assays in a region of a gene"

usage: aggregateAssayToRegion [hD] -a <annotation file> -d <gzip'd metylation data file> -s <combination score method: [ average ]>"

annotation file format:
 assayName,chrom,position,gene,label

 where the columns are
    assayName   [string] assay identifieer (e.g., cg00000029) 
    chrom       [int] chrmosome
    position    [int] nucleotide position on the chromosome
    gene        [string] label of the associated gene (e.g., TRPV1A)
    label       [string] label with which the regions will be guided (e.g., Promoter)
 """)
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    annotationFileName=""
    assayFileName=""
    combineScoreMethod=""

    try:
        opts, args = getopt.getopt(argv, "ha:d:s:D", ["help","annotationFileName", \
            "assayFileName","combineScoreMethod","debug"])
    except getopt.GetoptError as err:
        print("Unexpected options detected:",err)
        usage(20)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-a", "--annotationFileName"):
            annotationFileName = arg
        if opt in ("-d", "--assayFileName"):
            assayFileName = arg
        if opt in ("-s", "--combineScoreMethod"):
            combineScoreMethod = arg
        elif opt in ("-D", "--debug"):
            global DEBUG
            print("DEBUG set.")
            DEBUG = 1
   
    if annotationFileName == "":
        usage(202)

    if assayFileName == "":
        usage(203)
        
    if combineScoreMethod == "":
        usage(204)
    
    ## read in the  annotation file and collect the probe names for inclusion
    (assaysForGeneRegionH,probesAnnotH)=readAnnotatioFileToDict(annotationFileName)
    print("Found",len(probesAnnotH),"eCpG probes in annotation.")

    ## process each gene
    print("Processing labels",labelsToEvalH,"for",len(assaysForGeneRegionH),"genes")
    n=0
    probeMergeH={}
    for gene,dataH in assaysForGeneRegionH.items():
        if DEBUG: print("------------")
        if DEBUG: print(gene, dataH.keys())
        for label in labelsToEvalH.keys():
            if label not in dataH:
                if DEBUG: print(label,"Not found.")
                continue
            if DEBUG: print("Merging",label,dataH[label])
            ## check distance (sliding window?)
            if gene not in probeMergeH:
                probeMergeH[gene]={}
            ## name it
            rname=label+"_"+"-".join(dataH[label].keys())
            probeMergeH[gene][rname]=[]
            probeMergeH[gene][rname].append(dataH[label].keys())
            if DEBUG: print("Merged",rname,label,gene,probeMergeH[gene])
        n+=1
        #if n == 25:
        #    print("***WARNING. Pre-maturely terminating processing of annotation data. Are you debugging?")
        #    break
    print("Merged probes into",len(probeMergeH),"regions.")
            
        
    ## cache the assay data from file
    print("Caching array data & filtering out unused probes.")
    reader = pd.read_csv(assayFileName, sep='\t', chunksize=1024*50)
    c=0
    for chunk in reader:
        #print chunk
        chunk.rename(columns={ chunk.columns[0]: "probe" }, inplace=True)
        #print chunk
        #print c,chunk.shape
        #chunk.drop(['cg00000289'], inplace=True)
        #keeps=['cg00000289','foo','cg00002033']
        #chunk=chunk.loc[chunk['probe'].isin(keeps)]
        #chunk = chunk.loc[keeps,]
        #chunk = chunk.loc[probesAnnotH.keys(),]
        chunk = chunk.loc[chunk['probe'].isin(probesAnnotH.keys())]
        #print chunk
        #print "chunk",c,chunk.shape
        if c == 0:
            print("Creating arrayData from first chunk")
            arrayData = chunk
        else:
            print("Adding chunk to arrayData",c,chunk.shape)
            arrayData = arrayData.append(chunk)
        #print "arrayData",arrayData
        c+=1
        #if c == 5:
        #    print("***WARNING. Pre-maturely terminating processing of assay data. Are you debugging?")
        #    break
    #print "arrayData",arrayData
    print("arrayData",arrayData.shape)
    s=get_size(arrayData)
    print("Cached array data for eCpG probes (",s/1024,"Kbytes)."  )
    print("Removing unneeded columns.")
    arrayData = arrayData.loc[:,~arrayData.columns.str.startswith('Detection')]
    print("arrayData",arrayData.shape)
    print("Cached array data for eCpG probes (",s/1024,"Kbytes)."  )
    print("Reseting index.")
    arrayData.set_index('probe', inplace=True)
    #print arrayData.columns

    ## Merge regions for the array data
    regionData = pd.DataFrame(columns=arrayData.columns)
    regionAnnotation = pd.DataFrame(columns={'region','annotProbes','annotCount','arrayProbes','arrayCount'})
    for gene,regionH in probeMergeH.items():
        for region,probesA in regionH.items():
            #print("Collecting array data for",gene,region,probesA)
            if len(probesA) > 1:
                print("***NOTICE:",gene,region,"has mulitple probe sets:",probesA)
                print("        Using only the first.")
            regionAnnotation.at[gene,'region'] = region
            regionAnnotation.at[gene, 'annotProbes'] = ','.join(probesA[0])
            regionAnnotation.at[gene, 'annotCount'] = len(probesA[0])
            #print(regionAnnotation)
            regionData.loc[region] = [None] * len(regionData.columns)
            for pid in arrayData.columns:
                scoresA=[]
                arrayProbes=[]
                for probe in probesA[0]:
                    if not probe in arrayData.index:
                        #print "Probe not found:",probe
                        continue
                    #print "Score for",probe,":",arrayData.get_value(probe,pid)
                    scoresA.append(float(arrayData.at[probe,pid]))
                    arrayProbes.append(probe)
                if DEBUG: print(pid,probesA[0],scoresA)
                try: 
                    score=scoresA[0]
                except:
                    #print("***NOTICE: score not found for",probe,pid,scoresA)
                    score=None
                if len(scoresA) > 1:
                    score=stats.mean(scoresA)
                regionData.at[region, pid]=score
                regionAnnotation.at[gene, 'arrayProbes'] = ','.join(arrayProbes)
                regionAnnotation.at[gene, 'arrayCount'] = len(arrayProbes)
            #if len(scoresA) > 1:
                #print("Region",gene,region,"has >1 probe in the array data:",probesA[0])
            if len(scoresA) == 0:
                print("Region",gene,region,"has no probes in the array data:",probesA[0])
            #print regionData
            #print(regionAnnotation)
        #print("Completed region(s) for",gene,regionData.shape)
        #break
    print(regionData.head(5))
    print(regionAnnotation.head(5))
    

    ## Save the region information
    outFileName="regionAggregateInfo.tsv.gz"
    print("Saving aggregated information to:",outFileName)
    regionAnnotation.to_csv(outFileName, sep='\t', float_format='%.4f', index_label="gene", compression='gzip')

   
    ## export the data
    outFileName="regionAggregateData.csv.gz"
    print("Saving aggregated data to:",outFileName)
    regionData.to_csv(outFileName, sep=',', float_format='%.4f', index_label="region", compression='gzip')
    
    
#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])
