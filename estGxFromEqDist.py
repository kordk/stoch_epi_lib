#!/usr/bin/env python3.8

# kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string
#import statistics as stats
import pandas as pd
import pandas_profiling
import numpy as np
#import json
import simplejson as json
import requests
import matplotlib.pyplot as plt  

DEBUG=0

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

#### Give some usage information for this script #######################################
def usage(errorNum):
    print("""
estGxFromEqDist.py - create gene expression estimations from the equilibrium distributions, plot, and generate pickle for shared plots

usage: estGxFromEqDist.py [hDP] -l <file with list of Eq filenames for the patient> -p <plot directory> -m <merged json> -s <sample index in merged json>

 """)
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    eqListFileName=""
    plotDir="plots"
    mergedFileName=""
    sampleId=""
    profileJson=0

    try:
        opts, args = getopt.getopt(argv, "hj:l:p:m:i:PD", ["help",
            "eqListFileName",
            "plotDir",
            "mergedFileName",
            "sampleId",
            "profileJson",
            "debug"])
    except getopt.GetoptError as err:
        print("Unexpected options detected:",err)
        usage(20)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-l", "--eqListFileName"):
            eqListFileName = arg
        if opt in ("-p", "--plotDir"):
            plotDir = arg
        if opt in ("-m", "--mergedFileName"):
            mergedFileName = arg
        if opt in ("-i", "--sampleId"):
            sampleId = str(arg)
        if opt in ("-P", "--profileJson"):
            profileJson = 1
        elif opt in ("-D", "--debug"):
            global DEBUG
            print("DEBUG set.")
            DEBUG = 1

    if eqListFileName == "":
        usage(207)

    if mergedFileName == "":
        usage(209)

    if DEBUG:
        print("sampleId:",type(sampleId),sampleId)
    if sampleId == "":
        usage(211)

    eqA=[]
    eqFP = open(eqListFileName,"r") 
    for p in eqFP:
        eqA.append(p.strip('\n'))
    print("Looking for equilibrium estimations for",len(eqA),"biological samples.")

    print("Outputting plots to:",plotDir)
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)
    else:
        print("Directory already exists.")

    ## Load the matching data to obtain the observed Gx value
    with open(mergedFileName,"r") as json_file:
        matching = json.load(json_file)
    #if DEBUG: print("Loaded matching data (real) raw:\n", json.dumps(matching, indent=2) )
   
    ## NOTE: if this head() and tail() fail, try upgrading pandas and python before pulling your hair out
    matching_df = pd.DataFrame(matching) 
    if DEBUG: print("Loaded matching data (real) pd:",mergedFileName,matching_df.shape)
    if DEBUG: print("real:\n",matching_df.head())
    mdf=matching_df.set_index("Label")

    #if DEBUG: print("real:\n",matching_df.tail())
    if profileJson:
        print("Profiling JSON")
        profile = mdf.profile_report(title='Pandas Profiling Report')
        profile.to_file("testJson_dataProfiling.html")
    
    #testSymbol="CXCR5"
    testSymbol="C17orf75"
    print(mdf.head())
    gxmdf=mdf.loc[sampleId]
    print("gxmdf:\n",gxmdf)
    print("gxmdf['Transcript']:\n", dict(list(gxmdf['Transcript'].items())[0:3]))
    #print(gxmdf['Transcript'])
    #print(gxmdf['Transcript'].to_dict()[testSymbol])
    #print(gxmdf['Transcript'].to_dict()[testSymbol])
    gx=gxmdf['Transcript'][testSymbol]
    if DEBUG: print("sample",sampleId,"gene",testSymbol,":",gx,2**gx)
   
    #print("***WARNING. Exiting prematurely for debug.")
    #sys.exit(8)

    results = {}
    for estfile in eqA:
        if DEBUG: print("Opening json file:",estfile)
        with open(estfile) as f:
            est = json.load(f)
            if profileJson:
                print("Profiling JSON")
                profile = est.profile_report(title='Pandas Profiling Report')
                profile.to_file("estFileJson_dataProfiling.html")
                
            g=0
            results["GeneId"]=[]
            results["gxEst"]=[]
            results["gxEstSd"]=[]
            results["gxObs"]=[]
            results["gxRes"]=[]
    
            ## iterate through genes
            for e in est:
                
                gxEst=-1.0
                gxObs=-1.0

                ## collect the estimates
                if DEBUG: print(sampleId,g,e['GeneID'],len(e['XDomain']), e['MeanScaled'], e['StdScaled'])
                gxEst   = e['MeanScaled']
                gxEstSd = e['StdScaled']

                if DEBUG: print(sampleId,g,e['GeneID'],len(e['XDomain']), len(e['DistEst']),"-Est->", gxEst)
                if DEBUG: print(sampleId,g,e['GeneID'],len(e['XDomain']), len(e['DistEst']),"-Esd->", gxEstSd)

                ## collect the observed
                try:
                    gxObs=mdf.loc[sampleId]['Transcript'][e['GeneID']]
                except KeyError:
                    print("Unable to find observed Gx data for:",g,"Skipping.")
                    continue
                if DEBUG: print(sampleId,g,e['GeneID'],len(e['XDomain']), len(e['DistEst']),"-Obs->", gxObs)
                gxRes = gxObs-gxEst
                if DEBUG: print(sampleId,g,e['GeneID'],len(e['XDomain']), len(e['DistEst']),"-Res->", gxRes)

                results["GeneId"].append(e['GeneID'])
                results["gxObs"].append(gxObs)
                results["gxEst"].append(gxEst)
                results["gxEstSd"].append(gxEstSd)
                results["gxRes"].append(gxRes)

                print("Mean expression values:","gxObs",np.mean(results["gxObs"]),"gxEst", np.mean(results["gxEst"]))

                ## plot each one
                fig, ax = plt.subplots(nrows=1, ncols=1)
                ax.plot(e['XDomain'],e['DistEst'], label='Predicted')
                ax.set_title("Sample: "+sampleId+" Gene: "+str(e['GeneID']))
                ax.set_ylabel("Density")
                ax.set_xlabel("Expression")
                ax.axvline(x=gxEst, color='k', linestyle='--', label='MeanScaled')
                ax.axvline(x=gxObs, color='r', linestyle=':', label='Observed')
                ax.legend()
                fig.tight_layout()
                #fig.savefig(plotDir+'/'+str(a)+'_'+e['GeneID'])
                fig.savefig(plotDir+'/'+e['GeneID'])
                plt.close(fig)
                #break
                g+=1
        break
    if DEBUG: print("results:",results)
    results_df = pd.DataFrame(results) 
    if DEBUG: print(results_df.head())

    outFileName="eqDistGxSummary.pickle"
    print("Saving results to",outFileName)
    results_df.to_pickle(outFileName)


    
#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])
