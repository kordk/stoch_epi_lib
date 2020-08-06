#!/usr/bin/env python3.8

# kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string
import math
#import statistics as stats
import pandas as pd
import numpy as np
#import json
import simplejson as json
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
estGxFromEqDistPlotByGene.py - create plot of gx estimate performance by gene

usage: estGxFromEqDist.py [hD] -l <file with list of results from estGxFromEqDist.py> -p <plot directory> [-2]

 """)
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    #eqDistGxSummaryFileName="eqDistGxSummary.pickle"
    eqDistGxSummaryFileName=""
    plotDir="plots"
    log2TxfmGx=0

    try:
        opts, args = getopt.getopt(argv, "hl:p:2D", ["help",
            "eqDistGxSummary",
            "plotDir",
            "log2TxfmGx",
            "debug"])
    except getopt.GetoptError as err:
        print("Unexpected options detected:",err)
        usage(20)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-l", "--eqDistGxSummaryFileName"):
            eqDistGxSummaryFileName = arg
        if opt in ("-p", "--plotDir"):
            plotDir = arg
        if opt in ("-2", "--log2TxfmGx"):
            log2TxfmGx = 1 
        elif opt in ("-D", "--debug"):
            global DEBUG
            print("DEBUG set.")
            DEBUG = 1

    if eqDistGxSummaryFileName == "":
        usage(207)
   
    resA=[]
    resFP = open(eqDistGxSummaryFileName,"r") 
    for p in resFP:
        resA.append(p.strip('\n'))
    print("Looking for equilibrium estimations results from",len(resA),"biological samples:",eqDistGxSummaryFileName)

    print("Outputting plots to:",plotDir)
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)
    else:
        print("Directory already exists.")

    r=0
    results={}
    ## Iterate through each shuffle
    for rfile in resA:
        print("Loading file",rfile)
        df = pd.read_pickle(rfile)
        print(r,df.head())
        for g in df['GeneId']:
            #if DEBUG: print(g)
            gdf = df[ df['GeneId'] == g]
            #if DEBUG: print(gdf)
            if g not in results:
                results[g]={}
                results[g]['gxEst']=[]
                results[g]['gxObs']=[]
                results[g]['gxRes']=[]
            #print("->>",gdf['gxEst'].iloc[0])
            results[g]['gxEst'].append(float(gdf['gxEst'].iloc[0]))
            results[g]['gxObs'].append(float(gdf['gxObs'].iloc[0]))
            results[g]['gxRes'].append(float(gdf['gxRes'].iloc[0]))
            
            if DEBUG: print(r,g,results[g]['gxEst'], results[g]['gxObs'])    
        r+=1
        
        #if r == 3: break

    ## save this joined dataframe
    #if DEBUG: print("results:",results)
    results_df = pd.DataFrame(results)
    if DEBUG: print(results_df.head())

    outFileName="eqDistGxSummaryJoined.pickle"
    print("Saving combined results to:",outFileName)
    results_df.to_pickle(outFileName)

    ## plot predicted vs observed 
    ## https://stackoverflow.com/questions/25497402/adding-y-x-to-a-matplotlib-scatter-plot-if-i-havent-kept-track-of-all-the-data
    for g in results.keys():
        print("Plotting gene:",g)
        #print(g,results[g]['gxEst'])
        #lx = [math.log(x,2) for x in results[g]['gxEst']]
        #ly = [math.log(x,2) for x in results[g]['gxObs']]
        x = np.asarray(results[g]['gxEst'])
        y = np.asarray(results[g]['gxObs'])

        fig, ax = plt.subplots(nrows=1, ncols=1)

        #c = x**2 + y**2
        #c = np.asarray(results[g]['gxRes']) 
        #ax.scatter(np.log2(x),np.log2(y),
        #    label='Expected vs. Observed', cmap=plt.cm.coolwarm, c=c, s=25)
        if log2TxfmGx:
            ax.scatter(np.log2(x),np.log2(y),
                label='Expected vs. Observed', cmap=plt.cm.coolwarm, s=25)
        else:
            ax.scatter(x,y,
                label='Expected vs. Observed', cmap=plt.cm.coolwarm, s=25)
            
        ax.set_title(str(g))
        ax.set_xlabel("Predicted log2(Expression)")
        ax.set_ylabel("Observed log2(Expression)")

        #ax.set_xlim((0,15))
        #ax.set_ylim((0,15))
        #ax.axvline(x=gxEst, color='k', linestyle='--', label='Point estimate (dot*step)')
        #ax.axvline(x=gxObs, color='r', linestyle=':', label='Observed')
        #ax.legend()
        #ax.plot([0,1],[0,1], transform=ax.transAxes, color='r')

        # now plot both limits against eachother
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()])-1,  # min of both axes
            np.max([ax.get_xlim(), ax.get_ylim()])+1,  # max of both axes
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        fig.tight_layout()
        #fig.savefig(plotDir+'/'+str(a)+'_'+e['GeneID'])
        fig.savefig(plotDir+'/'+g)
        plt.close(fig)
        #break

    
#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])
