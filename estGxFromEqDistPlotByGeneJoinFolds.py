#!/usr/bin/env python3.8

# kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string
import math

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
estGxFromEqDistPlotByGeneJoinFolds.py - create plot of gx estimate performance by gene combining folds

usage: estGxFromEqDistPlotByGeneJoinFolds.py [hD] -l <file with list of results from estGxFromEqDistPlotByGene.py -p <plot directory>  -s

 """)
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    eqDistGxSummaryFoldsFileName=""
    plotDir="plots"
    smallRange=0

    try:
        opts, args = getopt.getopt(argv, "hl:p:sD", ["help",
            "eqDistGxSummaryFolds",
            "plotDir",
            "debug"])
    except getopt.GetoptError as err:
        print("Unexpected options detected:",err)
        usage(20)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-l", "--eqDistGxSummaryFoldsFileName"):
            eqDistGxSummaryFoldsFileName = arg
        if opt in ("-p", "--plotDir"):
            plotDir = arg
        if opt in ("-s", "--smallRange"):
            smallRange = 1
        elif opt in ("-D", "--debug"):
            global DEBUG
            print("DEBUG set.")
            DEBUG = 1

    if eqDistGxSummaryFoldsFileName == "":
        usage(207)


    #import statistics as stats
    import pandas as pd
    import numpy as np
    #import json
    import simplejson as json
    import matplotlib.pyplot as plt  
    import matplotlib.cm as cm
    # for RSME
    from sklearn.metrics import mean_squared_error
    from math import sqrt

   
    resA=[]
    resFP = open(eqDistGxSummaryFoldsFileName,"r") 
    for p in resFP:
        resA.append(p.strip('\n'))
    print("Looking for equilibrium estimations results from",len(resA),"folds:",eqDistGxSummaryFoldsFileName)

    print("Outputting plots to:",plotDir)
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)
    else:
        print("Directory already exists.")

    r=0
    foldsH = {}
    genesA = []
    for foldS in resA:
        foldA = foldS.split(',') 
        fold = foldA[0]
        rfile = foldA[1]
        print("Listing fold,",fold,"file",rfile)
        foldsH[fold]=rfile
        df = pd.read_pickle(rfile)
        if DEBUG: 
            print(r,df.head())
            print(list(df.columns))
        genesA.extend(df.columns)
        print("Found genes for,",fold,":",df.columns)
        #break
    print("Expecting data from",len(foldsH.keys()),"folds.")
    gA = np.unique(np.array(genesA))
    print("Found",gA.size,"unique genes across all folds.")
    if DEBUG: print(foldsH)

    ## collect RSME for each gene across shuffles
    dfRmse = pd.DataFrame(columns=['Gene', 'Shuffle', 'RMSE'])

    ## plot predicted vs observed 
    ## color by fold
    ## https://stackoverflow.com/questions/25497402/adding-y-x-to-a-matplotlib-scatter-plot-if-i-havent-kept-track-of-all-the-data
    from matplotlib.lines import Line2D
    ## Plot Observed (y-axis) vs. Predicted (x-axis)
    ## https://www.sciencedirect.com/science/article/pii/S0304380008002305
    print("Collecting data for plots and RMSE.")
    for g in gA.tolist():
        print("Plotting gene:",g,type(g))
        fig, ax = plt.subplots(nrows=1, ncols=1)    ## obs vs. pred
        figResPred, axResPred = plt.subplots(nrows=1, ncols=1)    ## res vs. pred

        ## colors per fold
        ## https://stackoverflow.com/questions/12236566/setting-different-color-for-each-series-in-scatter-plot-on-matplotlib
        colors = cm.rainbow(np.linspace(0, 1, len(foldsH.keys())))
        my_points=[]
        pntsResPred=[]
        for f,rfile in foldsH.items():
            if DEBUG: print("{0} => {1}".format(f, rfile))
            df = pd.read_pickle(rfile)
            print(df.head())
            print(df.columns)
            print(df.columns.map(type))
            print(df.loc[:,g])
            results = df.loc[:,g]
            print(results['gxEst'])
            
            ## add this fold to the plot
            x = np.asarray(results['gxEst'])
            y = np.asarray(results['gxObs'])
            v = np.asarray(results['gxRes'])

            c = colors[int(f)]
            ax.scatter(np.log2(x),np.log2(y), c=c, s=2)
            ## Residual = Observed – Predicted
            axResPred.scatter( 
                np.log2(x),
                np.log2(y)-np.log2(x),
                c=c, s=2)
            
            ## add to the legend
            my_points.append(Line2D([0], [0], marker='.', color=c, lw=2))
            pntsResPred.append(Line2D([0], [0], marker='.', color=c, lw=2))

            ## calculate RSME
            #rmse = sqrt(mean_squared_error(y_actual, y_predicted))
            rmse = sqrt(mean_squared_error(y, x))
            dfRmse = dfRmse.append({'Gene': g, 'Shuffle': f, 'RMSE': rmse}, ignore_index=True)
            print("RSME",g,"shuffle",f,rmse)

            #break

        #print(dfRmse)
        dfr = dfRmse.loc[dfRmse['Gene']==g]
        #print(dfr)
        print(dfr['RMSE'])
        print("Mean RSME for ",g,dfr['RMSE'].mean)

        print(my_points)
        ax.set_title(str(g))
        ax.set_xlabel("Predicted log2(Expression)")
        ax.set_ylabel("Observed log2(Expression)")
        if smallRange:
            ## use data to set limits
            lims = [
                np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
                np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
            ]
            figNameBase=g+'-small-'
        else:
            ## hard set limits
            lims = [
                -5.0,  # min of both axes
                15.0,  # max of both axes
            ]
            figNameBase=g

        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0, linestyle="--")
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        fig.tight_layout()
        #ax.legend(frameon=False)
        #ax.legend(my_points, foldsH.items(), frameon=False)
        #ax.plot([0,1],[0,1], transform=ax.transAxes, color='r')
        fig.savefig(plotDir+'/'+figNameBase)
        plt.close(fig)


        print(pntsResPred)
        axResPred.set_title(str(g))
        axResPred.set_xlabel("Predicted log2(Expression)")
        axResPred.set_ylabel("Residuals = log2(Obs) - Log2(Pred)")
        ## hard set limits
        if smallRange:
            ## use data to set limits
            lims = [
                np.min([axResPred.get_xlim(), axResPred.get_ylim()]),  # min of both axes
                np.max([axResPred.get_xlim(), axResPred.get_ylim()]),  # max of both axes
            ]
        else:
            lims = [
                -5.0,  # min of both axes
                15,  # max of both axes
            ]

        #axResPred.plot(lims, lims, 'k-', alpha=0.75, zorder=0, linestyle=":")
        axResPred.axhline(y=0.0, linestyle=":")
        axResPred.set_aspect('equal')
        axResPred.set_xlim(lims)
        axResPred.set_ylim(lims)
        figResPred.tight_layout()
        #ax.legend(frameon=False)
        #ax.legend(my_points, foldsH.items(), frameon=False)
        #ax.plot([0,1],[0,1], transform=ax.transAxes, color='r')
        figResPred.savefig(plotDir+'/'+figNameBase+'rp')
        plt.close(figResPred)
        

        #break

    
#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])
