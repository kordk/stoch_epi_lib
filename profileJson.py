#!/usr/local/bin/python3.8

## A program to profile a JSON file

## kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string
#import numpy as np


#### Give some usage information for this program #######################################
def usage(errorNum):
    a = """
usage: profileJson.py [hD] -j <JSON file>
"""
    print(a)
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    jsonFileName=""

    try:
        opts, args = getopt.getopt(argv, "hDj:", \
            ["help","debug","jsonFileName"])

    except getopt.GetoptError:
        usage(20)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-j", "--jsonFileName"):
            jsonFileName = arg
        elif opt in ("-D", "--debug"):
            global DEBUG
            DEBUG = 1
            print("[main] Debug set.")

    if jsonFileName == "":
        usage(202)

    htmlReport=jsonFileName+"_DataProfiling.html"

    import simplejson as json
    with open(jsonFileName,"r") as json_file:
        jf = json.load(json_file)
        
    import pandas as pd
    import pandas_profiling
    df = pd.DataFrame(jf)

    print("Loaded json pd:",jsonFileName,df.shape)
    print("DF",df.head())
    print("Profiling JSON")
    print("Saving profile as:",htmlReport)
    profile = df.profile_report(title='Pandas Profiling Report for '+jsonFileName)
    profile.to_file(htmlReport)

#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])

