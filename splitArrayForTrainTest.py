#!/usr/bin/env python3.8

# kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string
#import statistics as stats
import pandas as pd
#import numpy as np
#import json
import simplejson as json
from sklearn.model_selection import train_test_split

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
splitArrayForTrainTest.py - create training and testing datasets from input array of merged transcript and alpha data

usage: splitArrayForTrainTest.py [hD] -T -j <merged JSON file> -f <frequency of training data> ]

 """)
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    doSplit=1
    mergedFileName=""
    splitFreq=0.0
    doTest=0

    try:
        opts, args = getopt.getopt(argv, "hj:f:TD", ["help",
            "mergedFileName",
            "splitFreq",
            "doTest",
            "debug"])
    except getopt.GetoptError as err:
        print("Unexpected options detected:",err)
        usage(20)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-j", "--mergedFileName"):
            mergedFileName = arg
        if opt in ("-f", "--splitFreq"):
            splitFreq = float(arg)
        if opt in ("-T", "--doTest"):
            doTest = 1
        elif opt in ("-D", "--debug"):
            global DEBUG
            print("DEBUG set.")
            DEBUG = 1

    if mergedFileName == "":
        usage(207)

    if splitFreq == "":
        usage(208)

    ## load and convert it into pandas dataframe
    if DEBUG or doTest: 
        with open("./fake_matching.json") as json_file:
            fmatching = json.load(json_file)
        fmatching_df = pd.DataFrame(fmatching)
        print("Loaded matching data (fake):",fmatching_df.shape)
        print("fake:\n",fmatching_df.head())
   
        ## https://stackoverflow.com/questions/8230315/how-to-json-serialize-sets
        with open('fake.json', 'w') as handle:
                json.dump(fmatching_df, handle, iterable_as_array=True) 
        #fmatching_json = fmatching_df.to_json()
        #with open('fake.json', 'w') as handle:
        #        json.dump(fmatching_json, handle)

        with open("fake.json") as json_file:
            cmatching = json.load(json_file)
        print("Loaded matching data (fake - reload):",fmatching_df.shape)
        print("fake - reload:\n",fmatching_df.head())

    if doSplit:
        ## load and convert it into pandas dataframe
        with open(mergedFileName) as json_file:
            matching = json.load(json_file)
        matching_df = pd.DataFrame(matching) 
        if DEBUG: print("Loaded matching data (real):",matching_df.shape)
        if DEBUG: print("real:\n",matching_df.head())

        ## use sklearn to split
        test_df,train_df,=train_test_split(matching_df,test_size=splitFreq)
        train_df.reset_index(drop=True,inplace=True)
        test_df.reset_index(drop=True,inplace=True)
        if DEBUG: print("train:",train_df.shape)
        if DEBUG: print("train:",train_df.head())
        if DEBUG: print("train:",train_df.shape)
        if DEBUG: print("train:",train_df.head())
        if DEBUG: print("test:",test_df.shape)
        if DEBUG: print("test:",test_df.head())

        ## save the two datasets

        #train_json = train_df.to_json()
        #with open('train.json', 'w') as handle:
        #        json.dump(train_json, handle)

        #with open('train.json', 'w') as handle:
        #        json.dump(train_df, handle, iterable_as_array=True) 

        trainFP=open('train.json', 'w')
        trainFP.write(str(train_df.to_dict('records')).replace("'",'"'))

        if doTest:
            with open('train.json') as json_file:
                check = json.load(json_file)
            check_df = pd.DataFrame(check) 
            print("Loaded matching data (train):",check_df.shape)
            print("train:\n",check_df.head())
            print("Loaded matching data (train - reload):",check_df.shape)
            print("train - reload::\n",check_df.head())


        testFP=open('test.json', 'w')
        testFP.write(str(test_df.to_dict('records')).replace("'",'"'))

        if doTest:
            with open('test.json') as json_file:
                check = json.load(json_file)
            check_df = pd.DataFrame(check) 
            print("Loaded matching data (test):",check_df.shape)
            print("test:\n",check_df.head())
            print("Loaded matching data (test - reload):",check_df.shape)
            print("test - reload:\n",check_df.head())


    
#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])
