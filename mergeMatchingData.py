#!/usr/bin/env python3.8

# kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string,gzip
import statistics as stats
import math
import pandas as pd
import numpy as np

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

def readProbeAnnotToDict(fn):
    "process the ENTREZ annotation file for probes"
    annFP = open(fn,"r")
    geneH={}

    ## expected format
    #ILMN_Gene  RefSeq_ID       Entrez_Gene_ID  Probe_Id
    #LOC23117   XM_933824.1     23117           ILMN_1725881
    #FCGR2B     XM_938851.1     2213            ILMN_1804174
    #TRIM44     NM_017583.3     54765           ILMN_1796063
    #LOC653895  XM_936379.1     653895          ILMN_1811966

    print("[readProbeAnnotToDict] Mapping [Probe ID] -> ENTREZ Gene ID.")
    for line in annFP:
        annA = line.strip('\n').split('\t')
        if len(annA) == 4:
            geneH[str(annA[3]).strip()] = str(annA[2]).strip() # Entrez_Gene_ID
            #geneH[annA[3]] = annA[1] # RefSeq_ID
            #geneH[annA[3]] = annA[0] # ILMN_Gene

    #print("[readProbeAnnotToDict] Probe ID: ILMN_1762337 -> ENTREZ:", geneH['ILMN_1762337'])
    print("[readProbeAnnotToDict] Loaded",len(geneH.keys()),"Probe Ids.")
    print("[readProbeAnnotToDict] geneH:", dict(list(geneH.items())[:3]))
    #print("[readProbeAnnotToDict] Probe ID: ILMN_1762337 -> 920?", geneH['ILMN_1727284'])
    return geneH

def readSymbolAnnotToDict(fn):
    "process the symbol annotation file for ENTREZ IDs" 
    annFP = open(fn,"r")
    symH={}

    ## expected format
    #hgnc_id     symbol      entrez_id
    #HGNC:5      A1BG        1
    #HGNC:37133  A1BG_AS1    503538
    #HGNC:24086  A1CF        29974

    print("[readSymbolAnnotToDict] Mapping [ENTREZ Gene ID] -> Symbol.")
    for line in annFP:
        annA = line.strip('\n').split('\t')
        symH[str(annA[2]).strip()] = annA[1] # symbol

    
    print("[readSymbolAnnotToDict] keys:", list(symH.keys())[0:4])
    print("[readSymbolAnnotToDict] ENTREZ ID: Symbol", dict(list(symH.items())[:5]))
    print("[readSymbolAnnotToDict] ENTREZ ID: 346389 -> Symbol:", symH['346389'])
    print("[readSymbolAnnotToDict] Loaded,",len(symH.keys()),"ENTREZ IDs.")
    print("[readSymbolAnnotToDict] ENTREZ ID: 1 -> A1BG?", symH['1'])
    print("[readSymbolAnnotToDict] ENTREZ ID: 920 -> CD4?", symH['920'])

    return symH

def updateGeneNames(pfn,sfn,df):
    probeH = readProbeAnnotToDict(pfn)
    symbolH = readSymbolAnnotToDict(sfn)
    print("[updateGeneNames] Loaded,",len(symbolH.keys()),"ENTREZ IDs.")
    print("[updateGeneNames] ENTREZ ID: 1 -> A1BG?", symbolH['1'])

    print("[updateGeneNames] Mapping [ENTREZ Gene ID] -(Probes)-> Symbol.")
    p2sH={}
    ## get the ENTREZ ID for the probe
    m=0
    for p,e in probeH.items():
        #print("[updateGeneNames] Looking for probe",p,"ENTREZ",e, str(e), str(e).strip(),symbolH[str(e).strip()])
        if not e in symbolH:
            m+=1
            continue
        p2sH[p]=symbolH[e] ## set the HUGO symbol for the ENTREZ ID
        

    print("[updateGeneNames] Annotated",len(p2sH),"probes. Failed to annotate",m,"probes.")
    print("[updateGeneNames] p2sH:", dict(list(p2sH.items())[:3]))

    print("[updateGeneNames] Quick verification of lookup hash. If this fails, you may need to update the code to comment out the appropriate test.")

    # use with Illumina array
    print("[updateGeneNames] Probe ID: ILMN_1762337 -> Symbol:", p2sH['ILMN_1762337'])
    print("[updateGeneNames] Probe ID: ILMN_1727284 -> CD4?:", p2sH['ILMN_1727284'])

    # use with RNAseq
    #print("[updateGeneNames] Probe ID: 1 -> Symbol:", p2sH['1'])
    #print("[updateGeneNames] Probe ID: 7398 -> USP1?:", p2sH['7398'])

    print("[updateGeneNames] Annotating probes as HUGO symbols.")
    df.replace({'Probe_ID':p2sH}, inplace=True)
    #print("[updateGeneNames]",df.iloc[:15, :5] )
    print("[updateGeneNames]",df.head )
    return df

def readRegionAggregateInfoToDict(fn):
    "process the region annotation file for regions"
    annFP = open(fn,"r")
    raiH={}

    ## expected format
    # gene     arrayCount  annotProbes  region               arrayProbes  annotCount
    # CRYL1    1           cg23775119   TSS_cg23775119       cg23775119   1
    # METTL22  1           cg03698956   TSS_cg03698956       cg03698956   1
    # GBF1     1           cg06839377   Promoter_cg06839377  cg06839377   1
    # CAMK1    1           cg13896608   Promoter_cg13896608  cg13896608   1

    for line in annFP:
        annA = line.strip('\n').split('\t')
        raiH[region]= annA[3]
   
    print("[readRegionAggregateInfoToDict] Loaded,",len(raiH.keys()),"regions.")

    print("[readRegionAggregateInfoToDict] Found region Promoter_cg13512987-cg05044173?", \
        raiH["Promoter_cg13512987-cg05044173"])
    return raiH

def loadGxData(fn,detPval):
    "load gene expression data"
    reader = pd.read_csv(fn, sep='\t', chunksize=1024*50)
    print("[loadGxData] Loading Gx data.")

    #co = 0.05
    co = float(detPval)
    print("[loadGxData] Attempting to filter out probes with Detection pvalue <",co)

    c=0
    dropArrayNames = []
    for chunk in reader:
        print("[loadGxData] chunk",c,"shape:", chunk.shape)

        ## Check Detection Value
        ## collect the column names
        ## axis=1: evaluate across the rows (i.e., row minimum)
        dcols = list(chunk.filter(regex='^Detection pvalue_', axis=1))
        print(dcols)
        
        ## Filter out probes which were not expressed (detection < cutoff 'co') in any sample
        ## record and output?
        nfilt = chunk[dcols].min(axis=1) < co
        print("[loadGxData] Chunk",c,"detection filter (<",co,"):\n", nfilt.value_counts())
        if DEBUG: print(nfilt.head())

        ## Test for missing values
        ## all() checks wether all values in a list interpret to True, meaning, if at least one of them is False, it will return False
        print(all(nfilt))
        if all(nfilt):
            print("[loadGxData] Chunk",c,"detection filter columns found.")
            ## only drop if column is present for all samples
            cDrop = chunk.loc[chunk[dcols].min(axis=1) > co]
            #print(cDrop.columns.to_numpy().tolist())
            #dropArrayNames = cDrop[['Probe_ID']].values
            dropArrayNames = cDrop[['Probe_ID']].values.tolist()
            chunk = chunk.loc[chunk[dcols].min(axis=1) < co]
            print("[loadGxData]", dropArrayNames[1:10])
        else:
            print("[loadGxData] Chunk",c,"detection filter columns not found. Not performing detection filter.", chunk.shape)

        if c == 0:
            print("[loadGxData] Creating array df from first chunk")
            arrayData = chunk
        else:
            print("[loadGxData] Adding chunk to array df",c,chunk.shape)
            arrayData = arrayData.append(chunk) 
        c+=1
        #if c == 5:
        #    print("***\n***\n")
        #    print("***WARNING. Pre-maturely terminating processing of assay data. Are you debugging?")
        #    print("***\n***\n")
        #    break
    #print("[loadGxData] info",arrayData.info(memory_usage='deep'))

    fname="gxAssayDetPvalDrop.txt"
    print("[loadGxData] Writing excluded Gx probes failing Detection p-value to file:",fname)
    with open(fname, 'w') as f:
        for item in dropArrayNames:
            f.write("%s\n" % item)
    
    print("[loadGxData] loaded shape:", arrayData.shape)
    print("[loadGxData] Removing unneeded columns.")
    arrayData = arrayData.loc[:,~arrayData.columns.str.startswith('Detection')]
    #print("[loadGxData] info",arrayData.info(memory_usage='deep'))
    print("[loadGxData] new shape:", arrayData.shape)
    
    ## Filter out probes with missing data (i.e., NaN)
    print("[loadGxData] Removing probes with missing data.")
    arrayData_nona = arrayData.dropna()
    print("[loadGxData] final shape:", arrayData_nona.shape)
        
    print("[loadGxData] Fixing the column names.")
    if DEBUG: print("[loadGxData] Before:", arrayData_nona.columns.values)
    #print(type(arrayData.columns.values))
    new = np.array([re.sub('Average signal_','',a) for a in arrayData_nona.columns.values])
    arrayData_nona.columns = new.tolist()
    if DEBUG: print("[loadGxData] After:", arrayData_nona.columns.values)
    print("[loadGxData] arrayData_nona:\n", arrayData_nona.head(n=5))
    print("[loadGxData] arrayData_nona:\n", arrayData_nona.dtypes.head(n=5))

    print("[loadGxData] Forcing string cast on Probe_ID column.")
    arrayData_nona['Probe_ID'] = arrayData_nona['Probe_ID'].astype(str)
    print("[loadGxData] arrayData_nona:\n", arrayData_nona.dtypes.head(n=5))

    print("[loadGxData] [QC] ILMN_1727284/CD4:", arrayData_nona.loc[arrayData_nona['Probe_ID']=="ILMN_1727284"], flush=True)
    print("[loadGxData] [QC] 1/A1BG:", arrayData_nona.loc[arrayData_nona['Probe_ID']=="1"], flush=True)

    return arrayData_nona

def readSampleIdAnnotToDict(fn):
    "process the sample annotation file for methylation arrays" 
    annFP = open(fn,"r")
    a2sH={}

    ## expected format
    #!Sample_description  6164647041_R01C01  6164647041_R01C02  6164647041_R02C01
    #!Sample_title        6088               6094               6163

    print("[readSampleIdAnnotToDict] Collecting Sample ID for Methlylation Array ID.")

    df = pd.read_csv(fn, sep='\t', dtype=str)
    print("[readSampleIdAnnotToDict] keys:", list(df.columns.values)[0:4])

    ## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1868036
    ## Expect: 6164647041_R01C01 -> 6088
    #print("[readSampleIdAnnotToDict] 6164647041_R01C01 -> 6088?", df['6164647041_R01C01'].item())

    a=0
    for an in df.columns.values:
       a2sH[str(an)]=str(df[an].item())
       if DEBUG: print("[readSampleIdAnnotToDict] Loaded example:", an, df[an].item())
       #if a > 5:
       #     print("***\n***\n")
       #     print("***WARNING. Pre-maturely terminating processing of array sample ID data. Are you debugging?")
       #     print("***\n***\n")
       #     break
       a+=1

    if DEBUG: print("[readSampleIdAnnotToDict] a2sH:",a2sH)
    return a2sH
  
def loadMtData(mfn,afn):
    "load methylation data"

    reader = pd.read_csv(mfn, sep=',', chunksize=1024*50)
    print("[loadMtData] Loading Mt data.")

    c=0
    for chunk in reader:
        print("[loadMtData] chunk",c,"shape:", chunk.shape)
        if c == 0:
            print("[loadMtData] Creating array df from first chunk")
            arrayData = chunk
        else:
            print("[loadMtData] Adding chunk to array df",c,chunk.shape)
            arrayData = arrayData.append(chunk) 

        c+=1
        #if c == 5:
        #    print("***\n***\n")
        #    print("***WARNING. Pre-maturely terminating processing of assay data. Are you debugging?")
        #    print("***\n***\n")
        #    break

    print("[loadMtData] loaded shape:", arrayData.shape)
    print("[loadMtData] Annotating arrays with patient IDs.")
    if DEBUG: print("[loadMtData] Before:", arrayData.head)
    aidH = readSampleIdAnnotToDict(afn)
    if DEBUG: print("[loadMtData] aidH['6164647041_R01C01']",aidH['6164647041_R01C01'])
    arrayData.rename(columns=aidH, inplace=True)

    if DEBUG:
        print("[loadMtData] Promoter_cg13512987-cg05044173 (before):", \
            arrayData.loc[arrayData['region']=="Promoter_cg13512987-cg05044173"])
        print("[loadMtData] Promoter_cg21689902-cg03673190 (before):", \
            arrayData.loc[arrayData['region']=="Promoter_cg21689902-cg03673190"])
    
    ## Filter out probes with missing data (i.e., NaN)
    print("[loadMtData] Removing probes with missing data.") 

    nDrop=arrayData[arrayData.isna().any(axis=1)]
    print("[loadMtData] nDrop:", nDrop.shape)

    dropArrayNames = nDrop[['region']].values.tolist()
    print("[loadMtData]", dropArrayNames[1:10])

    fname="mtAssayNaDrop.txt"
    print("[loadMxData] Writing excluded Mt probes with missing data to file:",fname)
    with open(fname, 'w') as f:
        for item in dropArrayNames:
            f.write("%s\n" % item)

    arrayData_nona = arrayData.dropna(how='any', axis=0)
    print("[loadMtData] final shape:", arrayData_nona.shape)

    if DEBUG:
        print("[loadMtData] Promoter_cg13512987-cg05044173 (after):", \
            arrayData_nona.loc[arrayData_nona['region']=="Promoter_cg13512987-cg05044173"])
        print("[loadMtData] Promoter_cg21689902-cg03673190 (after):", \
            arrayData_nona.loc[arrayData_nona['region']=="Promoter_cg21689902-cg03673190"])
    
    print("[loadMtData]",arrayData_nona.head, flush=True)
    return arrayData_nona

#### Give some usage information for this script #######################################
def usage(errorNum):
    print("""
mergeMatchingData.py - merge epigenetic and gene expression data matched to an individual

usage: mergeMatchingData.py [hD2m:g:d:p:s:]
    -2 log2 transform gene expression data
    -m gzip'd methylation data file
    -a methylation sample name annotation
    -g gene expression data
    -d gene expression detection p-value cutoff
    -p Probe ID to ENTREZ annotation file name 
    -s ENTREZ to Symbol annotation file name

 """)
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    geneExpFileName=""
    probeAnnotationName=""
    symbolAnnotationName=""
    methylFileName=""
    sampleMethylAnnotationName=""
    mergedFileName=""
    log2TxfmGx=0
    gxDetectionPvalue=0.05

    try:
        opts, args = getopt.getopt(argv, "hg:p:s:m:a:o:2d:D", ["help",
            "geneExpFileName",\
            "probeAnnotationName", \
            "symbolAnnotationName", \
            "methylFileName",\
            "sampleMethylAnnotationName",\
            "mergedFileName",\
            "log2TxfmGx",\
            "gxDetectionPvalue",\
            "debug"])
    except getopt.GetoptError as err:
        print("Unexpected options detected:",err)
        usage(20)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-g", "--geneExpFileName"):
            geneExpFileName = arg
        if opt in ("-p", "--probeAnnotationName"):
            probeAnnotationName = arg
        if opt in ("-s", "--symbolAnnotationName"):
            symbolAnnotationName = arg
        if opt in ("-m", "--methylFileName"):
            methylFileName = arg
        if opt in ("-a", "--sampleMethylAnnotationName"):
            sampleMethylAnnotationName = arg
        if opt in ("-o", "--mergedFileName"):
            mergedFileName = arg
        if opt in ("-2", "--log2TxfmGx"):
            log2TxfmGx = 1
        if opt in ("-d", "--gxDetectionPvalue"):
            gxDetectionPvalue = arg
        elif opt in ("-D", "--debug"):
            global DEBUG
            print("DEBUG set.")
            DEBUG = 1
   
    if geneExpFileName == "":
        usage(202)

    if probeAnnotationName == "":
        usage(203)

    if symbolAnnotationName == "":
        usage(204)
        
    if methylFileName == "":
        usage(205)

    if sampleMethylAnnotationName == "":
        usage(206)

    if mergedFileName == "":
        usage(207)

    if log2TxfmGx:
        print("Gene expression data will be log2 transformed.")
  
    ## "Detection p value is a concept (mainly) used by Illumina. There is a set of negative 
    ## controls on Illumina BeadArrays that are used to compute the background intensity. A 
    ## small detection p value indicates that the measured intensity is very likely to be a 
    ## true (significant) signal and not background noise. The intensity associated with a 
    ## certain detection p value corresponds to the respective upper quantile (percentile) of 
    ## the (probability) distribution of the intensities of the negative controls."
    ## https://www.researchgate.net/post/What_is_the_difference_between_detection_p_value_and_differential_p_value_in_microarray_data_analysis

    print("Setting gene expression p-value to",gxDetectionPvalue)

    ## testing
    #mtdf = loadMtData(methylFileName,sampleMethylAnnotationName)
    #sys.exit(33)

    ## cache the regionAggregateInfo
    #regionInfoH = readRegionAggregateInfoToDict(regionAggregateInfoName)

    ## load the gene expression data (<50Mb - smaller than Mt data)
    ## TODO: Purge any rows with NaNs in Gx data
    gxdf = loadGxData(geneExpFileName,gxDetectionPvalue)

    ## annotate the Gx probes for genes names
    gxdf = updateGeneNames(probeAnnotationName,symbolAnnotationName,gxdf)

    ## deal with multiple measurements (probes) for one annotated symbol?
    ## by doing nothing, only the last one in the list will be retained

    ## iterate on the methylation data (1.5Gb raw)
    ## only include data for samples in the Gx dataset (ignore unmatched)
    ## TODO: Purge any rows with NaNs in Mt data
    mtdf = loadMtData(methylFileName,sampleMethylAnnotationName)

    ## create a quick hash of sample IDs for Mt data
    mtPidH={}
    for pid in list(mtdf.columns.values):
       mtPidH[pid]=1

    ## TODO: reduce the size of the df?
    ## https://www.dataquest.io/blog/pandas-big-data/

    ## merge the data by patient ID
    mergeA=[]
    p=0
    print("Merging gene expression and methylation data for sample common to both.")
    for mpid in gxdf:
        pid=str(mpid)
        if pid == "Probe_ID":
            continue

        if DEBUG: print("Merging for sample ID:",pid)

        ## check to see patient ID is in both datasets before moving on
        if not pid in mtPidH:
            if DEBUG: print("Sample ID",pid,"missing from Mt data. Ignoring.")
            continue
        s=0

        mergeH={}

        ## PID
        mergeH['Label']=pid
        if DEBUG: print("Started for sample ID:",pid)

        ## Gx data 
        mergeH['Transcript']={}
        g=0
        for i, row in gxdf.iterrows():
            #print("row:",pid,row.dtype,row['Probe_ID'],row[pid])
            iv = row[pid]
            if log2TxfmGx:
                iv = math.log2(row[pid])
            mergeH['Transcript'][row['Probe_ID']]=iv
            #if g > 5:
            #    print("***\n***\n")
            #    print("***WARNING. Pre-maturely terminating transcript. Are you debugging?")
            #    print("***\n***\n")
            #    break
            g+=1

        ## Mt data
        mergeH['Alpha']={}
        m=0
        for i, row in mtdf.iterrows():
            #print("row:",pid,row.dtype,row['Probe_ID'],row[pid])
            mergeH['Alpha'][row['region']]=row[pid]
            #if m > 5: 
            #    print("***\n***\n")
            #    print("***WARNING. Pre-maturely terminating methylation. Are you debugging?")
            #    print("***\n***\n")
            #    break
            m+=1

        #if p > 5:
        #    print("***\n***\n")
        #    print("***WARNING. Pre-maturely terminating merge. Are you debugging?")
        #    print("***\n***\n")
        #    break
        if DEBUG: print("Merged data for sample",pid,mergeH)
        mergeA.append(mergeH)
        ky = list(mergeH.keys())
        print("Completed merge for sample",p, pid, ky[0], len(mergeH[ky[0]]), ky[1], len(mergeH[ky[1]]), ky[2], len(mergeH[ky[2]]))
        sys.stdout.flush()
        p+=1
       
    if DEBUG: print("Final merged array:", mergeA)

    print("Merged",len(mergeA),"samples.")
    ## save the object
    f = open(mergedFileName,'w')
    f.write(str(mergeA).replace("'",'"'))
    f.close()


    ## save the pandas dataframe as json
    ## SEE: https://datatofish.com/export-pandas-dataframe-json/
    
#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])
