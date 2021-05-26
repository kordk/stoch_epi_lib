#####
#####   Script to produce RMSE of linear models using scikit learn, assuming 100 data splits located in
#####       the directory dataFolder (given as required command line argument)
#####       Each split should be located in a numbered subdirectory and include a train.json and test.json file.



import numpy as np
from sklearn import linear_model as lm
import json
import pandas as pd
import sys


existingTable = pd.read_csv("ModelRMSE.csv",index_col = 0)

NumFolds = 100
parameters = {}


dataFolder = sys.argv[1]
network_bindingSite_file = sys.argv[2] #should be a TSV with column for gene and binding region (as named in data file)

with open(dataFolder + "/0/train.json", "r") as file:
    TrainData = json.load(file)

with open(dataFolder + "/0/test.json", "r") as file:
    TestData = json.load(file)


FullData = TrainData + TestData
print("Total samples (training + testing)", len(FullData))

SiteOrder = list(FullData[0]["Alpha"].keys())#np.unique([ky for TD in FullData for ky in list(TD["Alpha"].keys())])
TranscOrder = list(existingTable.index)

trainDataIn = np.array([[dat["Alpha"][ky] if ky in dat["Alpha"].keys() else 0 for ky in SiteOrder] for dat in FullData])
trainDataOut = np.array([[dat["Transcript"][ky] if ky in dat["Transcript"].keys() else 0 for ky in TranscOrder] for dat in FullData])


Lasso = lm.MultiTaskLassoCV(max_iter = 3000)
Lasso.fit(trainDataIn,trainDataOut)
parameters["MTLASSOalpha"] = Lasso.alpha_

print("Fit MT LASSO alpha")

trainDataOut_byGene = trainDataOut.T

regionInfo = pd.read_csv(network_bindingSite_file,sep = '\t')


network_bygene = {}
for gn in TranscOrder:
    ### get the info for that gene
    sites = regionInfo[regionInfo.gene == gn]['region'].values
    # print([si in SiteOrder for si in sites])
    network_bygene[gn] = [SiteOrder.index(si) for si in sites if si in SiteOrder]



LassoSepParams = {}
NetLassoSepParams = {}
for i in range(len(trainDataOut_byGene)):
    #using all sites
    LassoSep = lm.LassoCV(max_iter = 3000)
    LassoSep.fit(trainDataIn,trainDataOut_byGene[i])
    LassoSepParams[i] = LassoSep.alpha_
    #using only connected sites
    if len(network_bygene[TranscOrder[i]]):
        trainIntmp = trainDataIn[:,network_bygene[TranscOrder[i]]]
        LassoNetSep = lm.LassoCV(max_iter = 3000)
        LassoNetSep.fit(trainIntmp,trainDataOut_byGene[i])
        NetLassoSepParams[i] = LassoNetSep.alpha_


    print("Fit LASSO alpha, gene" + str(i))


parameters["LASSOalpha"] = LassoSepParams
parameters["NetLASSOalpha"] = NetLassoSepParams


ElasticNet = lm.MultiTaskElasticNetCV(max_iter = 3000)
ElasticNet.fit(trainDataIn,trainDataOut)
parameters["MTENalpha"] = ElasticNet.alpha_
parameters["MTENl1R"] = ElasticNet.l1_ratio_

print("Fit MT EN alpha & l1 Ratio")

ENSep_Alaphas = {}
ENSep_l1Rs = {}
NetENSep_Alaphas = {}
NetENSep_l1Rs = {}

for i in range(len(trainDataOut_byGene)):
    ENSep = lm.ElasticNetCV(max_iter = 3000)
    ENSep.fit(trainDataIn,trainDataOut_byGene[i])
    ENSep_Alaphas[i] = ENSep.alpha_
    ENSep_l1Rs[i] = ENSep.l1_ratio_

    if len(network_bygene[TranscOrder[i]]):
        trainIntmp = trainDataIn[:,network_bygene[TranscOrder[i]]]
        NetENSep = lm.ElasticNetCV(max_iter = 3000)
        NetENSep.fit(trainIntmp,trainDataOut_byGene[i])
        NetENSep_Alaphas[i] = NetENSep.alpha_
        NetENSep_l1Rs[i] = NetENSep.l1_ratio_


    print("Fit EN alpha & l1 Ratio, gene" + str(i))


parameters["ENalpha"] = ENSep_Alaphas
parameters["ENl1R"] = ENSep_l1Rs

parameters["NetENalpha"] = NetENSep_Alaphas
parameters["NetENl1R"] = NetENSep_l1Rs

# print(parameters)
#######################################

LassoRMSE = pd.DataFrame(index = existingTable.index, columns = range(NumFolds))
LassoRMSESep = pd.DataFrame(index = existingTable.index, columns = range(NumFolds))
NetLassoRMSESep = pd.DataFrame(index = existingTable.index, columns = range(NumFolds))
ElasticNetRMSE = pd.DataFrame(index = existingTable.index, columns = range(NumFolds))
ElasticNetRMSESep = pd.DataFrame(index = existingTable.index, columns = range(NumFolds))
NetElasticNetRMSESep = pd.DataFrame(index = existingTable.index, columns = range(NumFolds))


sys.stdout.write("[%s]" % (" " * NumFolds))
sys.stdout.flush()
sys.stdout.write("\b" * (NumFolds+1)) # return to start of line, after '['

failedon = []

for j in range(NumFolds):

    try:

        with open("../eval/2021-05-13/28925389-PTSD+Ca_v_cntl/"+ str(j) + "/train.json", "r") as file:
            TrainData = json.load(file)

        with open("../eval/2021-05-13/28925389-PTSD+Ca_v_cntl/"+ str(j) + "/test.json", "r") as file:
            TestData = json.load(file)

        trainDataIn = np.array([[dat["Alpha"][ky] if ky in dat["Alpha"].keys() else 0 for ky in SiteOrder] for dat in TrainData])
        trainDataOut = np.array([[dat["Transcript"][ky] if ky in dat["Transcript"].keys() else 0 for ky in TranscOrder] for dat in TrainData])
        trainDataOut_byGene = trainDataOut.T


        testDataIn = np.array([[dat["Alpha"][ky] if ky in dat["Alpha"].keys() else 0 for ky in SiteOrder] for dat in TestData])
        testDataOut = np.array([[dat["Transcript"][ky] if ky in dat["Transcript"].keys() else 0 for ky in TranscOrder] for dat in TestData])

        Lasso = lm.MultiTaskLasso(alpha = parameters["MTLASSOalpha"])
        Lasso.fit(trainDataIn,trainDataOut)
        predictedTranscriptLasso = np.array(Lasso.predict(testDataIn))

        SqErrorsLasso =(testDataOut - predictedTranscriptLasso)**2
        rtmnSqErrorByGeneLasso = np.sum(SqErrorsLasso, axis = 0)**(0.5)
        rtmnSqDictLasso = dict([(TranscOrder[i],rtmnSqErrorByGeneLasso[i]) for i in range(len(rtmnSqErrorByGeneLasso))])

        for ky in LassoRMSE.index:
            LassoRMSE.loc[ky,j] = rtmnSqDictLasso[ky]

        predictedTranscriptLassoSep = np.empty_like(predictedTranscriptLasso)
        predictedTranscriptNetLassoSep = np.zeros_like(predictedTranscriptLasso)


        for i in range(len(trainDataOut_byGene)):
            LassoSep = lm.Lasso(alpha = parameters["LASSOalpha"][i])
            LassoSep.fit(trainDataIn,trainDataOut_byGene[i])
            predictedTranscriptLassoSep[:,i] = np.array(LassoSep.predict(testDataIn))
            if len(network_bygene[TranscOrder[i]]):
                trainIntmp = trainDataIn[:,network_bygene[TranscOrder[i]]]
                testTmp = testDataIn[:,network_bygene[TranscOrder[i]]]
                LassoNetSep = lm.Lasso(alpha = parameters["NetLASSOalpha"][i])
                LassoNetSep.fit(trainIntmp,trainDataOut_byGene[i])
                predictedTranscriptNetLassoSep[:,i] = np.array(LassoNetSep.predict(testTmp))

        SqErrorsLassoSep = (testDataOut - predictedTranscriptLassoSep)**2
        rtmnSqErrorByGeneLassoSep = np.sum(SqErrorsLassoSep, axis = 0)**(0.5)
        rtmnSqDictLassoSep = dict([(TranscOrder[i],rtmnSqErrorByGeneLassoSep[i]) for i in range(len(rtmnSqErrorByGeneLassoSep))])

        SqErrorsNetLassoSep = (testDataOut - predictedTranscriptNetLassoSep)**2
        rtmnSqErrorByGeneNetLassoSep = np.sum(SqErrorsNetLassoSep, axis = 0)**(0.5)
        rtmnSqDictNetLassoSep = dict([(TranscOrder[i],rtmnSqErrorByGeneNetLassoSep[i]) for i in range(len(rtmnSqErrorByGeneNetLassoSep))])

        for ky in LassoRMSESep.index:
            LassoRMSESep.loc[ky,j] = rtmnSqDictLassoSep[ky]

        for ky in NetLassoRMSESep.index:
            if len(network_bygene[ky]):
                NetLassoRMSESep.loc[ky,j] = rtmnSqDictNetLassoSep[ky]
            else:
                NetLassoRMSESep.loc[ky,j] = 'NA'

        ElasticNet = lm.MultiTaskElasticNet(alpha = parameters["MTENalpha"], l1_ratio = parameters["MTENl1R"])
        ElasticNet.fit(trainDataIn,trainDataOut)
        predictedTranscript_ElasticNet = np.array(ElasticNet.predict(testDataIn))

        SqErrors_EN = (testDataOut - predictedTranscript_ElasticNet)**2
        rtmnSqErrorByGene_EN = np.sum(SqErrors_EN, axis = 0)**(0.5)
        rtmnSqDict_EN = dict([(TranscOrder[i],rtmnSqErrorByGene_EN[i]) for i in range(len(rtmnSqErrorByGene_EN))])

        for ky in ElasticNetRMSE.index:
            ElasticNetRMSE.loc[ky,j] = rtmnSqDict_EN[ky]


        predictedTranscript_ENSep = np.empty_like(predictedTranscriptLasso)
        predictedTranscript_NetENSep = np.zeros_like(predictedTranscriptLasso)


        for i in range(len(trainDataOut_byGene)):
            ElasticNetSep = lm.ElasticNet(alpha = parameters["ENalpha"][i], l1_ratio = parameters["ENl1R"][i])
            ElasticNetSep.fit(trainDataIn,trainDataOut_byGene[i])
            predictedTranscript_ENSep[:,i] = np.array(ElasticNetSep.predict(testDataIn))

            if len(network_bygene[TranscOrder[i]]):
                trainIntmp = trainDataIn[:,network_bygene[TranscOrder[i]]]
                testTmp = testDataIn[:,network_bygene[TranscOrder[i]]]
                NetElasticNetSep = lm.ElasticNet(alpha = parameters["NetENalpha"][i], l1_ratio = parameters["NetENl1R"][i])
                NetElasticNetSep.fit(trainIntmp,trainDataOut_byGene[i])
                predictedTranscript_NetENSep[:,i] = np.array(NetElasticNetSep.predict(testTmp))

        SqErrors_ENSep = (testDataOut - predictedTranscript_ENSep)**2
        rtmnSqErrorByGene_ENSep = np.sum(SqErrors_ENSep, axis = 0)**(0.5)
        rtmnSqDict_ENSep = dict([(TranscOrder[i],rtmnSqErrorByGene_ENSep[i]) for i in range(len(rtmnSqErrorByGene_ENSep))])

        SqErrors_NetENSep = (testDataOut - predictedTranscript_NetENSep)**2
        rtmnSqErrorByGene_NetENSep = np.sum(SqErrors_NetENSep, axis = 0)**(0.5)
        rtmnSqDict_NetENSep = dict([(TranscOrder[i],rtmnSqErrorByGene_NetENSep[i]) for i in range(len(rtmnSqErrorByGene_NetENSep))])


        for ky in ElasticNetRMSESep.index:
            ElasticNetRMSESep.loc[ky,j] = rtmnSqDict_ENSep[ky]

        for ky in NetElasticNetRMSESep.index:
            if len(network_bygene[ky]):
                NetElasticNetRMSESep.loc[ky,j] = rtmnSqDict_NetENSep[ky]
            else:
                NetElasticNetRMSESep.loc[ky,j] = 'NA'

    except:
        failedon += [j]
        LassoRMSE.drop(j,axis = 1, inplace = True)
        LassoRMSESep.drop(j,axis = 1, inplace = True)
        NetLassoRMSESep.drop(j,axis = 1, inplace = True)
        ElasticNetRMSE.drop(j,axis = 1, inplace = True)
        ElasticNetRMSESep.drop(j,axis = 1, inplace = True)
        NetElasticNetRMSESep.drop(j,axis = 1, inplace = True)

    sys.stdout.write("-")
    sys.stdout.flush()

sys.stdout.write("]\n") # this ends the progress bar

if len(failedon):
    print("Failed on folds (possibly missing data): ",failedon)

####
for ky in NetLassoRMSESep.index:
    if not len(network_bygene[ky]):
        NetLassoRMSESep.drop(ky,axis = 0,inplace = True)
for ky in NetElasticNetRMSESep.index:
    if not len(network_bygene[ky]):
        NetElasticNetRMSESep.drop(ky,axis = 0,inplace = True)

LassoRMSE.to_csv("MTLassoRMSE.csv")
LassoRMSESep.to_csv("LassoRMSE.csv")
NetLassoRMSESep.to_csv("NetworkLassoRMSE.csv")
ElasticNetRMSE.to_csv("MTElasticNetRMSE.csv")
ElasticNetRMSESep.to_csv("ElasticNet.csv")
NetElasticNetRMSESep.to_csv("NetworkElasticNet.csv")



print("MT LASSO\n",LassoRMSE.describe())
print("LASSO\n",LassoRMSESep.describe())
print("Network LASSO\n",NetLassoRMSESep.describe())
print("MT Elastic Net\n",ElasticNetRMSE.describe())
print("Elastic Net\n",ElasticNetRMSESep.describe())
print("Network Elastic Net\n",NetElasticNetRMSESep.describe())
