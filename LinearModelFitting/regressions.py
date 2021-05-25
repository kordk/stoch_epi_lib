#####
#####   Script to produce RMSE of linear models using scikit learn, assuming 100 data splits located in
#####       the directory dataFolder (given as required command line argument)
#####       Each split should be located in a numbered subdirectory and include a train.json and test.json file.


import numpy as np
from sklearn import linear_model as lm
import json
import pandas as pd
import sys


existingTable = pd.read_csv("table.csv")

NumFolds = 100
parameters = {}


dataFolder = sys.argv[1]

with open(dataFolder + "/0/train.json", "r") as file:
    TrainData = json.load(file)

with open(dataFolder + "/0/test.json", "r") as file:
    TestData = json.load(file)

FullData = TrainData + TestData

SiteOrder = np.unique([ky for TD in FullData for ky in list(TD["Alpha"].keys())])
TranscOrder = list(existingTable.GeneSymbol)

trainDataIn = np.array([[dat["Alpha"][ky] if ky in dat["Alpha"].keys() else 0 for ky in SiteOrder] for dat in FullData])
trainDataOut = np.array([[dat["Transcript"][ky] if ky in dat["Transcript"].keys() else 0 for ky in TranscOrder] for dat in FullData])

Lasso = lm.MultiTaskLassoCV(max_iter = 3000)
Lasso.fit(trainDataIn,trainDataOut)
parameters["MTLASSOalpha"] = Lasso.alpha_

print("Fit MT LASSO alpha")

trainDataOut_byGene = trainDataOut.T
LassoSepParams = {}
for i in range(len(trainDataOut_byGene)):
    LassoSep = lm.LassoCV(max_iter = 3000)
    LassoSep.fit(trainDataIn,trainDataOut_byGene[i])
    LassoSepParams[i] = LassoSep.alpha_
    print("Fit LASSO alpha, gene" + str(i))


parameters["LASSOalpha"] = LassoSepParams


ElasticNet = lm.MultiTaskElasticNetCV(max_iter = 3000)
ElasticNet.fit(trainDataIn,trainDataOut)
parameters["MTENalpha"] = ElasticNet.alpha_
parameters["MTENl1R"] = ElasticNet.l1_ratio_

print("Fit MT EN alpha & l1 Ratio")

ENSep_Alaphas = {}
ENSep_l1Rs = {}

for i in range(len(trainDataOut_byGene)):
    ENSep = lm.ElasticNetCV(max_iter = 3000)
    ENSep.fit(trainDataIn,trainDataOut_byGene[i])
    ENSep_Alaphas[i] = ENSep.alpha_
    ENSep_l1Rs[i] = ENSep.l1_ratio_
    print("Fit EN alpha & l1 Ratio, gene" + str(i))


parameters["ENalpha"] = ENSep_Alaphas
parameters["ENl1R"] = ENSep_l1Rs

print(parameters)
#######################################

LassoRMSE = pd.DataFrame(index = existingTable.GeneSymbol, columns = range(NumFolds))
LassoRMSESep = pd.DataFrame(index = existingTable.GeneSymbol, columns = range(NumFolds))
ElasticNetRMSE = pd.DataFrame(index = existingTable.GeneSymbol, columns = range(NumFolds))
ElasticNetRMSESep = pd.DataFrame(index = existingTable.GeneSymbol, columns = range(NumFolds))

sys.stdout.write("[%s]" % (" " * NumFolds))
sys.stdout.flush()
sys.stdout.write("\b" * (NumFolds+1)) # return to start of line, after '['

for j in range(NumFolds):

    with open(dataFolder + "/"+ str(j) + "/train.json", "r") as file:
        TrainData = json.load(file)

    with open(dataFolder + "/"+ str(j) + "/test.json", "r") as file:
        TestData = json.load(file)

    SiteOrder = np.unique([ky for TD in TrainData for ky in list(TD["Alpha"].keys())])

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

    for i in range(len(trainDataOut_byGene)):
        LassoSep = lm.Lasso(alpha = parameters["LASSOalpha"][i])
        LassoSep.fit(trainDataIn,trainDataOut_byGene[i])
        predictedTranscriptLassoSep[:,i] = np.array(LassoSep.predict(testDataIn))

    SqErrorsLassoSep = (testDataOut - predictedTranscriptLassoSep)**2
    rtmnSqErrorByGeneLassoSep = np.sum(SqErrorsLassoSep, axis = 0)**(0.5)
    rtmnSqDictLassoSep = dict([(TranscOrder[i],rtmnSqErrorByGeneLassoSep[i]) for i in range(len(rtmnSqErrorByGeneLassoSep))])

    for ky in LassoRMSESep.index:
        LassoRMSESep.loc[ky,j] = rtmnSqDictLassoSep[ky]

    ElasticNet = lm.MultiTaskElasticNet(alpha = parameters["MTENalpha"], l1_ratio = parameters["MTENl1R"])
    ElasticNet.fit(trainDataIn,trainDataOut)
    predictedTranscript_ElasticNet = np.array(ElasticNet.predict(testDataIn))

    SqErrors_EN = (testDataOut - predictedTranscript_ElasticNet)**2
    rtmnSqErrorByGene_EN = np.sum(SqErrors_EN, axis = 0)**(0.5)
    rtmnSqDict_EN = dict([(TranscOrder[i],rtmnSqErrorByGene_EN[i]) for i in range(len(rtmnSqErrorByGene_EN))])

    for ky in ElasticNetRMSE.index:
        ElasticNetRMSE.loc[ky,j] = rtmnSqDict_EN[ky]


    predictedTranscript_ENSep = np.empty_like(predictedTranscriptLasso)

    for i in range(len(trainDataOut_byGene)):
        ElasticNetSep = lm.ElasticNet(alpha = parameters["ENalpha"][i], l1_ratio = parameters["ENl1R"][i])
        ElasticNetSep.fit(trainDataIn,trainDataOut_byGene[i])
        predictedTranscript_ENSep[:,i] = np.array(ElasticNetSep.predict(testDataIn))

    SqErrors_ENSep = (testDataOut - predictedTranscript_ENSep)**2
    rtmnSqErrorByGene_ENSep = np.sum(SqErrors_ENSep, axis = 0)**(0.5)
    rtmnSqDict_ENSep = dict([(TranscOrder[i],rtmnSqErrorByGene_ENSep[i]) for i in range(len(rtmnSqErrorByGene_ENSep))])

    for ky in ElasticNetRMSESep.index:
        ElasticNetRMSESep.loc[ky,j] = rtmnSqDict_ENSep[ky]


    sys.stdout.write("-")
    sys.stdout.flush()

sys.stdout.write("]\n") # this ends the progress bar

LassoRMSE.to_csv("MTLassoRMSE.csv")
LassoRMSESep.to_csv("LassoRMSE.csv")
ElasticNetRMSE.to_csv("MTElasticNetRMSE.csv")
ElasticNetRMSESep.to_csv("ElasticNet.csv")

print("MT LASSO\n",LassoRMSE.describe())
print("LASSO\n",LassoRMSESep.describe())
print("MT Elastic Net\n",ElasticNetRMSE.describe())
print("Elastic Net\n",ElasticNetRMSESep.describe())
