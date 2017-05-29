import sys, os
import json
import argparse
import copy
import time
import numpy, math
from pathlib import Path

import cwgutils
import cwgpatient

currWorkLoc = os.path.dirname(os.path.realpath(__file__))
workDir = os.getcwd()
pwd = workDir
holderPatients = {}

def updatePathwayScores(pathwayDir, initialPathwayFile):
    pathwayFile = os.path.join(pathwayDir, 'list')
    #initialPathwayFile = os.path.join(initialPathwayDir, 'list')

    print pathwayFile, initialPathwayFile

    if Path(pathwayFile).exists():
        print 'Pathway Relationship Score data found. Proceeding...'
        return
    else:
        print 'default pathway file', pathwayFile, ' does not exist - looking at initial pathway file:', initialPathwayFile

    if not Path(initialPathwayFile).exists():
        print 'No list found in', initialPathwayFile, 'to update pathways'
        return

    listFile = open(pathwayFile, 'w')

    pathways = cwgutils.readColumn(initialPathwayFile, 0)
    initalPathwayDir = os.path.join(currWorkLoc, 'Gene_Pathways_Scores_Backup')
    for pathway in pathways:
        genes, gofScores, lofScores = cwgutils.readColumns(os.path.join(initialPathwayFile, pathway), [0, 1, 2])
        # if lofScore is None - set it to gofScore
        gofScores = [float(score) for score in gofScores]
        lofScores = [float(score) for score in lofScores]
        cwgutils.pickleDump(os.path.join(pathwayDir, pathway + '-GOF'), dict(zip(genes, gofScores)))
        cwgutils.pickleDump(os.path.join(pathwayDir, pathway + '-LOF'), dict(zip(genes, lofScores)))
        listFile.write(pathway + '\n')

    listFile.close()

def categorizeDrugResponse(patients, scoresDir, trainingDir, patientMSource, patientResponse):
    scoresFile = os.path.join(scoresDir, 'list')
    print patientMSource

    pathways = cwgutils.readLines(scoresFile)
    drugSpecificPathways = cwgutils.arrayToDict(pathways, defaultValue=0)

    trngFile = open(os.path.join(trainingDir, '.kbplist'), 'w')

    noiseList = []
    for patientName in patients:
	print patientName
	print patientMSource[patientName]
        patientProfile = cwgpatient.loadPatientProfile(patientName, scoresDir+'/list', patientMSource[patientName])
        patientProfile["response"] = patientResponse[patientName]
        patientProfile = cwgpatient.genPathwaysSignature(patientProfile, scoresDir, trainingDir)
	#Sumanth : genPathwaysSignature takes 3 arguments 4 given; 'source' is not a list of pathways, it is location of source pathway scores

        # if response is None for patient in holderPatients  then the value can get modified
        holderPatients[patientName] = patientProfile

        cwgutils.pickleDump(os.path.join(trainingDir, patientName + '.clf'), patientProfile)
        trngFile.write(patientName + '\n')

	patientPathwayFile = os.path.join(trainingDir, patientName)
	print patientPathwayFile

        patientLoad = cwgutils.pickleLoad(os.path.join(trainingDir, patientName + '.clf'))
	patientPathways = patientLoad["pathways"] 
        noiseList.append(len(patientPathways))

# 0.2 is a random value
	#**Sumanth : Update drugSpecificPathways scores only if you find it in the patient profile not otherwise
        for pathway in patientPathways:
            if patientProfile["response"] == 'R':
		if patientProfile["pathways"][pathway] > 0:
                	drugSpecificPathways[pathway] += 0.2
		elif patientProfile["pathways"][pathway] < 0:
			drugSpecificPathways[pathway] -= 0.2
            elif patientProfile["response"] == 'N':
                if patientProfile["pathways"][pathway] > 0:
                	drugSpecificPathways[pathway] -= 0.2
		elif patientProfile["pathways"][pathway] < 0:
			drugSpecificPathways[pathway] += 0.2
    print drugSpecificPathways

    if len(noiseList) > 0:
        noise = math.sqrt(numpy.mean(noiseList))

    trngFile.close()
    print noise
    pathwaysForDrug = cwgutils.filter(drugSpecificPathways, lambda val: (val >= 0.08*noise or val <= -0.08*noise))
    cwgutils.pickleDump(os.path.join(trainingDir, 'DrugPathways'), pathwaysForDrug)
    print pathwaysForDrug
    return pathwaysForDrug

# returns score
def clusterPatient(profile):
   # first pull a subset out of allPatients
    prevSubset = copy.deepcopy(holderPatients.values())

    for pathway in profile["pathwaysmap"].keys():
        currSubset = [patient for patient in prevSubset if  pathway in patient["pathwaysmap"] and profile["pathwaysmap"][pathway] * patient["pathwaysmap"][pathway] > 0]
        if (len(currSubset) > 0):
            prevSubSet = currSubset
        else:
            break

    score = 0
    similarity = []
    for patient in prevSubSet:
        score += (1 if patient["response"] == 'R' else -1 if patient["response"] == 'N' else 0)
        if len(patient["pathwaysmap"]) > len(profile["pathwaysmap"]): 
            update_score = float(score*len(apatient.pathwaysmap)/(0.01+len(element.pathwaysmap)))
            similarity.append(update_score)
        
    return score


def calibrateResponseCoefficients(patients, scoresDir, trainingDir, pathwaysForDrug, patientResponse, patientMsource):
    # drug pathway should exist but not the calibration data
    calibDir = os.path.join(trainingDir, 'CalibrationData')
    scoresFile = os.path.join(scoresDir, 'list')
    pathways = cwgutils.readLines(scoresFile)
    cwgutils.mkpdir(calibDir)

    if pathwaysForDrug == None:
        pathwaysForDrug = cwgutils.pickleLoad(os.path.join(trainingDir, 'DrugPathways'))

    responseProbability = cwgutils.arrayToDict(patients, 0.0)
    causalDepth = cwgutils.arrayToDict(patients, 0)

    nDrug = len(pathwaysForDrug)
    profiles = {}

    # assume holdePatients is global variable, which was set in previous call
    for patientName in patients:
	print 'Processing', patientName
        profile = cwgpatient.loadPatientProfile(patientName, scoresFile, patientMsource[patientName])
        profile["response"] = patientResponse[patientName]
        profile = cwgpatient.genPathwaysSignature(profile, scoresDir, calibDir)
        profiles[patientName] = profile
        nPatient = len(pathways)

        for drugPath in pathwaysForDrug.keys():
            if drugPath in profiles[patientName]["pathwaysmap"]:
                responseProbability[patientName] += ((2.0 / nPatient) * (2.0 / nDrug) * pathwaysForDrug[drugPath] * profile["pathwaysmap"][drugPath])
	print responseProbability[patientName]
	print profile["response"]

#        if len(profile["pathwaysmap"]) > 2:
#            causalDepth[patientName] = clusterPatient(profile)

    predictedResponse = {}
    accuracy = {}

    PCRList = [float(responseProbability[patientName]) for patientName in patients if profiles[patientName]["response"] == 'R']
#        PARList = [float(causalDepth[patientName]) for patientName in patients if profiles[patientName]["response"] == 'R']
    NCRList = [float(responseProbability[patientName]) for patientName in patients if profiles[patientName]["response"] == 'N']
#        NARList = [float(causalDepth[patientName]) for patientName in patients if profiles[patientName]["response"] == 'N']

    pcr = sum(PCRList) / len(PCRList)
#        par = sum(PARList) / len(PARList)
    ncr = sum(NCRList) / len(NCRList)
#        nar = sum(NARList) / len(NARList)

    xThresh = 0.5 * (min(PCRList) + max(NCRList))
#        yThresh = 0.5 * (par + nar)
    yThresh = 0
    print pcr, ncr

    return [xThresh, yThresh]

# diff is the difference between currProbableResponse and threshold

def optimizePathwayWeightagePA(patientPathways, pathwaysForDrug, beta, skew, patientResponse):
    if beta == 1:
        return

    for pathway in patientPathways.keys():
        if pathway in pathwaysForDrug:
            if pathwaysForDrug[pathway] >= 0.1 and pathwaysForDrug[pathway] < 0:
                pathwaysForDrug[pathway] = 0.1
            elif pathwaysForDrug[pathway] <= 0.1 and pathwaysForDrug[pathway] > 0:
                pathwaysForDrug[pathway] = -0.1

            if skew < 0 and patientResponse == 'R':
                if patientPathways[pathway] > 0 and pathwaysForDrug[pathway] > 0:
                    pathwaysForDrug[pathway] = pathwaysForDrug[pathway] * (1.25* float(1 + abs(skew)))
                elif patientPathways[pathway] > 0 and pathwaysForDrug[pathway] < 0:
                    pathwaysForDrug[pathway] = pathwaysForDrug[pathway] / (1.1 * float(1 + abs(skew)))
                elif patientPathways[pathway] < 0 and pathwaysForDrug[pathway] > 0:
                    pathwaysForDrug[pathway] = pathwaysForDrug[pathway] / (1.1 * float(1 + abs(skew)))
                elif patientPathways[pathway] < 0 and pathwaysForDrug[pathway] < 0:
                    pathwaysForDrug[pathway] = pathwaysForDrug[pathway] * (1.25* float(1 + abs(skew)))
            elif skew > 0 and patientResponse == 'N':
                if patientPathways[pathway] > 0 and pathwaysForDrug[pathway] > 0:
                    pathwaysForDrug[pathway] = pathwaysForDrug[pathway] / (1.1 * float(1 + abs(skew)))
                elif patientPathways[pathway] > 0 and pathwaysForDrug[pathway] < 0:
                    pathwaysForDrug[pathway] = pathwaysForDrug[pathway] * (1.25 * float(1 + abs(skew)))
                elif patientPathways[pathway] < 0 and pathwaysForDrug[pathway] > 0:
                    pathwaysForDrug[pathway] = pathwaysForDrug[pathway] * (1.25 * float(1 + abs(skew)))
                elif patientPathways[pathway] < 0 and pathwaysForDrug[pathway] < 0:
                    pathwaysForDrug[pathway] = pathwaysForDrug[pathway] / (1.1 * float(1 + abs(skew)))                        

def testPatientResponse(patients, patientResponse, patientMsource, scoresDir, trainingDir, pathwaysForDrug, xThresh, yThresh, update_mode=0):
    predictionDir = os.path.join(pwd, 'PredictionData')
    scoresFile = os.path.join(scoresDir, 'list')
    pathways = cwgutils.readLines(scoresFile)
    cwgutils.mkpdir(predictionDir)
    predFile = open(os.path.join(predictionDir, 'cluster.dat'), 'w')


    # if len(holderPatients.keys()) == 0:
        # patientNames = cwg.readLines(os.path.join(trainingDir, '.kbplist'))
        # for name in patientNames:
            # holderPatients[name] = pickleLoad(trainingDir, name + '.clf')


    kbpFile = open(os.path.join(trainingDir, '.kbplist'), 'a')
    testPatientResponseProbability = cwgutils.arrayToDict(patients, 0.0)
    testPatientCausalDepth = cwgutils.arrayToDict(patients, 0)
    predictedResponse = {}

    nDrug = len(pathwaysForDrug)
    profiles = {}

    for patientName in patients:
        profile = cwgpatient.loadPatientProfile(patientName, os.path.join(scoresDir, 'list'), patientMsource[patientName])
        profile = cwgpatient.genPathwaysSignature(profile, scoresDir, trainingDir)
	print patientResponse[patientName]
        profile["response"] = patientResponse[patientName]

        nPatient = len(pathways)
        for drugPath in pathwaysForDrug.keys():
            if drugPath in profile["pathwaysmap"]:
                testPatientResponseProbability[patientName] += ((2.0 / nPatient) * (2.0 / nDrug) * pathwaysForDrug[drugPath] * profile["pathwaysmap"][drugPath])
	print testPatientResponseProbability[patientName]
	print profile["response"]

        if testPatientResponseProbability[patientName] >= xThresh:
            predictedResponse[patientName] = 'R'
        elif testPatientResponseProbability[patientName] < yThresh:
            predictedResponse[patientName] = 'N'

        if patientName not in holderPatients.keys():
            holderPatients[patientName] = profile

#        if len(profile["pathwaysmap"].keys()) > 2:
#            testPatientCausalDepth[patientName] = clusterPatient(profile)

        # suffix = (response == None) ? '' : (response == 'R') ? '1' : '0'
        suffix = '' if profile["response"] == None else '1' if profile["response"] == 'R' else '0'
#        ss = patientName + ' ' + ('10.2f' % (testPatientResponseProbability[patientName] - xThresh)) + ' ' + ('%02d' % testPatientCausalDepth[patientName])
#        if len(suffix):
#            ss = ss + ' ' + suffix
#        predFile.write(ss + '\n')

        diff = testPatientResponseProbability[patientName] - xThresh
        if update_mode == 0:
            continue
        elif update_mode == 1:
            if profile["response"] != None:
                if predictedResponse[patientName] == profile["response"]:
                    optimizePathwayWeightagePA(pathways, pathwaysForDrug, 1, diff, profile["response"])
                else:
                    optimizePathwayWeightagePA(pathways, pathwaysForDrug, -1, diff, profile["response"])
            else:
                pass
        elif update_mode == 2:
            if profile["response"] != None:
                if predictedResponse[patientName] == profile["response"]:
                    OptimizePathwayWeightageHyp(pathways, pathwaysForDrug, 1, diff, profile["response"])
                else:
                    OptimizePathwayWeightageHyp(pathways, pathwaysForDrug, -1, diff, profile["response"])
            else:
                pass

        if patientName not in patientNames:
            pickleDump(os.path.join(trainingDir, patientName + '.clf'))
            kbpFile.write(patientName + '\n')

    return predictedResponse

def genIcogScore(fileName, prediction, drug):
    with open(fileName, 'a') as f:
        for k,v in prediction.item():
            val = '5' if v == 'R' else '-5' if v == 'N' else '0' 
            if runningMode:
                f.write(drug + ',' + val + '\n')
            else:
                f.write(k + ',' + val + '\n')

def processDrug(drug, mode):
    if mode != 0:
        return

    print 'process drug:', drug
    drugDir = os.path.join(currWorkLoc, drug)

    # cwgutils.mkpdir(drugDir)
    pathwayScoresDir = os.path.join(drugDir, 'UpdatedPathwaysScores')
    pathwayTrainingDir = os.path.join(drugDir, 'UpdatedPathwaysTraining')

    print 'drug dir:', drugDir, 'scores dir:', pathwayScoresDir, 'training dir:', pathwayTrainingDir

    cwgutils.mkpdir(pathwayScoresDir)
    cwgutils.mkpdir(pathwayTrainingDir)

    #### - here supply the initial pathway dir rather than asking for it using raw input
    updatePathwayScores(pathwayScoresDir, './Gene_Pathway_Score_Backup/list')
    #Sumanth : I am not sure about the objective of pathwayTrainingDir, I am assuming this is where the initial score files for pathways and list is kept, and this is parsed into the updatedPatwayScores.
    #Have modified the input to this function accordingly

    testFile        = os.path.join(drugDir, drug + '_test.csv')
    calibrationFile = os.path.join(drugDir, drug + '_calib.csv')
    trainingFile    = os.path.join(drugDir, drug + '_training.csv')

    print 'testFile:', testFile, 'calibrationFile:', calibrationFile, 'trainingFile:', trainingFile
    # have to initialize these variables before, otherwise their scopes gets bound to if statement below
    patients        = []
    patientResponse = {}
    patientMSource  = {}
    #Sumanth: patients is a list, response and Msource are dicts; else throws error that list indicies cannot be anything other than int

    if Path(trainingFile).exists():
        patients, response, msource = cwgutils.readColumns(trainingFile, [0, 1, 2])

#        if msource[0] == None:
#            msource = copy.deepcopy(response)

        patientResponse = dict(zip(patients, response))
        patientMsource  = dict(zip(patients, msource))
	print patientResponse, patientMsource
    else:
        print 'training data file:', trainingFile, ' does not exist'

    if not Path(os.path.join(pathwayTrainingDir, 'DrugPathways')).exists():
        pathwaysForDrug = categorizeDrugResponse(patients, pathwayScoresDir, pathwayTrainingDir, patientMsource, patientResponse)
    else:
        print 'Drug Pathways for ', drug, ' already exists: ', os.path.join(pathwayTrainingDir, 'DrugPathways')
	pathwaysForDrug = cwgutils.pickleLoad(os.path.join(pathwayTrainingDir,'DrugPathways'))
	print 'Drug Pathways data loaded'
        
    if Path(calibrationFile).exists():
        calibPatients, calibResponse, calibMsource = cwgutils.readColumns(calibrationFile, [0, 1, 2])

        calibPatientResponse = dict(zip(calibPatients, calibResponse))
        calibPatientMsource  = dict(zip(calibPatients, calibMsource))
        thresholds = calibrateResponseCoefficients(calibPatients, pathwayScoresDir, pathwayTrainingDir, pathwaysForDrug, calibPatientResponse, calibPatientMsource)
    else:
        print 'Calibration data file:', calibrationFile, ' does not exist'

# else use drugName, patientName to find pres-stored thresholds - [thx, thy] value
# explore reading from database

    if Path(testFile).exists():
        testPatients, testResponse, testMsource = cwgutils.readColumns(testFile, [0, 1, 2])

        testPatientActualResponse = dict(zip(testPatients, testResponse))
        testPatientMsource  = dict(zip(testPatients, testMsource))

        # -- make sure  that thresholds exist

        # update modes

        prediction = testPatientResponse(testPatients, testPatientActualResponse, testPatientMsource, pathwayScoresDir, pathwayTrainingDir, pathwaysForDrug, thresholds[0], thresholds[1], 0)
	print prediction

#        genIcogScore(os.path.join(workDir, 'icog_result.csv'), prediction, drug)


def parseArgs():
    parser = argparse.ArgumentParser(description='iCOG Execution...')
    parser.add_argument('-testlist',    type=str, help='Path of test_data.csv optional for mode 0', default=None)
    parser.add_argument('-deplocation', type=str, help='Dependency location',     required=True)
    # parser.add_argument('-runmode',     type=int, help='Running mode 0-R&D,1-MS', required=True)
    parser.add_argument('-runmode',     type=int, help='Running mode 0-R&D,1-MS', default=0)
    parser.add_argument('-drugfile',    type=str, help='Path of drug file',       required=True)
    parser.add_argument('-upmode',      type=int, help='Update Mode 0/1 drefault 0', default=0)

    # print help if program is called without arguments
    if len(sys.argv) == 1:
        sys.argv.append('-h') 

    return parser.parse_args()

if __name__ == '__main__':
    args = parseArgs()
    print 'Running in ', args.runmode, ' mode...'

    if Path(args.deplocation).is_dir() == False:
        print 'Error: Unable to find the dependency location'
        sys.exit(-1)

    if args.runmode and not Path(args.testlist).exists():
        print 'Unable to find input file at ', testlist
        sys.exit(-1)

    if not Path(args.drugfile).exists():
        print 'Drug names file does not exist .. exitting'
        sys.exit(-1)

    if args.runmode:
        with open(os.path.join(workDir, 'icog_moonshot_test_data.csv'),'w') as fh:
            fh.write("MOONSHOT_PATIENT,,%s\n" %(testlist))
        args.testlist = os.path.join(workDir, 'icog_moonshot_test_data.csv')

    currWorkLoc = args.deplocation

    icogOutputDir = os.path.join(workDir, 'icog_' + time.strftime('%Y_%m_%d_%H_%M'))

    print workDir, currWorkLoc, icogOutputDir

    cwgutils.mkpdir(icogOutputDir)


    # read the historical data first, called holderPatients here, which is an array of patients
    # each patient is an object described and handled in cwgpatient.py
    # this list contains the patients for which prediction program has already been run

#    holderPatients = cwgutils.jsonLoad(os.path.join(workDir, 'holder_patients.json'))

#    print json.dumps(holderPatients[0]), len(holderPatients)
    #Sumanth : Script breaking because it is not able to find holder_patients in the primary run, check and proceed...
    # the read the list of input drugs
    drugList = cwgutils.readLines(args.drugfile)
    print drugList

    for drug in drugList:
        processDrug(drug, args.runmode)
