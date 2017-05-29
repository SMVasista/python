import cwgutils
import os
from pathlib import Path

# ** Sanjay - the expression >>for x in dict<< is same as >>for x in dict.keys()<< - but second way is preferable as it clarifies that
#             it is a dict not an array

targetLOFGOF = {}

def loadPatientProfile(name, pathwaysFile, geneFile):
    # initialize a dict

    profile = { "type": None, "id": name, "gender": None, "response": None, "disease": None}
    #if pathways == None:
    profile["pathways"] = cwgutils.arrayToDict(cwgutils.readLines(pathwaysFile), defaultValue=0)
    #else:
     #   profile["pathways"] = pathways

    # load mutaions in LOF/GOF format and construct disease object
    geneList, mutationList = cwgutils.readColumns(geneFile, [0, 1])
    transformedMutations = map(lambda mut: 'GOF' if mut == 'OE' else 'LOF' if mut == 'KD' else None, mutationList)

    profile["disease"]   = dict(zip(geneList, mutationList))
    profile["mutations"] = dict(zip(geneList, transformedMutations))
    return profile


def loadLOFGOFData(profile, source):
    # load it only once
    if len(targetLOFGOF.keys()) > 0:
        return targetLOFGOF

    print source
    for pathway in profile["pathways"].keys():
        targetLOFGOF[pathway] = {
                "LOF": cwgutils.pickleLoad(os.path.join(source, pathway+'-LOF')),
                "GOF": cwgutils.pickleLoad(os.path.join(source, pathway+'-GOF'))
                }
    return targetLOFGOF

def dumpProfile(profile, dest, mode):
    if mode == 0 or re.match(r'PredictionData', dest):
        if len(profile["pathwaysmap"]) == 0:
            print "Patient does not express any strong pathways signatures... Proceeding with weak signal..."
            cwgutils.pickleDump(profile["pathways"], os.path.join(dest, profile["id"]))
        else:
            cwgutils.pickleDump(profile["pathwaysmap"], os.path.join(dest, profile["id"]))

    if mode == 0:   # not Moonshot
        with open(os.path.join(dest, 'list'), 'a') as f:
            f.write(profile["id"] + '\n')

# generate pathways signature
# first load all the target LOF and GOF from the files

def genPathwaysSignature(profile, source, dest):
    target = loadLOFGOFData(profile, source)
#** Sumanth: Logic was wrong here, corrected
    for gene in profile["mutations"].keys():
        for pathway in profile["pathways"].keys():
            if gene in target[pathway]["GOF"].keys() and profile["mutations"][gene] == 'GOF':
                profile["pathways"][pathway] += target[pathway]["GOF"][gene]
            elif gene in target[pathway]["LOF"].keys() and profile["mutations"][gene] == 'LOF':
                profile["pathways"][pathway] -= target[pathway]["LOF"][gene]

    profile["pathwaysmap"] = cwgutils.filter(profile["pathways"], lambda x: x >= 0.75 or x <= -0.5)
    print profile["pathwaysmap"]
    return profile
