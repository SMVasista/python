import re
import json
import pickle
import os
import csv
from pathlib import Path

# intialize a dict with array elements as keys and value as 
# either defaultValue or a function call to that element

# ** Sanjay - when you use with open(fileName, 'r') clause - python closes the file when done - you do not need to close it explicitly
# ** move all utility kind of functions, which are not central to ICOG functionality to utility library like this - it makes the main
#    code smaller and more readble

def arrayToDict(keylist, defaultValue):
    if type(defaultValue).__name__ == 'function':
        return dict((key, defaultValue(key)) for key in keylist)
    else:
        return dict((key, defaultValue) for key in keylist)

def readLines(fileName):
    with open(fileName, 'r') as f:
        return  [line.strip() for line in f.readlines() if len(line.strip()) > 0]

# read lines from a file, and split each line using sep - the tokens are stripped of leading/traing whitespaces

def readLinesAndSplit(fileName, sep=","):
    return [map(lambda token: token.strip(), line.split(sep)) for line in readLines(fileName)]

# used for reading a CSV like file - returns array of attribute values
def readColumns(fileName, columns):
    tokensList = readLinesAndSplit(fileName)
    columnData = []
    for column in columns:
        data = [tokens[column] if column < len(tokens) else None for tokens in tokensList]
        columnData.append(data)
    return columnData

def readColumn(fileName, column):
    tokensList = readLinesAndSplit(fileName)
    return [tokens[column] for tokens in tokensList]

# pickle convenience functions - uses fileName instead of file handle
def pickleDump(fileName, obj, fileMode='w'):
    with open(fileName, fileMode) as f:
        return pickle.dump(obj,f)

def pickleLoad(fileName):
    with open(fileName, 'r') as f:
        return pickle.load(f)

# val could be a value - in that case equality == is used for filtering
# if it is a function - function(elem) value is checked
# in case of dict collection - value is checked for (key, value) pairs

def filter(collection, val):
    if isinstance(collection, dict):
        if type(val).__name__ == 'function':
            return dict((k, v) for k,v in collection.iteritems() if val(v))
        else:
            return dict((k, v) for k,v in collection.iteritems() if v == val)
    elif isinstance(collection, list):
        if type(val).__name__ == 'function':
            return  [elem for elem in collection if val(elem)]
        else:
            return  [elem for elem in collection if val == elem]

# compact - will remove all the null elements of array or all the key,value pairs from dict where value is null (None)

def compact(collection):
    if isinstance(collection, list):
        return [ elem for elem in collection if elem != None ]
    elif isinstance(collection, dict):
        return dict((k, v) for k,v in collection.iteritems() if v != None)

# take only the key, value pairs from given dict where key is in the given list

def pick(obj, keylist):
    return dict((k, v) for k,v in obj.iteritems() if k in keylist)


# check if dir exists, if not create it - functionality akin to mkdir -p (available in python 3 but not in 2)

def mkpdir(pathName):
    if not Path(pathName).exists():
        os.makedirs(pathName)

# assume UTF is default encoding and indent value is always 2

def jsonLoad(fileName):
    with open(fileName, 'r') as fp:
        return json.loads(fp.read())

def jsonDump(fileName, obj):
    with open(fileName, 'w') as fp:
        fp.write(json.dumps(obj, indent=2))

# converts ab_cd_ef_gh or AB_CD_EF_GH to abCdEfGh
def snakeToCamel(s):
    ss = ''.join([x[0].upper() + x[1:].lower() for x in s.split('_')])
    return ss[0].lower() + ss[1:]

# replaces only null (None) items in dst list with values from src at the same index - effectively dst = dst || src
def replaceNullWithValues(dst, src):
    for idx, value in enumerate(dst):
        if value == None and src[idx] != None:
            dst[idx] = src[idx]


def countColumnsInCSV(filename):
	with open(filename) as f:
		reader = csv.reader(f, delimiter=',', skipinitialspace=True)
		first_row = next(reader)
		num_cols = len(first_row)
	return num_cols

#This function specifically arrays: converts the 'array' into 'value' of key which is the first element of the inner array
def reHeadArray(array):
	key = array.pop(0)
	return key, array

def mapToFloat(array):
	mod = []
	for value in array:
		value = float(value)
		mod.append(value)
	return mod

#Will take column 1 and 2 as a dict
def readFiletoDict(location):
	key = readColumn(location, 0)
	value = readColumn(location, 1)
	return dict(zip(key, value))


#EOF
