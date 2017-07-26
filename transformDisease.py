from __future__ import division
import math as mt
import numpy as NP
import sys, re, os, random, pickle
import cwgutils

with open('./NVD.data','r+') as f:
  W = pickle.load(f)
with open('./COM.data','r+') as f:
  C = pickle.load(f)
with open('./1711.bdt','r+') as f:
  BDT = pickle.load(f)

def captureBiomarkerData(mxpfile):
  

def parseExpressionData(mxpfile):
  dataPts = cwgutils.countCSV() - 1
  expVars = cwgutils.readColumn(mxpfile, 1, ' ')
  for var in expVars:
      

if __name__=="__main__":
  script, mxpFile, indication = sys.argv
  
