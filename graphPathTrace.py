from __future__ import division
from sys import argv
from pathlib import Path
import re, json, pickle, os, csv, math, numpy
import cwgutils

def orderNodes(nodes, connections):
	cnxrank = {}
	rank = {}
	for elem in nodes:
		cnxrank[elem] = [0, 0]	
		for line in connections:
			bitQuery = re.split(' -- ',line)[0]
			if elem in bitQuery:
				cnxrank[elem][0] += 1
		for line in connections:
			bitQuery = re.split(' -- ',line)[1]
			if elem in bitQuery:
				cnxrank[elem][1] += 1
		rank[elem] = int(cnxrank[elem][0] * cnxrank[elem][1])
	print rank

def updateCnxWeight(nodes, matrix):
	for i in range(len(nodes)):
		for j in range(len(nodes)):
			pins = []
			if matrix[nodes[i]][nodes[j]] > 0:
				for k in range(len(nodes)):
					if matrix[nodes[j]][nodes[k]] > 0:
						pins.append(nodes[k])
			print pins
			for elem in pins:
				if nodes[i] != elem:
					if matrix[nodes[i]][elem] != 1:
						matrix[nodes[i]][elem] = matrix[nodes[j]][elem] + 1
#	print matrix
				

def writeDictList(data):
	lines = cwgutils.readLines(data)
	cnxlist = list(set(cwgutils.readColumn(data, 0)))
	nodelist = []
	for elem in cnxlist:
		nodelist.append(re.split(' -- ', elem)[0])
	nodelist = list(set(nodelist))
	cnxarray = {}
	for elem in nodelist:
		cnxarray[elem] = {}
	for line in lines:
		element1 = re.split(' -- ', line)[0]
		element2 = re.split(' -- ', line)[1]
		cnxarray[element1][element2] = 1
	for elem1 in nodelist:
		for elem2 in cnxarray:
			try:
				if cnxarray[elem1][elem2] == 1:
					continue
			except:
				if elem1 == elem2:
					 cnxarray[elem1][elem2] = 1
				else:
					cnxarray[elem1][elem2] = 0
	return nodelist, cnxlist, cnxarray

if __name__=="__main__":
	script, datafile = argv
	Nodes, Connections, Matrix = writeDictList(datafile)
#	print Matrix
	orderedNodes = orderNodes(Nodes, Connections)
#	updateCnxWeight(Nodes, Matrix)
