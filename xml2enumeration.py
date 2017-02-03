#!/usr/bin/python3
""" Creates gene feature enumeration file for any given gene from a .xml file of IMGT.
Create a folder and inside of it add this script + the .xml file from IMGT.
On the command line change to the directory with the files and type:
python3 scriptName xmlFile genename(for example DPB1)"""

import xml.etree.ElementTree as ET
import re
import os
import sys

getdirectory = os.getcwd()
os.chdir(getdirectory)

if len(sys.argv) != 3:
        print('usage:python3 scriptName xmlFile genename(for example DPB1)')
        sys.exit()

xmlFile = sys.argv[1]
geneName = sys.argv[2]

def initialize_dict(header_line, separator="\t"):
    header_names = header_line.strip().split(separator)
    data_dict = {}
    for column in header_names:
        data_dict[column] = {}
    return data_dict

def parse_header(header_line, separator="\t"):
    return header_line.strip().split(separator)

name = 'name="HLA-'+ geneName
filein = open(xmlFile, 'r')
gene = re.compile(r'\b'+name+'\S+')#finds all alleles in the xmlFile for the gene you examine
alleles = []
for line in filein:
    line = line.strip().split('\t')
    mo = gene.search(str(line))
    if mo:
        so = mo.group()
        alleles.append(so.split('"')[1])
filein.close()
alSeqName = geneName+'_SeqPerFeature.txt'#creates _SeqPerFeature.txt file with the alleles and the sequence of their gene features. 0 is used if a feature is not sequenced.
alleleSeq = open(alSeqName,'w')
tree = ET.parse(xmlFile)
allele = tree.findall('{http://hla.alleles.org/xml}allele')
genefeatures =[]
for x in allele:
    if x.attrib['name'] == alleles[0]:
        feature = x.findall('{http://hla.alleles.org/xml}sequence/{http://hla.alleles.org/xml}feature')
        for f in feature:
            coord= f.find('{http://hla.alleles.org/xml}SequenceCoordinates')
            genefeatures.append(f.attrib['name'])
for z in alleles:
    for x in allele:
        if x.attrib['name'] == z:
            feature = x.findall('{http://hla.alleles.org/xml}sequence/{http://hla.alleles.org/xml}feature')
            testfeatures = []
            for f in feature:
                testfeatures.append(f.attrib['name'])
                if len(testfeatures) > len(genefeatures):
                    genefeatures = testfeatures
if genefeatures[-1:] == ['Translation']:
      del genefeatures[-1:]

alleleSeq.write("alleleName"+'\t')
for item in genefeatures:
    alleleSeq.write(item+'\t')
alleleSeq.write('\n')
print('alleles removed from earlier versions'+'\n'+str(0))
for z in alleles:
    for x in allele:
        dictionary = {}
        if x.attrib['name']==z:
            try:
                seq = x.find('{http://hla.alleles.org/xml}sequence/{http://hla.alleles.org/xml}nucsequence').text
                alleleSeq.write(z+'\t')
            except AttributeError:
                print(z)
                continue
            feature = x.findall('{http://hla.alleles.org/xml}sequence/{http://hla.alleles.org/xml}feature')
            for f in feature:
                coord = f.find('{http://hla.alleles.org/xml}SequenceCoordinates')
                try:
                    dictionary[f.attrib['name']] = str(seq[(int(coord.attrib['start'])-1):int(coord.attrib['end'])])
                except AttributeError:
                    continue
            for y in genefeatures:
                if y in dictionary:
                    alleleSeq.write(dictionary[y]+'\t')
                else:
                    alleleSeq.write(str(0)+'\t')
            alleleSeq.write('\n')
alleleSeq.close()
enumerName = geneName +'enumeration.txt'
input_file = open(alSeqName, 'r')
output_file = open(enumerName, 'w')
header_line = input_file.readline()
header_reference = parse_header(header_line)
data_dict = initialize_dict(header_line)
output_file.write(header_line)
for line in input_file:#creates the enumeration.txt file with all the alleles and their gene features enumerated. _SeqPerFeature file is used.
        line_data = line.strip().split("\t")
        output_file.write(line_data[0]+'\t')
        for index, element in enumerate(line_data[1:]):
            try:
                if int(element) == 0:
                    output_file.write(str(0)+'\t')
            except ValueError:
                column_name = header_reference[index]
                if element in data_dict[column_name]:
                    output_file.write(str(data_dict[column_name][element])+'\t')
                else:
                    value = len(data_dict[header_reference[index]]) + 1
                    output_file.write(str(value)+'\t')
                    data_dict[column_name][element] = value
        output_file.write('\n')
output_file.close()
