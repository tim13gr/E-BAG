#!/usr/bin/python3
"""This script calculates the possible AmbiguousGenotypes for a given pair of alleles and a set of given features.
It uses the enumerationfile created from xml2enumeration.py.
On the command line change to the directory with the files, then type:
python3 scriptName alleleName1 alleleName2 features(example:5UTR,E1,I1,E2,3UTR) inputfile outputfile
type allele names without the HLA-. example DPB1*13:01:01')"""
import itertools
from itertools import combinations
import sys
import os

getdirectory = os.getcwd()
os.chdir(getdirectory)

if len(sys.argv) != 6:
        print('usage: python3 scriptName alleleName1 alleleName2 features(example:5UTR,E1,I1,E2,3UTR) inputfile outputfile' +'\n'+'type allele names without the HLA-. for example DPB1*13:01:01')
        sys.exit()
allname1 = sys.argv[1]
allname2 = sys.argv[2]
inputfeat = sys.argv[3].split(',')
infile = sys.argv[4]
outname = sys.argv[5]

enumerationfile = open(infile, 'r')
outfile = open(outname,'w')
al1 = 'HLA-'+ allname1
al2 = 'HLA-'+ allname2
outfile.write('InputGenotype'+'\t'+'AmbiguousGenotypes'+'\t'+'GeneFeaturesExamined'+'\t'+'InputEnumerationUsed'+'\t'+'AmbiguousGenotypesEnumeration'+'\n')
featDic = {"5' UTR":'5UTR',"3' UTR":'3UTR', 'Exon 1':'E1', 'Intron 1':'I1','Exon 2':'E2', 'Intron 2':'I2','Exon 3':'E3', 'Intron 3':'I3','Exon 4':'E4', 'Intron 4':'I4','Exon 5':'E5', 'Intron 5':'I5','Exon 6':'E6', 'Intron 6':'I6','Exon 7':'E7', 'Intron 7':'I7','Exon 8':'E8', 'Intron 8':'I8','Exon 9':'E9', 'Intron 9':'I9','Exon 10':'E10', 'Intron 10':'I10','Exon 11':'E11', 'Intron 11':'I11'}

def allele_list(elements): # calculates the possible alleles that match the input alleles gfe for the features provided
        for element in elements:
                a = int(line[(int(element)+1)])
                b = int(allel1[(int(element)+1)])
                c = int(allel2[(int(element)+1)])
                if (
                        (a == b) or
                        (a == c) or
                        (a == 0)
                        ) is False:
                        return False
                return True

def possible_pairs(elements):# applies specific criteria for all the possible combinations of allele pairs
        for element in elements:
            a = int(enumerationofd[f1][int(element)])
            b = int(enumerationofd[f2][int(element)])
            c = int(allel1[int(element)+1])
            d = int(allel2[int(element)+1])
            if (c == 0 or d == 0):
                if a==0 and (b==c or b==d or b==0) is False:
                    return False
                elif b==0 and (a==c or a==d or a==0) is False:
                    return False
                elif (a!=0 or b!=0) and (a==c or a==d or b==c or b==d) is False:
                    return False
            elif (c != 0 and d !=0):
                if a==0 or b==0:
                    if a == 0 and (b==c or b==d or b==0) is False:
                        return False
                    elif b==0 and (a==c or a==d or a==0) is False:
                        return False
                elif a !=0 and b != 0:
                    if ((a+b==c+d) and (a==c or a==d or b==c or b==d)) is False:
                        return False
        return True

# from the input sequenced features(example:E1,E2,E3) extracts a list with integers for these features(1,3,5 for E1,E2,E3)
featurelist = []
for line in enumerationfile:
    line = line.strip().split('\t')
    featurelist = line[1:]
    break
enumerationfile.seek(0,0)
sortFlist = []
for feat in featurelist:
    sortFlist.append(featDic[feat])
featNumb = {}
for index,element in enumerate(sortFlist):
    featNumb[element] = int(index)
finalFeatures = []
for x in inputfeat:
    finalFeatures.append(str(featNumb[x]))

# finds the enumeration of the input alleles
allel1 = []
allel2 = []
for line in enumerationfile:
    line = line.strip().split('\t')
    if al1 == line[0]:
        allel1 = line[:]
    elif al2 == line[0]:
        allel2 = line[:]
enumerationfile.seek(0,0)

# evaluates the features selected to check for ambiguity. If all features selected are not sequenced it returns an error message
finalFeatures2 = finalFeatures.copy()
count1=0
count2=0
count3=0
for x in finalFeatures:
        a = int(allel1[(int(x)+1)])
        b = int(allel2[(int(x)+1)])
        if (a == 0 and b != 0) or ( a != 0 and b == 0):
                count1 += 1
        elif ((a == 0) and (b == 0)):
                count2 += 1
                finalFeatures2.remove(x)
        elif a != 0 and b != 0 :
                count3 += 1
evalFeatures=[]
if count1 == len(finalFeatures):
        print('too many not sequenced features. Choose different set of features')
        sys.exit()
elif count2 == len(finalFeatures):
        print('too many not sequenced features. Choose different set of features')
        sys.exit()
elif (count2 < len(finalFeatures) and count2 != 0):
        evalFeatures = finalFeatures2
elif count3 == len(finalFeatures):
        evalFeatures = finalFeatures
else:
    evalFeatures = finalFeatures

# creates lists with the features removed or remained from the previous block of code. Only to print them in the output file.
fremoved=[]
fremaining = []
for index,element in enumerate(sortFlist):
    if str(index) not in evalFeatures:
        fremoved.append(element)
    elif str(index) in evalFeatures:
        fremaining.append(element)

# Finds all acceptable alleles according the inputs alleles enumeration,creates all possible genotypes
# based on them and keeps only the acceptable genotypes that satisfy the criteria used.Finally creates the output file
listofalleles=[]
enumerationofd = {}
for line in itertools.islice(enumerationfile,1,None):
    line = line.strip().split('\t')
    if allele_list(evalFeatures) is True:
        listofalleles.append(line[0])
        enumerationofd[line[0]]=(line[1:])
comb = combinations(listofalleles,2)
for f1,f2 in comb:
    if possible_pairs(evalFeatures) is True:
            if f1!= al1 or f2!= al2:
                                outfile.write(al1+'+'+al2+'\t')
                                outfile.write(f1+'+'+ f2 +'\t')
                                outfile.write(','.join(fremaining)+'\t')
                                outfile.write('-'.join([allel1[int(x)+1] for x in evalFeatures])+'/')
                                outfile.write('-'.join([allel2[int(x)+1] for x in evalFeatures])+'\t')
                                outfile.write('-'.join([enumerationofd[f1][int(x)] for x in evalFeatures])+'/')
                                outfile.write('-'.join([enumerationofd[f2][int(x)] for x in evalFeatures]))
                                outfile.write('\n')
outfile.close()
