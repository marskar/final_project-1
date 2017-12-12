#!/bin/usr/python 2.7.8

## Raw Stats Interpreter
## This script is designed to process the raw stats from sambamba
## Created by: Harper Fauni
## Created on: 9/5/2017

## Modification(s):

import argparse
import os

readCountList = []
noCovList = []
lowCovList = []
goodCovList = []
regionCoverageSummaryList = []

def statsProcessor(sambambaOut,outDir,runID):
    with open(sambambaOut, 'r') as input:
        sambambaLines = input.readlines()[1:]
        for line in sambambaLines:
            splitLine = line.strip('\n').split(' ')
            geneName = splitLine[3]
            transcriptName = splitLine[4]
            readCount = splitLine[6]
            meanCov = splitLine[7]
            zeroUpCovRegion = splitLine[8]
            goodCovRegion = splitLine[9]
            if float(meanCov) == 0:
                noCovList.append(splitLine)
            if float(meanCov) > 0 and float(meanCov) < 20:
                lowCovList.append(splitLine)
            if float(meanCov) >= 20:
                goodCovList.append(splitLine)
            readCountList.append(readCount)
            zeroCovPercentage = 100-float(zeroUpCovRegion)
            lowCovPercentage = 100-(float(zeroUpCovRegion)+float(goodCovRegion))

            regionCoverageSummaryList.append('{}\t{}\t{}\t{}\t{}\n'.format(geneName,transcriptName,zeroCovPercentage,lowCovPercentage,goodCovRegion))
            
    with open('{}/{}{}'.format(outDir,runID,'_Summary.txt'), 'w') as out:
        totalReadCount = str(sum(int(readCount) for readCount in readCountList))
        out.write('{}{}\n'.format('Total Read Count: ',totalReadCount))
        
    print "Total Read Count: "+str(sum(int(readCount) for readCount in readCountList))

    transcriptCoverageCounter(noCovList,lowCovList,goodCovList,outDir,runID)
    coverageBreakdownWriter(regionCoverageSummaryList,outDir,runID)
    
def transcriptCoverageCounter(noCovList,lowCovList,goodCovList,outDir,runID):
    noCovElements = len(noCovList)
    lowCovElements = len(lowCovList)
    goodCovElements = len(goodCovList)
    print "The number of transcript(s) that have no coverage: "+str(noCovElements)
    print "The number of transcript(s) that have low coverage: "+str(lowCovElements)
    print "The number of transcript(s) that have good coverage: "+str(goodCovElements)
    with open('{}/{}{}'.format(outDir,runID,'_Summary.txt'), 'a') as out:
        out.write('{}{}\n'.format('The number of transcript(s) that have no coverage: ',noCovElements))
        out.write('{}{}\n'.format('The number of transcript(s) that have low coverage: ',lowCovElements))
        out.write('{}{}\n'.format('The number of transcript(s) that have good coverage: ',goodCovElements))

def coverageBreakdownWriter(regionCoverageSummaryList,outDir,runID):    
    with open('{}/{}{}'.format(outDir,runID,'_Summary.txt'), 'a') as out:
        for regionCoverageSummary in regionCoverageSummaryList:
            out.write(regionCoverageSummary)

def baseCategorizer(sambambaOut,outDir,runID):

    totalBaseCounter=0
    noCovCounter=0 
    oneCovCounter=0 
    tenCovCounter=0 
    twentyCovCounter=0 
    thirtyCovCounter=0 

    with open(sambambaOut, 'r') as input:
        sambambaLines = input.readlines()[1:]
        for line in sambambaLines:
            totalBaseCounter+=1
            splitLine = line.strip('\n').split('\t')
            baseCov = int(splitLine[2])
            if baseCov == 0:
                noCovCounter+=1
            if baseCov >= 1:
                oneCovCounter+=1
            if baseCov >= 10:
                tenCovCounter+=1
            if baseCov >= 20:
                twentyCovCounter+=1
            if baseCov >= 30:
                thirtyCovCounter+=1

    noCovPercentage=float(noCovCounter)/float(totalBaseCounter)*100
    oneCovPercentage=float(oneCovCounter)/float(totalBaseCounter)*100
    tenCovPercentage=float(tenCovCounter)/float(totalBaseCounter)*100
    twentyCovPercentage=float(twentyCovCounter)/float(totalBaseCounter)*100
    thirtyCovPercentage=float(thirtyCovCounter)/float(totalBaseCounter)*100

    print "Percentage of bases with 0 coverage: "+str(noCovPercentage)+"%"
    print "Percentage of bases with coverage >= 1: "+str(oneCovPercentage)+"%"
    print "Percentage of bases with coverage >= 10: "+str(tenCovPercentage)+"%"
    print "Percentage of bases with coverage >= 20: "+str(twentyCovPercentage)+"%"
    print "Percentage of bases with coverage >= 30: "+str(thirtyCovPercentage)+"%"
    print "Total bases: "+str(totalBaseCounter)
    
    with open('{}/{}{}'.format(outDir,runID,'_Summary.tsv'), 'w') as out:
        out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('0 coverage (%)','>=1 coverage (%)','>=10 coverage (%)','>=20 coverage (%)','>=30 coverage (%)','Total bases'))
        out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(str(noCovPercentage),str(oneCovPercentage),str(tenCovPercentage),str(twentyCovPercentage),str(thirtyCovPercentage),str(totalBaseCounter)))

def depthDriver(sambambaOut,outDir,runID):
    print "Processing your stats data...\n"
    statsProcessor(sambambaOut,outDir,runID)

def baseDriver(sambambaOut,outDir,runID):
    print "Processing your sambamba depth base data...\n"
    baseCategorizer(sambambaOut,outDir,runID)
    


def main():

    # Sets the arguments using argparse

    argpar = argparse.ArgumentParser(
        description=(
            'This script is designed to process the raw stats from sambamba. this script is should produced a summarized version of the run.'
            'Assumes that this script is being ran on Biowulf'
        ),
        usage=(
            '\n-l\t< REQUIRED: Sambamba output >' 
            '\n-r\t< REQUIRED: Run Identifier, ex: ALS_COhort >' 
            '\n-o\t< OPTIONAL: Output pathway, default = current working directory >'
            '\n-b\t< OPTIONAL: Base flag, use this flag if the sambamba output is from "sambamba depth base" >'
        )
    )

    argpar.add_argument('-l', nargs='?', help='Sambamba output', required=True)
    argpar.add_argument('-r', nargs='?', help='Run identifier, ex: ALS_Cohort', required=True)
    argpar.add_argument('-o', nargs='?', help='Output pathway, default = current working directory', default=os.getcwd())
    argpar.add_argument('-b', help='Sbatch flag for big transcript list', action='store_true')

    allArgs = argpar.parse_args()

    # Parse args

    sambambaOut = allArgs.l
    runID = allArgs.r
    outDir = allArgs.o
    baseFlag = allArgs.b 

    if baseFlag:
        baseDriver(sambambaOut,outDir,runID)
    else:
        depthDriver(sambambaOut,outDir,runID)

if __name__ == "__main__":
    main()


















