#!/home/preston/anaconda2/envs/binf2/bin/python
#Function Script
#Pre.L

import pysam, sys, os, re, glob, subprocess
import numpy as np
import scipy as sp
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.SeqIO.FastaIO import SimpleFastaParser
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import fcluster

#==========================="simpleParse"===========================#
# Opens a file and parses it. For parsing simple Fasta files this
# approach is much faster than using SeqIO.parse

def simpleParse(filename):
    readID = []
    readSeq = []
    try:
        OpenFastaFile = open(filename)
    except IOError:
            print("Error: cannot find {}".format(filename))
    else:
        for line in SimpleFastaParser(OpenFastaFile):
            readID.append(line[0])
            readSeq.append(line[1])
        OpenFastaFile.close()
        return readID, readSeq

#==========================="Find Insertions"===========================#
# seqList is the sequence excluding softclipping ####WORKS
# identifying insertions and their quality scores in each entry of a bamfile
# returns a dictionary with aligned sequence position as KEY
# and a length = 2 list containing nucleotide and quality score

def findInsertions(seqList, cigarList, qualList):    
    seqIndex = -1 #set at -1 to deal with stupid fence post issues 
    softClipCount = -1
    insertionDict = {}
    for cTuple in cigarList:        
        # is an insertion
        if(cTuple[0] == 1):
            numInsertions = cTuple[1]
            for i in range(0,numInsertions):
                seqIndex += 1
                insertionDict[seqIndex] = [seqList[seqIndex]]
        # Cigar shows a deletion in sequence (that is nt does not exist in seq)
        # 2 is deletion, 4 is soft padding, 6 is padded alignment (to other reads)
        elif(cTuple[0] == 2 or cTuple[0] == 4 or cTuple[0] == 6):
            # Don't move the seqIndex
            pass
        else:            
            seqIndex += cTuple[1]    
    for key in list(insertionDict.keys()):        
        insertionDict[key].append(list(qualList)[key])
    return insertionDict


#==========================="Find Deletions"===========================#
# Returns a dictionary with a tuple as key, containing aligned positions of where deletions occur
# Tuple length >1 indicates consecutive deletions
# Dictionary value stores what nucleotide found on ref but not on read

def findDeletions(alignedPosList, refSeq):
    
    deletionDict = {}
    refIndex = 0
    for posIndex in range(0,len(alignedPosList)):        
        if(posIndex == 0):
            refIndex = alignedPosList[posIndex]
            continue
        if(refIndex > len(refSeq)):
            break     
        if(alignedPosList[posIndex] - refIndex != 1):
            refIndexList = []
            refSeqItem = []
            for i in range(refIndex+1,alignedPosList[posIndex]):
                # i+1 to fix the fact that alignment starts at pos 1 rather than pos 0                
                refIndexList.append(i+1)
                refSeqItem.append(refSeq[i])
            deletionDict[tuple(refIndexList)] = refSeqItem 
        refIndex = alignedPosList[posIndex]
    return deletionDict
        
#delDict = findDeletions(testObj.get_reference_positions(), testObj.get_reference_sequence())
        
#==========================="Find Substitution"===========================#
# Returns a dictionary with read nt index as key of where substitions occur
# Dictionary value stores an array containing in order:
# Reference nt, ref index of nt, read nt, read qual 

def findSubstitutions(seqList, cigarList, qualList,alignedPosList, refSeq):
    seqIndex = -1 # set at -1 to deal with stupid fence post issues     
    refIndex = 0 
    softClipCount = -1
    subsDict = {}
    for cTuple in cigarList:        
        # is an match or mismatch
        if(cTuple[0] == 0):
            numMatch = cTuple[1]
            for i in range(0,numMatch):
            # check if it is match or a mismatch            
                seqIndex += 1
                if(str(seqList[seqIndex]).upper() != refSeq[alignedPosList[refIndex]].upper()):
                    # mismatch - store it
                    subsDict[seqIndex] = [refSeq[alignedPosList[refIndex]].upper(),
                                          alignedPosList[refIndex],
                                          seqList[seqIndex],
                                          list(qualList)[seqIndex]]
                else:
                    # match -> good. Nothing to be done.
                    pass
                refIndex += 1
        # Cigar shows a deletion in sequence (that is, nt does not exist in seq)
        # 2 is deletion, 4 is soft padding, 6 is padded alignment (to other reads)
        elif(cTuple[0] == 2 or cTuple[0] == 4 or cTuple[0] == 6):
            # Don't move read index
            pass
        else:            
            seqIndex += cTuple[1]
    return subsDict


#==========================="Find ReadStartCodon"===========================#
# Input seqList should be from read.query_alignment_sequence
# mappedToRef is from read.get_reference_positions()
# returns the expected start codon position of the read based on alignment with
# reference.
# Whether or not it actually codes M (ATG) is not the purpose of this function.

def readStartCodon(refStartCodon,seqList, cigarList, mappedToRef):
    readSC = refStartCodon    
    readIndex = 0
    prevMapVal = 0
    insCount = 0    
    fixedFencePost = False    
    for i in range(0,len(mappedToRef)):
        if(mappedToRef[i] > refStartCodon):            
            break
        if(i == 0):
            prevMapVal = mappedToRef[i]            
            continue
        if(mappedToRef[i] - prevMapVal != 1):            
            readSC -= (mappedToRef[i] - prevMapVal-1)
        prevMapVal = mappedToRef[i]    
    for cTuple in cigarList:
        if(mappedToRef[readIndex] > refStartCodon):
            break
        # is an insertion
        # cannot be seen in read.get_erference_positions()
        if(cTuple[0] == 1):            
            insCount += cTuple[1]
        # deletions already considered, ignore
        # ignore softclipping because read.query_alignment_sequence excludes soft clips    
        elif(cTuple[0] == 2 or cTuple[0] == 4):            
            pass
        elif(cTuple[0] == 0):            
            if(not fixedFencePost):                
                readIndex += (cTuple[1]-1)
                fixedFencePost = True
            else:                
                readIndex += cTuple[1]                        
    readSC += insCount
    return readSC


#==========================="findStopCodons"===========================#
# Takes in a read sequence. Looks within the sequence for stop
# codons. If one is identified, replace it with the triplet from 
# reference. Returns a hashtable (same structure) just with 
# stop codons removed. Return read sequence.

def findStopCodons(readSeq, refSeq, codingPos):
    stopCodons = ['TAA','TAG','TGA']
    codingIndex = codingPos - 1    
    firstFlag = True
    listSeq = list(readSeq)
    codonLocation = []
    codonString = ''        
    for i in range(codingIndex,len(listSeq)):
        if((i - codingIndex) % 3 == 0 and firstFlag == False):
            assert(len(codonString) == 3)
            assert(len(codonLocation) == 3)
            if(codonString in stopCodons):
                #replace it with the triplet from the same spot in the reference sequence
                listSeq[codonLocation[0]] = refSeq[codonLocation[0]]
                listSeq[codonLocation[1]] = refSeq[codonLocation[1]]
                listSeq[codonLocation[2]] = refSeq[codonLocation[2]]
            codonString = ''
            codonLocation = []        
        codonString = codonString + listSeq[i]
        codonLocation.append(i)        
        if(firstFlag == True):
            firstFlag = False        
    strSeq = ''.join(listSeq)
    return strSeq

#==========================="trimReads"===========================#
# This function is to trim reads to the shortest read found int the alignment
# so that there won't be reads having tail end with blanks when there are 
# longer reads present in the alignment.

def trimReads(readHash, trimLen):
    
    trimmedRead = ''
    newHash = {}
        
    for key in list(readHash.keys()):               
        if(len(key) > trimLen):            
            lenDiff = len(key) - trimLen
            trimmedRead = ''.join(list(key)[:-lenDiff])            
            
            if(trimmedRead in newHash):
                newHash[trimmedRead].append(readHash[key])
            else:
                newHash[trimmedRead] = readHash[key]
        else:
            if(key in newHash):
                newHash[key].append(readHash[key])
            else:
                newHash[key] = readHash[key] 
    return newHash


#==========================="do some subprocess"===========================#
# Splits the 0.5(n*(n-1)) hamming distance calculation into different instances
# so it can process the work in parallel.

def doSomeSubprocess(filePath, numReads, toolpath, jump = 2):
    
    procs = []
    recordJump = -1 # this is to remember current i+chunk, -1 because i starts at 0    
    subjectFaFile = filePath.split('/')[1]
    
    for i in range(0,numReads):        
        if(recordJump >= i):
            continue        
        recordJump = i+jump
        if(recordJump >= numReads):
            recordJump = numReads            
        assert(recordJump <= numReads)        
        proc = subprocess.Popen([sys.executable, toolpath + '/getHam.py', 
                                '{}'.format(i), '{}'.format(recordJump),
                                '{}'.format(numReads),'{}'.format(filePath)])        
        procs.append(proc)    
    for proc in procs:
        proc.wait()    
    subProcessCleanUp(filePath)
    print("Subprocessing Completed ({})\n".format(subjectFaFile))   

#==========================="clean up subprocess"===========================#
# the doSomeSubprocess function produces a lot of files. This function merges
# all the files into one .txt and then removes the rest.

def subProcessCleanUp(filePath):
    subjectFolder = filePath + "_HamDist"    
    try:
        os.makedirs('{}/tempFolder'.format(subjectFolder))        
    except OSError:        
        pass # already exists 
    mvCmd = ['mv'] + (glob.glob("{}/*_*.txt".format(subjectFolder)))+ ['{}/tempFolder'.format(subjectFolder)]
    subprocess.call(mvCmd)        
    # potentially dangerous since string is interpreted as shell cmd.
    # but should be ok because all file naming are controlled. 
    # users cannot control the naming of any of these files. Unless
    # they fiddle with the code.
    subprocess.call('cat  {}/tempFolder/*_*.txt > {}/All.txt'.format(subjectFolder,subjectFolder), shell=True) 
    subprocess.call(['rm', '-r','{}/tempFolder'.format(subjectFolder)])

#==========================="Get hamming distance"===========================#
# get some hamming distance given two strings

def getSomeHam(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))


#==========================="cleanString"===========================#
# Clean file names where characters like bracers, brackets and end
# and spaces. Space are made into under scores '_'

def cleanString(aList):
    newList = []
    for stringItem in aList:
        step1 = re.sub('\s$','',re.sub('[\*\[\]\(\)-]','',stringItem))
        step2 = re.sub('\s+','_',step1)
        newList.append(step2)
    return newList

#==========================="Parse edge data as np array"===========================#
# Parses the file which describes the edges in to a distance matrix as np array.
# Typically edges have information about the 2 nodes it is connected to.
# This function does not parse in information about edge direction.
# Refer to iGraph manual for python if directed edges are required.

def buildDM(hamDistFile, numReads):
    
    twoDArray = np.zeros((numReads,numReads))

    with open(hamDistFile,'r') as fileHandle:
        for line in fileHandle:
            r1 = line.split('\t')[0]
            r2 = line.split('\t')[1]
            r1Num = int(r1.split('HAP')[1])            
            r2Num = int(r2.split('HAP')[1])
                        
            if(r1Num >= numReads):
                print(f"This is r1Num: {r1Num}")
            if(int(r2Num) >= numReads):
                print(f"This is r2Num: {r2Num}, {type(r2Num)}, NumReads = {numReads}")

            assert(r1Num < numReads)
            assert(int(r2Num) < numReads)            
            twoDArray[r1Num][int(r2Num)] = int(line.split('\t')[2]) 
            twoDArray[int(r2Num)][r1Num] = int(line.split('\t')[2])             
        fileHandle.close()
    return twoDArray

#==========================="Perform Hierarchical Clustering"===========================#
# Takes in a distance matrix and returns clusters based on user given max_d
# if toDraw is turned on, it will make an image of the dendrogram to help
# users to determine what value to use for parameter max_d

def hierarchClust(dist_matrix, max_d, toDraw):
    Xcondensed = sp.spatial.distance.squareform(dist_matrix)
    Z = linkage(Xcondensed, method="complete")
    clusters = fcluster(Z, max_d, criterion='distance')    
    if(toDraw):
        fig = plt.figure(figsize=(25, 10)) 
        dendrogram(Z)
        plt.show()
    return clusters

#==========================="write clustered sequences into file"===========================#
# Exactly as the title says.

def printClusteredHap(cList, clusterID, faFile):
    clusterIndexes = [i for i, value in enumerate(cList) if value == clusterID]
    readID, readSeq = simpleParse(faFile)    
    results = faFile.split('/')[0]
    fName = faFile.split('/')[1].split('.')[0]+f"_Cluster{clusterID}.fa"
    clustFileName = results +'/Clusters/' + fName     
    try:
        os.makedirs(results+'/Clusters')        
    except OSError:        
        pass # already exists 
    
    assert(len(readID) == len(cList))

    with open(clustFileName, 'w') as writeHandle:
        for hapNum in clusterIndexes:
            writeHandle.write(f">{readID[hapNum]}\n{readSeq[hapNum]}\n")
        writeHandle.close()
    print(f"Writing {fName} to {results + '/Clusters/'} done!")
    return clustFileName
            

#===========================" create consensus sequence"===========================#
# This method uses alignIO to make a stupid consensus. The clusters given
# to this functions should already be relatively similar.

def makeConsensus(clusterFaFile, cThreshold):    
    cons = AlignInfo.SummaryInfo(AlignIO.read(clusterFaFile,"fasta")).dumb_consensus(threshold=cThreshold,ambiguous='N')
    if('X' in list(str(cons))):
        numX = list(str(cons)).count('N')
        print (f"Warning: ambiguous nucleotides found (represented as character 'N') in {clusterFaFile}.")
    return cons

#==========================="delete All.txt"===========================#
# Removes folder. Especially All.txt because it could take up a lot of 
# space. Once the code is done We don't really need that file anymore 
# unless the user specifies to retain Hamming Distance information, 
# which All.txt stores.

def removeFolder(filePath):    
    subprocess.call(['rm', '-r',f'{filePath}'])

#==========================="main function"===========================#
if(__name__ == '__main__'):   
    
    program_function = """
    *****
    
    Hellluuuuuu EveryNyaaaan!

      /\_/\
     {-o.o-}
      >   <    
    
    scribbles.py contains functions for indelRemover003 series.

    -Fluff

    *****
    """

    print(program_function)