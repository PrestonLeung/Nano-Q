#!/usr/bin/python
#Ver. 0.03H
#Pre.L

import argparse, os, sys, re, pysam
import CommonFunctions as CF
import scribbles as sc

#========================"Global Variables"=======================#

OutputFileName = ''
StartCodon = 0
Jump = None
HD_Threshold = 234 # (Approximately 166+67.76, based on 1800 nt of HCV ONT sequencing)
ClusterMinRead = 30 # Conserative value based on error rate of ONT sequencing
ConsensusThreshold = 0.40
ShiftStart = 0

#==========================="iterate through alignment"===========================#
# This function is to go through the alignment and then 
# fiddle with insertions, deletions and snps of each read
# returns a dictionary with cleaned reads as well as the minimum length found 
# out of all the cleaned reads.

def iterAlign(alignment, refName, refLength, lenCutOff, qualThresh):    
    
    cleanReads = {}        
    global StartCodon
    minReadLen = lenCutOff
    subsDistribution = {}    
    rPos = 0
    tempReadCount = 0
    totalReadCount = 0
    for read in alignment.fetch(refName,1,refLength):        
        totalReadCount += 1
        # Checks the read length requirement to be accepted
        if(read.query_alignment_length < lenCutOff):
            tempReadCount += 1
            continue        
        readMappedPos = read.get_reference_positions()        
        # Makes sure that the read is starting at the rPos position
        if(readMappedPos[0] != rPos):            
            continue
        # Check read flag status        
        if(read.is_supplementary):
            print(f"Sup. Alignment: {str(read.query_name)}\tBit Flag:{read.flag}. Skipped.")
            continue
        if(read.is_secondary):
            print(f"Sec. Alignment: {str(read.query_name)}\tBit Flag:{read.flag}. Skipped.")
            continue
        # if(read.is_reverse):
        #     print(f"Reverse Complement: {str(read.query_name)}\t{read.flag}")
        
        refSeq = read.get_reference_sequence()
        readSeq = list(read.query_alignment_sequence)
        readCigars = read.cigartuples
        readQual = read.query_alignment_qualities
        errorCount = 0
        
        # grab insertions
        insDict = sc.findInsertions(readSeq, readCigars, readQual)        
        # grab deletions                
        delDict = sc.findDeletions(readMappedPos, refSeq)
        # grab substitutions
        subsDict = sc.findSubstitutions(readSeq,readCigars,readQual,readMappedPos,refSeq)
        # grab the "expected" starting codon position based on alignment
        readSCPos = sc.readStartCodon(StartCodon,readSeq,readCigars,readMappedPos)

        # Need to fiddle with the insertions before read start 
        # codon. If the start codon is the first position in the
        # given alignment, this short code essentially does nothing.        
        insCount = 0
        for key in sorted(insDict.keys()):            
            if(key > readSCPos):                
                break           
            readSeq[key] = '&'
            insCount += 1 
        readSCPos_copy = readSCPos - insCount 
        
        # Check substitions first because it doesn't change read length
        # If quality is below qualThresh, we replace it with ref
        for key in sorted(subsDict.keys()):
            if(subsDict[key][3] < qualThresh):
                readSeq[key] = subsDict[key][0]
                errorCount += 1

        # Count the number of snps that is below quality theshold in one read. 
        # These are considered as the number of "errors" per read.
        # Store it inside a distribution table.
        
        if(errorCount in subsDistribution):
            subsDistribution[errorCount] += 1
        else:        
            subsDistribution[errorCount] = 1

        # deal with insertions after the read start codon
        # note key is in index form and readSCPos requires -1
        # to match with index        
        for key in sorted(insDict.keys()):
            if(key <= readSCPos_copy):
                continue
            # Removing insertions indiscriminately
            # Does not care about quality score
            readSeq[key] = '&'
        
        # remove all occurrences of '&' 
        strlist = list(''.join([a for a in readSeq if a != '&']))

        # check out the deletions. For all dels 
        # insert ref nucleotide into the sequence        
        for key in sorted(delDict.keys()):            
            for i in range(0,len(key)):
                strlist.insert(key[i]-1, delDict[key][i])
       
        strKey = ''.join(strlist)        
        strKey2 = sc.findStopCodons(strKey, refSeq,StartCodon)
         
        # This values looks for the read that's the shortest after 
        # cleaning / trimming (e.g insertions removed). This value
        # is used to make all reads align with the same length within
        # it's reference scope. If there are multiple references, for
        # each reference used, the minReadLen may be different for 
        # that particular alignment
        if(len(strKey2) < minReadLen):
            minReadLen = len(strKey2)

        if(strKey2 in cleanReads):
            cleanReads[strKey2].append(read.query_name)
        else:
            cleanReads[strKey2] = [read.query_name]
    print(f"Total reads skipped = {tempReadCount}")
    print(f"Total reads looked at = {totalReadCount}")
    return cleanReads, minReadLen, subsDistribution

#==========================="main function"===========================#

if __name__ == '__main__':   
    
    program_function = """
    *****
    Nano-Q Tool:
        
        This tool takes in a bam file performs a conservative cleaning procedure and
        then uses a hierachical clustering method based on hamming distance to identify
        potential variants. For each cluster formed using hierachical clustering, a 
        consensus sequence is produced to represent that cluster and frequency of 
        occurrence is measured based on reads per cluster over total reads extracted
        from the alignment.

        Process of algorithm:            
            1) For SNPS, it will check the base quality scores against user specified
               cutoff and substitute it with reference sequence used for the alignment
               in the bamfile if it is below the cutoff.
            2) It will remove all indels and substitutes it to the reference sequence
               used for the alignment in the bamfile.
            3) Stop codons check after cleaning is also performed. If a stop codon is 
               identified the entire triplet of nucleotide is replaced with the 
               reference used for the alignment. **Note - coding region is controlled 
               by option '-c'.
            4) The reads are then trimmed to the read with the shortest length to 
               keep output read lengths identical.
            5) Pairwise hamming distance are then considered for all reads, which
               are then used in a distance matrix for hierarchical clustering.
        

    *****
    """   
    parser = argparse.ArgumentParser()
    
    # Required Arguments
    parser.add_argument('-b', '--bamfile', help='Filename of bam file.',  required = True)
    parser.add_argument('-c', '--code_start', help='Start codon position in the reference sequence',required = True, type = int)
    parser.add_argument('-l','--read_length', help="Length cut off for read size", required = True, type = int)
    parser.add_argument('-nr', '--num_ref', help="'Number of references used in the alignment.", required = True, type = int)
    parser.add_argument('-q', '--qual_threshold', help = "Base quality score cut off.", required = True, type = int)
    
    
    # Optional Arguments
    parser.add_argument('-j','--jump', help="Increase this to make larger read intervals.", type = int)
    parser.add_argument('-ht','--HD_Threshold', help="Hamming distance threshold used to call clusters. [Default = 234]", type = float)
    parser.add_argument('-mc','--MinRead_Cluster', help="Minimum no. of reads to accept as a cluster. [Default = 30]", type = int)
    parser.add_argument('-ct','--Consensus_Threshold', help="Threshold to call a nucleotide to be consensus. [Default = 0.40]", type = float)
    parser.add_argument('-d','--dendrogram', help="Draw a dendrogram to help determine HD_Threshold (-ht) cutoff.", action='store_true')
    parser.add_argument('-hd','--keep_hdFile', help="Retain the Hamming Distances calculated. **Note: Could take up a lot of space**", action='store_true')
    parser.add_argument('-kc','--keep_clusters', help="Retain the clustered reads. **Note: Could take up a lot of space**", action='store_true')
    # parser.add_argument('-ss','--shift_start', help="Shifts the starting position of the reads aligned to reference if reads don't span reference. [Default = 1]", type = int)
    
    if(len(sys.argv) < 2):        
        print(program_function)
        parser.print_help()
        sys.exit()
    args=parser.parse_args()    
    
    if(not CF.okFile(args.bamfile)):
        raise CF.InputError(args.bamfile,  "Invalid bam file: ")
    if(args.jump):
        Jump = args.jump
    else:
        Jump = 10
    if(args.HD_Threshold):
        HD_Threshold = args.HD_Threshold
    if(args.MinRead_Cluster):
        ClusterMinRead = args.MinRead_Cluster
    if(args.Consensus_Threshold):
        ConsensusThreshold = args.Consensus_Threshold
    # if(args.shift_start):
    #     ShiftStart = (args.shift_start-1)
    StartCodon = args.code_start #set to global var because this doesn't change
    alignment = pysam.AlignmentFile(args.bamfile, "rb")
    refName = [alignment.get_reference_name(i) for i in range(0,args.num_ref)]
    refNameCleaned = sc.cleanString([alignment.get_reference_name(i) for i in range(0,args.num_ref)])
    refLen = [alignment.get_reference_length(refID) for refID in refName] 
    
    print("###########################################")
    print("# Stage 1.0 - Read cleaning and selection #")
    print("###########################################")
    
    fileList = []
    numCleanedList = []
    testTotalCleanReads = 0
    try:
        os.makedirs("Results")
    except OSError:
        pass # already exists
    
    for refIndex in range(0,args.num_ref):        
        cleaned, minReadLen, subsDistribution = iterAlign(alignment, refName[refIndex], refLen[refIndex], args.read_length, args.qual_threshold)
        trimCleaned = sc.trimReads(cleaned, minReadLen)
        
        if(len(trimCleaned) > 1 ):
            print(f"{refNameCleaned[refIndex]}: Number of trimmed Reads = {len(trimCleaned)}.")
            numCleanedList.append(len(trimCleaned))
            # testTotalCleanReads += len(trimCleaned)
            tempFileName = "Results/"+re.sub('\/','',refNameCleaned[refIndex])+".fa"
            fileList.append(tempFileName)
            hapNum = 0
            openFile = open(tempFileName, 'w')    
            
            for key in trimCleaned:
                hapName = 'testHAP' + str(hapNum)
                openFile.write(f">{hapName}\n{key}\n")
                hapNum += 1
            openFile.close()
        else:
            print(f"{refNameCleaned[refIndex]}: Number of trimmed Reads = {len(trimCleaned)}. This data set is excluded.")
        testTotalCleanReads += len(trimCleaned)
    print(f"Total number of accepted cleaned reads: {sum(numCleanedList)} and total cleaned reads: {testTotalCleanReads}.")
    print("Stage 1.0 Complete.\n")

    print("####################################################")
    print("# Stage 2.0 - Calculate Pairwise Hamming Distances #")
    print("####################################################")
    
    assert((len(numCleanedList)+len(fileList)) / 2 == len(fileList))    

    for i in range (0,len(numCleanedList)):
        print(f"Subprocessing {refNameCleaned[i]}.")

        if(numCleanedList[i] > 0):
            toolPath = '/'.join(re.split('/', sys.argv[0])[:-1])
            if(Jump):
                sc.doSomeSubprocess(fileList[i],numCleanedList[i], toolPath, jump = Jump)
            else:
                sc.doSomeSubprocess(fileList[i],numCleanedList[i], toolPath)
        else:            
            print(f"{refNameCleaned[i]} has {numCleanedList[i]} reads. Skipped.\n")
            
    print("Stage 2.0 Complete.\n")

    print("############################################")
    print("# Stage 3.0 - Constructing distance matrix #")
    print("############################################")

    dmList = []
    for i in range(0,len(fileList)):
        fName = re.split("/", fileList[i])[1]
        print(f"\nConstructing distance matrix for {fName}. There are {numCleanedList[i]} reads.")
        hamFilePath = fileList[i] + "_HamDist"
        hamFile = fileList[i] + "_HamDist/All.txt"                
        
        if(numCleanedList[i] > 1):
            dmList.append(sc.buildDM(hamFile,numCleanedList[i]))
        else:            
            print("Skipped. \n")
        
        if(not args.keep_hdFile):
            sc.removeFolder(hamFilePath) 
        
    print("Stage 3.0 Complete.\n")
    
    print("##################################################")
    print("# Stage 4.0 - Performing Hierarchical Clustering #")
    print("##################################################")
    
    # using index here because gives us the same index for using fileList
    clusterFiles = []
    clusterReadCount = []
    clusterTotal = []
    tempDict = {}
    clusterGrid = None
    for i in range(0,len(dmList)):
        clusterPath = []
        readCount = []
        dist_matrix = dmList[i]        
        print(f"Analysing {fileList[i].split('/')[1].split('.')[0]}...")
        clusterGrid = sc.hierarchClust(dist_matrix,HD_Threshold, args.dendrogram)
        print(f"Total number of clusters identified: {len(set(clusterGrid))}")
        print(f"Retrieving clusters with more than {ClusterMinRead} reads...")
        retainedReads = 0

        for clusterID in clusterGrid:
            if(clusterID not in tempDict):
                tempDict[clusterID] = 1
            else:
                tempDict[clusterID] += 1

        for key in sorted(tempDict.keys()):
            if(tempDict[key] > ClusterMinRead):                
                print(f"ClusterID {key} has {tempDict[key]} reads.")
                clusterPath.append(sc.printClusteredHap(clusterGrid,key,fileList[i]))
                readCount.append(tempDict[key])
                retainedReads += tempDict[key]
        print(f"Remaining clusters = {len(clusterPath)} and a total of {retainedReads} reads.")
        print("--------------------------------------\n")
        clusterFiles.append(clusterPath)
        clusterReadCount.append(readCount)
        clusterTotal.append(retainedReads)
        tempDict.clear() # This line empties out the dictionary

    print("Stage 4.0 Complete.\n")
    del tempDict, clusterGrid, dmList

    print("############################################################")
    print("# Stage 5.0 - Calcuating consensus sequences from clusters #")
    print("############################################################")
    print("Creating consensus sequences for each cluster...\n")
    
    allCons = []
    
    for i in range(0,len(clusterFiles)):
        consDict = {}
        for j in range(0,len(clusterFiles[i])):
            freq = round(clusterReadCount[i][j]/float(clusterTotal[i]),4)
            conSeq = sc.makeConsensus(clusterFiles[i][j],ConsensusThreshold)            
            if(str(conSeq) not in consDict):                
                consDict[str(conSeq)] = freq
            else:
                consDict[str(conSeq)] += freq
        allCons.append(consDict)
        
    for i in range(0,len(allCons)):
        if(len(allCons[i]) < 1):
            continue
        refID = fileList[i].split('/')[1].split('.')[0]
        print("Writing to file....")
        with open(f"Results/{refID}_ClusterConsensus.fa", 'w') as writeHandle:
            hapCount = 0
            for seqKey in allCons[i].keys():
                clusterID = f"consHap{hapCount}_{allCons[i][seqKey]}"
                writeHandle.write(f">{clusterID}\n{seqKey}\n")
                hapCount += 1
        writeHandle.close()
    
    print("Finished writing.")
    
    # Re open the file and find all the sequences that are identical. We then
    # put them together and add the already normalised frequencies
    for i in range(0,len(allCons)):
        if(len(allCons[i]) < 1):
            continue
        refID = fileList[i].split('/')[1].split('.')[0]
        print ("Consolidating consensus.....")
        seqHash = sc.consolidateConsensus(f"Results/{refID}_ClusterConsensus.fa")        
        with open(f"Results/{refID}_ClusterConsensusFinal.fa", 'w') as writeHandle:
            hapCount = 0
            for seqKey in seqHash.keys():                
                readName = f"consHap{hapCount}_{seqHash[seqKey]}"
                writeHandle.write(f">{readName}\n{seqKey}\n")
                hapCount += 1
        writeHandle.close()
        sc.removeFile(f"Results/{refID}_ClusterConsensus.fa")
    
    print("Finished writing.")
    print("Removing unnecessary files...")
    
    if(not args.keep_clusters and os.path.exists('Results/Clusters')):
        sc.removeFolder('Results/Clusters')
    
    print("Done!\n")
    print("Stage 5.0 Complete.\n")

# python indelRemover003H.py -b BC03_unique_sort.bam -c 1 -l 1800 -nr 5 -q 5 -j 10          

