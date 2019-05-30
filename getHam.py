#!/home/preston/anaconda2/envs/binf2/bin/python
#Pre.L

from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import scribbles as sc
import argparse, sys, os

#==========================="main function"===========================#

if __name__ == '__main__':   

    program_function = """
    *****
    getHam.


    *****
    """   

    parser = argparse.ArgumentParser()
    
    #Required Arguments
    parser.add_argument('start', type = int)
    parser.add_argument('end', type = int)
    parser.add_argument('numReads',type=int)
    parser.add_argument('faFile',type = str) 
    
    args = parser.parse_args()    

    if len(sys.argv) < 2:
        print(program_function)
        parser.print_help()
        sys.exit()
    
    recordSeqList = []
    recordNameList = []
    # partner_i = []
    # partner_j = []
    # hamList = []    
    # seqRecord = sc.getFile(args.faFile)

    recordNameList, recordSeqList = sc.simpleParse(args.faFile)    

    hamPath = args.faFile + "_HamDist"
    try:
        os.makedirs(hamPath)
        print("Creating directory: {}".format(hamPath))
    except OSError:
        #print "Using existing directory: {}".format(hamPath)
        pass # already exists
    
    # for string in seqRecord    :
    #     recordNameList.append(string.id)
    #     recordSeqList.append(string.seq)
        
    loopName = "{}_{}.txt".format(args.start,args.end)
    writeHandle = open(hamPath+"/"+loopName,'w')
    # assert((len(hamList)+len(partner_i)+len(partner_j)) / 3 == len(hamList))

    for i in range(args.start,args.end+1):      
        assert(i < args.end+1)
        for j in range(i+1, args.numReads): 
            assert(j < args.numReads)            
            writeHandle.write(str(recordNameList[i]) +"\t")
            writeHandle.write(str(recordNameList[j]) +"\t")
            writeHandle.write(str(sc.getSomeHam(str(recordSeqList[i]), str(recordSeqList[j]))) +"\n")


            # hamList.append(sc.getSomeHam(str(recordSeqList[i]), str(recordSeqList[j])))
            # partner_i.append(recordNameList[i])
            # partner_j.append(recordNameList[j])

    # for i in range(0,len(hamList)):
    #     writeHandle.write(str(partner_i[i]) +"\t")
    #     writeHandle.write(str(partner_j[i]) +"\t")
    #     writeHandle.write(str(hamList[i]) +"\n")
    
    writeHandle.close()
        
