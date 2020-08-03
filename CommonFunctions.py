#!/usr/local/bin/python3
# Common Functions
# Ver. 0.02
# Pre.L

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import CommonFunctions as CF

#==========================="okFile"===========================#        
#Check if the given file is a file. 1 means yes, 0 means no.

def okFile(file):
    if(os.path.isfile(file)):
        return 1
    else:
        return 0

#==========================="getFile"===========================#
#Opens a file and parses it as a seqRecord file and returns the record.

def getFile(filename):
    try:
        OpenFastaFile = open(filename)    
    except IOError:
            # print "Error: cannot find {}".format(filename)
            print(f"Error: cannot find {filename}.")
    else:
        record = SeqIO.parse(OpenFastaFile,  "fasta")
        OpenFastaFile.close        
        return record
#==========================="getFile"===========================#

def createRecord(seq, identity):
    
    newRecord = SeqRecord(seq)
    newRecord.id = identity
    
    return newRecord
    

#========================"Self defined Errors and Exceptions"=======================#

class Error(Exception):
    """Base class for exceptions."""
    pass
    
#Error classes for handling input errors
class InputError (Error):
    """Exception for errors in the input.
        Attributes:
        expr -- the input that caused the error
        msg  -- message giving details of the error
    """
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg    
    def __str__(self):
        return repr(str(self.msg)+str(self.expr))

class RefSeqError(Error):
    def __init__(self, filename, reason):
        self.filename = filename
        self.reason = reason
    def __str__(self):
        return repr(str(self.reason) + str(self.filename))    
        

#========================"EpitopeInfo Structure"=======================#
class Epitope_Info:
    
    def __init__(self,allele,rank, start, end):  #parameter sequence omitted for now
        self.allele = allele
        #self.sequence = sequence
        self.rank = rank
        self.start = start
        self.end = end
    def HLA_type(self):
        return self.allele
    #def epitope(self):
    #    return self.sequence
    def percentile_rank(self):
        return self.rank
    def end_position(self):
        return self.end
    def start_position(self):
        return self.start

#==========================="Check if there are 1+ Ref sequences"===========================#

def recordCheck(record, refSeq_Filename):
    
    count = 0
    
    for item in record:
        if(count > 0):
            raise CF.RefSeqError(refSeq_Filename, "There is more than one sequence in ")
        count += 1



#========================"Main"=======================#

if __name__ == '__main__':
    print "This python script is used to store common functions as an import module."
        
