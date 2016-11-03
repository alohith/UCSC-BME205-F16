#!/usr/bin/env python3
########################################################################
# File:fastqPrep.py
#  executable: fastqPrep.py [format in (required)] -P33out[default] <inputfile >outputfile
# Purpose: read in fastq file from stdin and convert the quality string
#       to user specified format.
#
# Author:       Akshar Lohith
# Group(sounding boards): Brandon Saint-John, Lauren Sanders, Roger Volden
#                   Kishwar Shafin, Alison Tang
# History:      AL 09/26/2016 Created
#
#Program flow method
# 1) read in fastq file for parsing
#       Inputs: Fastq file/stdin
#       Outputs: generator that gives a tuple containing all read lines
# 2) Transform Q score
#       Inputs: User conversion options, plus quality string from parsed tuple
#       outputs: converted quality score
# 3) write new fastq file
#       Inputs: parsed header and sequence information from provided fastq file (stdin)
#           (aka first 3 elements from tuple), plus converted quality score
#       Outputs: write, in standard fastq format, the new fastq file with user specified quality scoring (to stdout)
# 4) In cases of wrong formatting, empty fastq records, or incomplete records, warn the user
#       how they are handled and handle them. Improper headers lead to program exit
########################################################################
import argparse
import sys
class FastQreader():
    '''
    Pulled from D.Bernick's BME160 lab4 files.
    modified for handling fastq files.


    Class to provide reading of a file containing one or more FASTQ
    formatted sequences:
    object instantiation:
    FastQreader(<file name>,<userOptions>):

    object attributes:
    fileIn: the initial file name
    usrOptionIn: user's input Option
    usrOptionOut: user's output Option
    recordCount: a count of Fastq records processed
    badSequenceCount: a count of bad fastq records and lines missing parts
    emptyLine: a count of empty lines found in the parsing of the fastq file
    totalLineCount: a count of the total lines in the input fastq file

    methods:
    readChunk() : Returns header, sequence,quality header and quality score as strings.
    giveChunks(): Return to stderr the metrics for the fastq file processed.
    Author: Akshar Lohith
    '''

    def __init__ (self, fileIn, usrOptions):
        '''contructor: saves attribute fileIn, user defined options, and initiates counters for # records processed, and # empty records (anything missing from a canonical Fastq record format)'''
        self.fileIn = fileIn
        self.usrOptionIn, self.usrOptionOut = usrOptions
        self.recordCount = 0
        self.badSequenceCount = 0
        self.emptyLine = 0
        self.totalLineCount = 0

    def giveChunks(self):
        """Return to the user metrics of how the Fastq file parsing."""
        print("Input Sequence processed with {} bad sequences or lines. FastQ records and lines have been skipped.".\
        format(self.badSequenceCount), file= sys.stderr)
        print("Input has {}/{} empty lines".format(self.emptyLine, self.totalLineCount), \
        file = sys.stderr)
        return "Processed {} complete FastQ records.".format(self.recordCount)

    def checkHeaderLineBadSequence(self, headerIsSeqLine):
        '''Check if header line is a sequence line'''

        isSequence = {True:0, False:0}
        for base in headerIsSeqLine:
            if base in set("ACTGUN"):
                isSequence[True] +=1
            else:
                isSequence[False] +=1

        if isSequence[True] > isSequence[False]:
            print("Sequence line found at line {} without header. Added to bad line count.".format(self.totalLineCount), file= sys.stderr)
            self.badSequenceCount +=1

        else:
            self.badSequenceCount +=1

    def readChunk (self):
        '''
        Using filename given to Class object, read through the fastq file and returns each
        included FastQ record as 4 strings - header, sequence, quality header and quality.
        Whitespace is removed, sequence contents are made upper, quality score is modified
        as user dictates.
        A counter is incremented to keep track of number of FastQ records, lines in the file,
        empty lines, and bad sequences.
        '''
        lineCounter = 0
        for line in self.fileIn:
            #Check for the fastq record state (header, sequence,quality header, quality sequence)?
            lineCounter += 1
            #Increment the counter for lines in the file.
            self.totalLineCount +=1

            if lineCounter == 1:
                if line.startswith('@'):
                    #Found a header line, check if ':' is in the first element
                    #and merge the first 2 elements if not with a ':' separator
                    if ':' not in line.rstrip().split(' ')[0]:
                        header = ':'.join(line.rstrip().split(' ')[0:2])
                    else: #Hold the headerline without the comment
                        header = line.rstrip().split(' ')[0]
                else: #catch odd cases where a nucleotide sequence is present instead of header
                    lineCounter = 0
                    self.checkHeaderLineBadSequence(line.rstrip())
                    continue

            elif lineCounter == 2:
                #Check if nucleotide sequence line is valid after replacing ambigious bases
                checkSeq = line.rstrip().upper().replace('*','N').replace('.','N')

                if checkSeq.isspace():
                    print("WARNING: No sequence associated with record {}, Record # {}. Record will be skipped.".\
                    format(header[1:], self.recordCount), file = sys.stderr)
                    self.badSequenceCount +=1
                    lineCounter = 0
                    continue;
                else:
                    sequence = checkSeq

            elif lineCounter ==3 and line.startswith('+'):
                qualHeaderLine = line.rstrip()
                 # Found a potential quality header line
                if qualHeaderLine.startswith('+') and len(qualHeaderLine) ==1 :
                    qualHead = qualHeaderLine
                elif qualHeaderLine.startswith('+') and len(qualHeaderLine) == len(header):
                    qualHead = qualHeaderLine[0] + header[1:]
                else: # store the header if it is of expected length
                    print("Expected quality header line (starting with '+'). Found: {}. Fastq record count up to this read is: {}.".\
                    format(qualHeaderLine,self.recordCount), file=sys.stderr)
                    self.badSequenceCount +=1
                    lineCounter = 0
                    continue;

            elif lineCounter == 4:
                quality = line.rstrip()
                lineCounter = 0

                if len(quality) == len(sequence):
                    self.recordCount += 1
                    conversion = PhredQuality(typeIn=self.usrOptionIn,\
                    typeOut=self.usrOptionOut, qualChunk=(header,sequence,qualHead,quality))
                    yield (conversion.evalTranslate())

                else:
                     #Increment bad sequence counter if the input file has messed up
                     #lengths for quality and nucleotide sequences and continue to next
                     #fastq record block
                    self.badSequenceCount += 1
                    print("Read Quality Score length does not match Sequence Read length for {}. Record # {}. Bad Sequence count is: {}.".\
                    format(header,self.recordCount+1, self.badSequenceCount), file = sys.stderr)
                    continue;

            elif line.isspace(): #check if a line is empty and increment the empty line counter
                self.emptyLine +=1
                print("WARNING: Empty line found at line: {}".format(self.totalLineCount), file= sys.stderr)
                continue;

            else: #crash and burn if really out of phase/a header line does not have a '@' charcter
                raise sys.exit("ERROR: '@' not in first column for '{}'. Line #: {} ".\
                format(line, self.totalLineCount))

    def fastaOut(self, header, sequence):
        '''Return to caller modified header line for fasta header type and the
            nucleotide sequence, unchanged '''

        return ('>'+header[1:],sequence)

class PhredQuality(object):
    """
        PhredQuality takes in the user's input and output options, passed to it
        from the fastq reader, as well as a tuple of 4 elements containing the
        header, sequence, quality header and quality sequence.

        This object has a main evaluate translate function and functions for each translation
        case. The evaluate translate function (which should be the only function called)
        evaluates the user options for input and output type and, after checking that the quality
        sequence is within the expected range of ASCII values for a defined input type, passes
        the quality sequence, and maybe the nucleotide sequence for translation processing.
    """

    def __init__(self, typeIn, typeOut, qualChunk):
        self.outputType = typeOut
        self.inputType = typeIn
        self.readHeader = qualChunk[0]
        self.qualHeader = qualChunk[2]
        self.inQualSeq = qualChunk[3]
        self.inReadSeq = qualChunk[1]

    def evalTranslate(self):
        '''Evaluations occur for P33in, P64in, P64SOLin, P64Bin
            followed by a check for quality encodings in user defined range,
            followed by evaluation of output types (P33out or P64out)
            before returning to call function, tuple with corrected quality sequence, and
            (if applicable) corrected nucleotide sequence.'''
        #Make a list of the quality scores from the encoding to confirm the user defined input format.
        qual = [ord(encode) for encode in list(self.inQualSeq)]

        #Evaluate the different cases of user defined input and fastq record type, and output to caller
        if self.inputType == '-P33in' or self.inputType == '--PHRED33input':
            if min(qual) not in range(33,98) or max(qual) not in range(33,98):
                raise Exception("ERROR: User specified PHRED-33 (Sanger/Illumina 1.8+) encoding, but file is not encoded as such!")
            elif self.outputType == '-P64out' or self.outputType == '--PHRED64output':
                return self.P33toP64(self.inQualSeq)

            elif self.outputType == '-P33out' or self.outputType == '--PHRED33output':
                return self.P33toP33(self.inQualSeq)

        elif self.inputType == '-P64in' or self.inputType == '--PHRED64input':
            if min(qual) not in range(64,129) or max(qual) not in range(64,129):
                raise Exception("ERROR: User specified PHRED-64 encoding (Illumina 1.3+), but file is not encoded as such!")
            elif self.outputType == '-P64out' or self.outputType == '--PHRED64output':
                return self.P64toP64(self.inQualSeq)

            elif self.outputType == '-P33out' or self.outputType == '--PHRED33output':
                return self.P64toP33(self.inQualSeq)

        elif self. inputType == '-P64Bin' or self.inputType == '--PHRED64-B-offset':
            if min(qual) not in range(66,129) or max(qual) not in range(66,129):
                raise Exception("ERROR: User specified PHRED-64-B-offset (Illumina 1.5+) encoding, but file is not encoded as such!")
            elif self.outputType == '-P64out' or self.outputType == '--PHRED64output':
                return self.P64BtoP64(self.inQualSeq, self.inReadSeq)

            elif self.outputType == '-P33out' or self.outputType == '--PHRED33output':
                return self.P64BtoP33(self.inQualSeq, self.inReadSeq)

        elif self.inputType == '-P64SOLin' or self.inputType == '--PHRED64-SolexaIn':
            if min(qual) not in range(59,127) or max(qual) not in range(59,127):
                raise Exception("ERROR: User specified SOLEXA-64 (Solexa) encoding, but file is not encoded as such!")
            elif self.outputType == '-P64out' or self.outputType == '--PHRED64output':
                return self.S64toP64(self.inQualSeq)

            elif self.outputType == '-P33out' or self.outputType == '--PHRED33output':
                return self.S64toP33(self.inQualSeq)

    def P33toP64(self, Qseq):
        '''Return tuple with corrected quality sequence line for P33 input and P64 output.'''

        translator = str.maketrans("".join([chr(x) for x in range(33,98)]),\
                    "".join([chr(y) for y in range(64,129)]))
        return (self.readHeader,self.inReadSeq, self.qualHeader,str(Qseq.translate(translator)))

    def P33toP33(self, Qseq):
        '''Return tuple with corrected quality sequence line for P33 input and P33 output. (No Changes)'''
        # no change to the quality score return it
        return (self.readHeader,self.inReadSeq, self.qualHeader,Qseq)

    def P64toP33(self, Qseq):
        '''Return tuple with corrected quality sequence line for P64 input and P33 output.'''

        translator = str.maketrans("".join([chr(x) for x in range(64,129)]),\
                     "".join([chr(y) for y in range(33,98)]))
        return (self.readHeader,self.inReadSeq, self.qualHeader,str(Qseq.translate(translator)))

    def P64toP64(self, Qseq):
        '''Return tuple with corrected quality sequence line for P64 input and P64 output. (No Changes)'''
        # no change to the quality score return it
        return (self.readHeader,self.inReadSeq, self.qualHeader,Qseq)

    def P64BtoP64(self, Qseq, readSeq):
        '''Return tuple with corrected quality sequence line for P64 with B offset input and P64 output.'''

        corReadSeq = list()
        corQseq = list()
        for base, score in zip(list(readSeq), list(Qseq)):
            #go through the nucleotide sequence and quality sequence at the same time
            #evaluate for min score for this encodingif score == 'B':
            if score == 'B':
                #floor the nucleotide sequence and quality score if min score is found
                corReadSeq.append('N')
                corQseq.append('@')
            else:
                #otherwise take the existing called base and quality score
                corReadSeq.append(base)
                corQseq.append(score)
        return (self.readHeader,''.join(corReadSeq), self.qualHeader, ''.join(corQseq))

    def P64BtoP33(self, Qseq, readSeq):
        '''Return tuple with corrected quality sequence line for P64 with B offset input and P33 output.'''

        corReadSeq = list()
        corQseq = list()
        for base, score in zip(list(readSeq), list(Qseq)):
            #go through the nucleotide sequence and quality sequence at the same time
            if score == 'B': #evaluate for min score for this encoding
                corReadSeq.append('N') #floor the nucleotide sequence if min score is found
                corQseq.append('@')#floor the quality score to P64out format
            else:
                corReadSeq.append(base)#otherwise take the existing called base
                corQseq.append(score)#otherwise take the existing quality score
        translator = str.maketrans("".join([chr(x) for x in range(64,129)]),\
                     "".join([chr(y) for y in range(33,98)]))
        return (self.readHeader,''.join(corReadSeq), self.qualHeader, ''.join(corQseq).translate(translator))

    def S64toP33(self, Qseq):
        '''Return tuple with corrected quality sequence line for Solexa 64 input and P33 output.'''
        from math import log10
        solexaScore =list() # create list to contain converted Solexa to PHred scores
        for code, score in zip(''.join([str(chr(x)) for x in range(59,127)]), range(-5,63)):
            #evaluate all possible solexa scores
            solexaScore.append(str(chr(round(10*log10(pow(10,(score/10))+1))+33))) #calculate the Phred score and add it to the conversion list
        translator = str.maketrans("".join([str(chr(x)) for x in range(59,127)]), "".join(solexaScore))
        return (self.readHeader,self.inReadSeq, self.qualHeader,Qseq.translate(translator))

    def S64toP64(self, Qseq):
        '''Return tuple with corrected quality sequence line for Solexa 64 input and P64 output.'''
        from math import log10
        solexaScore =list()# create list to contain converted Solexa to PHred scores
        for code, score in zip(''.join([str(chr(x)) for x in range(59,127)]), range(-5,63)):
            #evaluate all possible solexa scores
            solexaScore.append(str(chr(round(10*log10(pow(10,(score/10))+1))+64)))#calculate the Phred score and add it to the conversion list
        translator = str.maketrans("".join([str(chr(x)) for x in range(59,127)]), "".join(solexaScore))
        return self.readHeader,self.inReadSeq, self.qualHeader,str(Qseq.translate(translator))

####
# Use a CommandLine object to parse the user input from the command line
####
class CommandLine(object) :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']

    methods:
    do_usage_and_die()
    prints usage and help and terminates with an error.
    '''

    def __init__(self, inOpts=None) :
        '''Generate the parameter parser'''
        self.cmdparser = argparse.ArgumentParser(description="Handle user input for fastqPrep", add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s [format in (required)] -P33out[default] <inputfile >outputfile'
                                                  )
        self.groupIN = self.cmdparser.add_mutually_exclusive_group(required =True) #use a group for inputs so only 1 option is taken
        # This is a required option the user must provide
        self.groupOUT = self.cmdparser.add_mutually_exclusive_group() #use a group for outputs
        self.groupIN.add_argument('-P33in', '--PHRED33input', action = 'store', dest='usrInputType', nargs = '?', help='Defined user input fastq quality format PHRED-33 (Sanger format)', const='-P33in')
        self.groupIN.add_argument('-P64in', '--PHRED64input', action = 'store', dest='usrInputType', nargs = '?', help='Defined user input fastq quality format PHRED-64 (Illumina format)',const='-P64in')
        self.groupIN.add_argument('-P64Bin', '--PHRED64-B-offset', action = 'store', dest='usrInputType', nargs = '?', help='Defined user input fastq quality format PHRED-64 with B offset in Qscore', const='-P64Bin')
        self.groupIN.add_argument('-P64SOLin', '--PHRED64-SolexaIn', action = 'store', dest='usrInputType', nargs = '?', help='Defined user input fastq quality format SOLEXA-64 (prior Illumina 1.3+)', const='-P64SOLin')
        # output group argument addition
        self.groupOUT.add_argument('-P33out', '--PHRED33output', action = 'store', dest='usrWriteType', nargs ='?', help='Defined user output fastq quality format PHRED-33 (DEFAULT)', const ='-P33out', default ='-P33out')
        self.groupOUT.add_argument('-P64out', '--PHRED64output', action = 'store', dest='usrWriteType', nargs ='?', help='Defined user output fastq quality format PHRED-64', const='-P64out')
        #optional argument to output a fasta file instead of converted fastq
        self.cmdparser.add_argument('-faOut', '--fasta-out', action='store_true', dest='usrFastaOpt', help='Defined user output as fasta file format.')
        if inOpts is None :
            self.args = vars(self.cmdparser.parse_args()) # parse the CommandLine options
        else :
            self.args = vars(self.cmdparser.parse_args(inOpts)) # parse the input options
#####
# Create the main fucntion to handle the argument parsing from the user
#####
def main(myCommandLine=None):
    '''
    Implements the commandline options and writes to stdout.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(['-h'])

    try:
        usrOptions = myCommandLine.args
        #Instantiate the generator object to print out fastq reads
        fileInObj = FastQreader(sys.stdin,(usrOptions['usrInputType'],usrOptions['usrWriteType']))
        for record in fileInObj.readChunk():
            if usrOptions['usrFastaOpt']: # if user specifies to output a fasta file, fasta file will be generated instead of converted fastq
                print('\n'.join(fileInObj.fastaOut(record[0],record[1])))
            else:
                print(str(record[0]))
                print(str(record[1]))
                print(str(record[2]))
                print(str(record[3]))

        print(fileInObj.giveChunks(), file= sys.stderr)

    except Exception as errorInstance:
        print(errorInstance, file= sys.stderr)
        raise SystemExit
if __name__ == "__main__":
    main();
    raise SystemExit
