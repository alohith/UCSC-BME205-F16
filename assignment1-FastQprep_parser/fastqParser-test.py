import sys
import argparse
class FastQreader():
    '''
    Pulled from D.Bernick's BME160 lab4 files.
    modified for handling fastq files.


    Class to provide reading of a file containing one or more FASTQ
    formatted sequences:
    object instantiation:
    FastQreader(<file name>):

    object attributes:
    fname: the initial file name

    methods:
    readChunk() : returns header, sequence,quality header and quality score as strings.
    Author: Akshar Lohith
    '''

    def __init__ (self, fname, usrOpts):
        '''contructor: saves attribute fname, user defined options, and initiates counters for # records processed, and # empty records (anything missing from a canonical Fastq record format)'''
        self.fname = fname
        self.usrOptIn, self.usrOptOut = usrOpts
        self.recordCount = 0
        self.badSequenceCount = 0
        self.emptyLine = 0
        self.totLineCount = 0
        self.badLineCount = 0
    def giveChunks(self):
        print("Input Sequence processed with {} bad sequences. These fastQ records have been skipped.".format(self.badSequenceCount), file= sys.stderr)
        print("Input has {} empty lines".format(self.emptyLine), file = sys.stderr)
        return print("Processed {} FastQ records.".format(self.recordCount),file=sys.stderr)
    def readChunk (self):
        '''
        using filename given in init, returns each included FastQ record
        as 4 strings - header, sequence, quality header and quality.
        whitespace is removed, sequence contents are made upper, quality score is untouched.
        A count is initiated to keep track of number of FastQ records
        '''
        fileH = self.fname
        # skip to first fasta header
        #line = fileH.readline()
        lineCounter = 1
        for line in fileH:
            if lineCounter ==1 :
                if line.startswith('@'):
                    if ':' not in line.rstrip().split(' ')[0]:
                        header = ':'.join(line.rstrip().split(' ')[0:2])
                    else:
                        header = line.rstrip().split(' ')[0]
                    #if line.rstrip().split(' ')[1]:
                        #print("Header comment found in read.
            if lineCounter == 2:
                checkSeq = line.rstrip().upper().replace('*','N').replace('.','N')
                if checkSeq.isspace():
                    print("WARNING: No sequence associated with record {}, Record # {}. Record will be skipped.".format(header[1:], self.recordCount), file = sys.stderr)
                    self.badSequenceCount +=1
                    lineCounter = 1
                    pass;
                else:
                    sequence = checkSeq
                    lineCounter = 3
            if lineCounter == 3:
                yo = line.rstrip()
                 # Some comment here.
                if yo.startswith('+') and len(yo) ==1 :
                    qualHead = yo
                    lineCounter = 4
                elif yo.startswith('+') and len(yo) == len(header):
                    qualHead = yo
                    lineCounter = 4
                else:
                    print("Expected quality header line (starting with '+'). Found: {}. Fastq record count up to this read is: {}. Record will be skipped.".format(yo, self.recordCount), file = sys.stderr)
                    self.badSequenceCount +=1
                    lineCounter = 1
                    #fileH.readline()
                    pass;
            if lineCounter == 4:
                quality = fileH.readline().rstrip()
                if len(quality) == len(sequence):
                    self.recordCount += 1
                    #print(header, sequence,qualHead,quality)
                    #conversion = PhredQuality(typeIn=self.usrOptIn,typeOut=self.usrOptOut, qualChunk=(header,sequence,qualHead,quality))
                    #(corHead, corSeq, corQhead, newQualseq) = conversion.evalTranslate()
                    #yield (conversion.evalTranslate())  #(str(corHead),str(corSeq),str(corQhead),str(newQualseq))
                    lineCounter = 1
                    yield (header,sequence,qualHead,quality)
                else: #crash and burn if the input file is messed up/lengths don't match
                    print("Read Quality Score length does not match Sequence Read length for {}. Record # {}: {}.".format(header,self.recordCount+1), file = sys.stderr)
                    lineCounter = 1
                    continue;
            elif line.isspace():
                print("WARNING: Empty line found at fastq record block {}".format(self.recordCount), file= sys.stderr)
                self.badSequenceCount +=1
                pass;
            else:
                #print("ERROR: '@' not in first column for '{}'.".format(line),file=sys.stderr)
                self.badSequenceCount +=1
                raise sys.exit("ERROR: '@' not in first column for '{}'.".format(line))
        # header is saved, get the rest of the sequence
        # up until the next header is found
        # then yield the results and wait for the next call.
        # next call will resume at the yield point
        # which is where we have the next header

        # final header and sequence will be seen with an end of file
        # return whatever is left..
        #self.chunkCount += 1
        #conversion = PhredQuality(typeIn=self.usrOptIn,typeOut=self.usrOptOut, qualChunk=(header,sequence,qualHead,quality))
        #(corHead, corSeq, corQhead, newQualseq) = conversion.evalTranslate()
        #yield (conversion.evalTranslate())#(str(corHead),str(corSeq),str(corQhead),str(newQualseq))

    def testGenerator(self):
        lineCounter = 0
        for line in sys.stdin:
            lineCounter += 1
            self.totLineCount +=1
            if lineCounter == 1 and line.startswith('@'):
                if ':' not in line.rstrip().split(' ')[0]:
                    header = ':'.join(line.rstrip().split(' ')[0:2])
                else:
                    header = line.rstrip().split(' ')[0]
                #print("header found!!")
            elif lineCounter == 2:
                checkSeq = line.rstrip().upper().replace('*','N').replace('.','N')
                if checkSeq.isspace():
                    print("WARNING: No sequence associated with record {}, Record # {}. Record will be skipped.".format(header[1:], self.recordCount), file = sys.stderr)
                    self.badSequenceCount +=1
                    lineCounter = 0
                    continue;
                else:
                    sequence = checkSeq
                #print("sequence found!!")
            elif lineCounter ==3 and line.startswith('+'):
                yo = line.rstrip()
                 # Some comment here.
                if yo.startswith('+') and len(yo) ==1 :
                    qualHead = yo
                    #print("quality header found!!")
                elif yo.startswith('+') and len(yo) == len(header):
                    qualHead = yo
                    #print("quality header found!!")
                else:
                    print("Expected quality header line (starting with '+'). Found: {}. Fastq record count up to this read is: {}.".format(yo, self.recordCount), file=sys.stderr)
                    self.badSequenceCount +=1
                    lineCounter = 0
                    sys.stdin.readline()
                    print("quality header NOT FOUND!!", file= sys.stderr)
                    #fileH.readline()
                    continue;
                #lineCounter =4
            elif lineCounter == 4:
                quality = line.rstrip()
                lineCounter = 0
                if len(quality) == len(sequence):
                    self.recordCount += 1
                    #print(header, sequence,qualHead,quality)
                    #conversion = PhredQuality(typeIn=self.usrOptIn,typeOut=self.usrOptOut, qualChunk=(header,sequence,qualHead,quality))
                    #(corHead, corSeq, corQhead, newQualseq) = conversion.evalTranslate()
                    #yield (conversion.evalTranslate())  #(str(corHead),str(corSeq),str(corQhead),str(newQualseq))
                    #print("Muhahahah!")
                    print (header)
                    print(sequence)
                    print(qualHead)
                    print(quality)
                else: #crash and burn if the input file is messed up/lengths don't match
                    self.badSequenceCount += 1
                    print("Read Quality Score length does not match Sequence Read length for {}. Record # {}. Bad Sequence count is: {}.".format(header,self.recordCount+1, self.badSequenceCount), file = sys.stderr)
                    continue;
            elif line.isspace():
                self.emptyLine +=1
                print("WARNING: Empty line found at fastq record block", file= sys.stderr)
                #self.badSequenceCount +=1
                continue;
            elif lineCounter==1:
                isSequence = {True:0, False:0}
                for base in line.rstrip():
                    if base in set("ACTGUN"):
                        isSequence[True] +=1
                    else:
                        isSequence[False] +=1
                if isSequence[True] > isSequence[False]:
                    print("Sequence line found at line {} without header. Added to bad line count.".format(self.totLineCount), file= sys.stderr)
                    self.badLineCount +=1
                    lineCounter = 0
                    continue;
                else:
                    lineCounter =0
                    continue
            else:
                raise sys.exit("ERROR: '@' not in first column for '{}'. Line #: {} ".format(line, self.totLineCount))

    def fastaOut(self, header, sequence):
        return ('>'+header[1:],sequence)

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
        #generate the parameter parser
        self.cmdparser = argparse.ArgumentParser(description="Handle user input for fastqPrep", add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s [format in (required)] -P33out[default] <inputfile >outputfile'
                                                  )
        self.groupIN = self.cmdparser.add_mutually_exclusive_group() #use a group for inputs so only 1 option is taken
        self.groupOUT = self.cmdparser.add_mutually_exclusive_group() #use a group for outputs
        self.groupIN.add_argument('-P33in', '--PHRED33input', action = 'store', dest='usrInputType', nargs = '?', help='Defined user input fastq quality format PHRED-33 (Sanger format)', const='-P33in')
        self.groupIN.add_argument('-P64in', '--PHRED64input', action = 'store', dest='usrInputType', nargs = '?', help='Defined user input fastq quality format PHRED-64 (Illumina format)',const='-P64in')
        self.groupIN.add_argument('-P64Bin', '--PHRED64-B-offset', action = 'store', dest='usrInputType', nargs = '?', help='Defined user input fastq quality format PHRED-64 with B offset in Qscore', const='-P64Bin')
        self.groupIN.add_argument('-P64SOLin', '--PHRED64-SolexaIn', action = 'store', dest='usrInputType', nargs = '?', help='Defined user input fastq quality format SOLEXA-64 (prior Illumina 1.3+)', const='-P64SOLin')

        self.groupOUT.add_argument('-P33out', '--PHRED33output', action = 'store', dest='usrWriteType', nargs ='?', help='Defined user output fastq quality format PHRED-33 (DEFAULT)', const ='-P33out', default ='-P33out')
        self.groupOUT.add_argument('-P64out', '--PHRED64output', action = 'store', dest='usrWriteType', nargs ='?', help='Defined user output fastq quality format PHRED-64', const='-P64out')
        if inOpts is None :
            self.args = vars(self.cmdparser.parse_args()) # parse the CommandLine options
        else :
            self.args = vars(self.cmdparser.parse_args(inOpts)) # parse the input options

class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg

def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(['-h'])


    #print(myCommandLine.args)  # print the parsed argument string .. as there is nothing better to do
    usrOptions = myCommandLine.args
    #print(usrOptions)
    #Calling Generator function to print out fastq reads
    #print(usrOptions['usrInputType'],usrOptions['usrWriteType'])
    fileInObj = FastQreader(sys.stdin,(usrOptions['usrInputType'],usrOptions['usrWriteType']))

    #for record in fileInObj.testGenerator():
    #    print(str(record[0]))
    #    print(str(record[1]))
    #    print(str(record[2]))
    #    print(str(record[3]))
    #print(fileInObj.giveChunks())
        #print(record)
        #print(p.evalTranslate())
    print("blah", file=sys.stderr)
    fileInObj.testGenerator()
    print(fileInObj.giveChunks(), file= sys.stderr)


        #raise SystemExit
if __name__ == "__main__":
    main();
    raise SystemExit
