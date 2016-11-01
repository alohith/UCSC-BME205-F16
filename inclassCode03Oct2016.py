import sys
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
        self.chunkCount = 0
        self.badSequnceCount = 0
    def giveChunks(self):
        print("WARNING: Input Sequence processed with {} empty sequnces. These fastQ records have been skipped.".format(self.badSequnceCount), file= sys.stderr)
        return print("Processed {} FastQ records.".format(self.chunkCount),file=sys.stderr)
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
        for line in fileH:
            if line.startswith('@') :
                if ':' not in line.rstrip().split(' ')[0]:
                    header = ':'.join(line.rstrip().split(' ')[0:2])
                else:
                    header = line.rstrip().split(' ')[0]
                #if line.rstrip().split(' ')[1]:
                    #print("Header comment found in read.
                checkSeq =fileH.readline().rstrip().upper().replace('*','N').replace('.','N')
                if checkSeq.isspace():
                    print("WARNING: No sequence associated with record {}, Record # {}. Record will be skipped.".format(header[1:], self.chunkCount), file = sys.stderr)
                    self.badSequnceCount +=1
                    #fileH.readline()
                    #fileH.readline()
                    pass;
                else:
                    sequence = checkSeq
                yo = fileH.readline().rstrip()
                 # Some comment here.
                if yo.startswith('+') and len(yo) ==1 :
                    qualHead = yo
                elif yo.startswith('+') and len(yo) == len(header):
                    qualHead = yo
                else:
                    print("Expected quality header line (starting with '+'). Found: {}. Fastq record count up to this read is: {}. Record will be skipped.".format(yo, self.chunkCount), file = sys.stderr)
                    self.badSequnceCount +=1
                    #fileH.readline()
                    pass;
                quality = fileH.readline().rstrip()
                if len(quality) == len(sequence):
                    self.chunkCount += 1
                    #print(header, sequence,qualHead,quality)
                    conversion = PhredQuality(typeIn=self.usrOptIn,typeOut=self.usrOptOut, qualChunk=(header,sequence,qualHead,quality))
                    #(corHead, corSeq, corQhead, newQualseq) = conversion.evalTranslate()
                    yield (conversion.evalTranslate())  #(str(corHead),str(corSeq),str(corQhead),str(newQualseq))
                else: #crash and burn if the input file is messed up/lengths don't match
                    print("Read Quality Score length does not match Sequence Read length for {}. Record # {}: {}.".format(header,self.chunkCount+1), file = sys.stderr)
                    continue;
            elif line.isspace():
                print("WARNING: Empty line found at fastq record block {}".format(self.chunkCount), file= sys.stderr)
                self.badSequnceCount +=1
                pass;
            else:
                #print("ERROR: '@' not in first column for '{}'.".format(line),file=sys.stderr)
                raise sys.exit("ERROR: '@' not in first column for '{}'.".format(line))
                self.badSequnceCount +=1
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


def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(['-h'])

    try:

        #print(myCommandLine.args)  # print the parsed argument string .. as there is nothing better to do
        usrOptions = myCommandLine.args
        #print(usrOptions)
        #Calling Generator function to print out fastq reads
        #print(usrOptions['usrInputType'],usrOptions['usrWriteType'])
        fileInObj = FastQreader(sys.stdin,(usrOptions['usrInputType'],usrOptions['usrWriteType']))
        for record in fileInObj.readChunk():
            print(str(record[0]))
            print(str(record[1]))
            print(str(record[2]))
            print(str(record[3]))
        #print(fileInObj.giveChunks())
            #print(record)
            #print(p.evalTranslate())
        raise Usage('Error')

    except:
        print(fileInObj.giveChunks(), file= sys.stderr)
        #raise SystemExit
if __name__ == "__main__":
    main();
    raise SystemExit
