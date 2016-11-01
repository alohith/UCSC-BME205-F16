class FastAreader :
    '''
    Class to provide reading of a file containing one or more FASTA
    formatted sequences:
    object instantiation:
    FastAreader(<file name>):

    object attributes:
    fname: the initial file name

    methods:
    readFasta() : returns header and sequence as strings.
    Author: David Bernick
    Date: April 19, 2013
    '''
    def __init__ (self, fname):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def readFasta (self):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''

        with open(self.fname) as fileH:
            # initialize return containers
            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('@') :
                print("hello")
                line = fileH.readline()
                print(line)
            header = line[1:].rstrip()
            print("I has header")

            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            count = 0
            for line in fileH:
                print(count)
                if line.startswith ('@'):
                    print("count")
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                    count += 1
                else :
                    print("yo")
                    sequence += ''.join(line.rstrip().split()).upper()
                    count -= 1
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence

def main():
    myReader = FastAreader('biopython_2.fq')
    for head, seq in myReader.readFasta():
        print(head, seq)
if __name__ == "__main__":
    main();
    raise SystemExit
