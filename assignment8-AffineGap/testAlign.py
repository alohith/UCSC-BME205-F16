#!/usr/bin/env python3
########################################################################
# File:testAlign.py
#  executable: testAlign.py
# Purpose:
#
# Author:       Akshar Lohith
# Group(sounding boards):
# History:      AL 11/14/2016 Created
#
########################################################################
import sys
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
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import random
        import argparse
        self.parser = argparse.ArgumentParser(description = 'This program will run a randomized motif search on a set of fasta records to find the motif with the lowest entropy score, as defined by Shannon\'s entropy formula', add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s [options] <input >output'
                                                  )
        self.parser.add_argument('-i', '--iterations', type=int, action = 'store', default=1000, dest='numIterate', help='Number of times to run the motif search before returning a consensus sequence. Default is 1000 iterations.')
        self.parser.add_argument('-p', '--pseudoCount', type=float, action = 'store', default=1, dest='definedPseudoCounts', help='Define a pseudo count to use for evaluation of the kmers in the sequence. Default is 1 pseudo count.')
        self.parser.add_argument('-k', '--kmerLength', type=int, action='store', default=13, dest='usrKlength',help='Set what kmer size to use to find the motif. This is a required parameter.')
        self.parser.add_argument('-r', '--randomShuffle', action='store_true', dest='shuffleOpt', help='Run the random motif search algorithim on shuffled sequences.')

        if inOpts is None :
            self.args = vars(self.parser.parse_args()) # parse the CommandLine options
        else :
            self.args = vars(self.parser.parse_args(inOpts)) # parse the input options


    def __del__ (self) :
        '''
        CommandLine destructor.
        '''
        # do something if needed to clean up before leaving
        pass

    def do_usage_and_die (self, str) :
        '''
        If a critical error is encountered, where it is suspected that the program is not being called with consistent parameters or data, this
        method will write out an error string (str), then terminate execution of the program.
        '''
        import sys
        print(str, file=sys.stderr)
        self.parser.print_usage()
        return 2

class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg
        
def main():
    pass

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    raise SystemExit
