#!/usr/bin/env python3
########################################################################
# File:program.py
#  executable: program.py
# Purpose:scaffold file for python executables

#   stderr: errors and status
#   stdout:
#
# Author: David Bernick, Verena Friedl
# History:      dlb 08/20/2011 Created
#               vf 09/24/2016 updated for python3.5.2
#
########################################################################

########################################################################
# CommandLine
########################################################################
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

        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does',
                                             epilog = 'Program epilog - some other stuff you feel compelled to say',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                                  )
        self.parser.add_argument('-o', '--option', action = 'store', help='foo help')
        self.parser.add_argument('-i', '--integer', type=int, choices=range(5, 10), action = 'store', help='help for a boolean option')
        self.parser.add_argument('-c', '--character', choices ='abcdef', action = 'store', help='help for a charcter option')
        self.parser.add_argument('-b', '--bool', action = 'store', nargs='?', const=True, default=False, help='boolean switch')
        self.parser.add_argument('-r', '--requiredBool', action = 'store', nargs='?', required=True,const=True, default=False, help='required boolean switch')
        self.parser.add_argument('-l', '--list', action = 'append', nargs='?', help='list help') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
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
########################################################################
# Main
# Here is the main program
#
#
########################################################################


def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(['-r'])

    try:

        print(myCommandLine.args)  # print the parsed argument string .. as there is nothing better to do

        if (myCommandLine.args['requiredBool']) :
            print('requiredBool is', str(myCommandLine.args['requiredBool']))
        else :
            pass
#        raise Usage('testing')

    except Usage as err:
       myCommandLine.do_usage_and_die(err.msg)

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    raise SystemExit
