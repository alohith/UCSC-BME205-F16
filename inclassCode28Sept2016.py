#!bin/
class switcher():
    def __init__(self, pick = True):
        if pick:
            self.writer = self.fred
        else:
            self.writer = self.mary
    def writer (self):
        pass
    def fred(self):
        print("fred")
    def mary(self):
        print("mary")

obj = switcher(True)
obj.writer()
obj = switcher(False)
obj.writer()

obj = switcher()
obj.writer()

''''
prints:
fred
mary
fred
'''


args = cmdparser.parse_args()
#print(args.usrInputType)
#print(args.usrWriteType)
print(args)
# User input type check
if args.usrInputType == '-P33in' or '--PHRED33input':
    #Pass this to generator?
    print("This")
    pass
elif args.usrInputType == '-P64in' or '--PHRED64input':
    print("IS")
    pass
elif args.usrInputType == '-P64Bin' or '--PHRED64-B-offset':
    print("Very")
    pass
elif args.usrInputType == '-P64SOLin' or '--PHRED64-SolexaIn':
    print("frustrating")
    pass


# user output type
if args.usrWriteType == '-P64out' or '--PHRED64output':
    print('Crazy')
    pass
elif args.usrWriteType == '-P33out'or '--PHRED33output':
    print("shit")
    #usrWriteType = '-P33out'
#Testing Generator function to print out fastq reads
testObj = FastQreader(sys.stdin)
for blah in testObj.readFastqChunk():
    print(blah)
#    print("yo")
