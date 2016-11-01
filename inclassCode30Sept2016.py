for x,y in zip(a,b)
# a way to iterate through a pair of containers that can be zipped together!!
fname.open(path, mode)
#mode ={'r','w',etc...}
#generator expression
fred = (k*2 for k in myList)#genertor function?
print(fred)
print(list(fred))
print(list(fred))



header = "@HWI-ST611_0189:1:1101:1413:2137#0/1"
readSeq ="TTGCGATATTGATGGTGCCTACAGGGAGCTGAAGTCCCTCCCTGATGAATTACCCTGCCACTCGGGCACCAAGGACATCGCAAAGATCGGAAGAGCGGGT"
qualhead ="+HWI-ST611_0189:1:1101:1413:2137#0/1"
quality = "^^^cZ_acccccacd[ddeYY`Y^ddZbbY^cdaYbcdR^ccXXX^^accdbcdda_chY_\``cccSKZ__]_^^TRZ\[W\^BBBBBBBBBBBBBBBB"

corReadSeq = list()
for base, score in zip(readSeq, Qseq):
    if score == 'B':
        corReadSeq.append('N')
    else:
        corReadSeq.append(base)
