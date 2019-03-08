import itertools

filename = str(input("filename:"))

startkmers = []
endkmers = []
DNAbases = ["A", "T", "C", "G"]



while True:
    instartkmers = str(input("Start kmers (enter blank to finish):"))
    if instartkmers == "":
        break
    startkmers.append(instartkmers)

while True:
    inendkmers = str(input("End kmers (enter blank to finish):"))
    if inendkmers == "":
        break
    endkmers.append(inendkmers)

krange = []

start = {}
end = {}
other = {}

genetable = {}
genecounts = {}
genelengths = {}



def initks():
    for i, v in enumerate(startkmers):
        l = len(startkmers[i])
        if l not in krange:
            krange.append(l)
initks()

def initkmerdicts():
    for j, v in enumerate(startkmers):
        k = len(startkmers[j])
        if k not in start.keys():
            start[k] = []
        start[k].append(v)

    for i, v in enumerate(endkmers):
        l = len(endkmers[i])
        if l not in end.keys():
            end[l] = []
        end[l].append(v)

    for k in krange:
        genecounts[k] = []
        genelengths[k] = []
        other[k] = None
initkmerdicts()

def initgenetypes():
    for k in krange:
        genetable[k] = {}
        for s in start[k]:
            for e in end[k]:
                genetable[k][s+" "+e] = 0
initgenetypes()

def initother():
    for k in krange:
        all = [''.join(c) for c in itertools.product(DNAbases, repeat=k)]
        alls = []
        others = []
        for x in all:
            if x not in startkmers:
                alls.append(x)
        for i in alls:
            if i not in endkmers:
                others.append(i)
        other[k] = others
initother()

fasta = [">"]

def sequencefinder(filename):
    file = open(filename, 'r')
    igenome = "".join(line.strip() for line in file if line[0] not in fasta)
    return igenome

genome = sequencefinder(filename)


def geneparser(filename):
    genome = sequencefinder(filename)
    genomesize = len(genome)
    print(genome[0:21])
    print("Total genome length: "+str(genomesize))
    for k in krange:
        genes = 0
        for x in range(0, genomesize-k):
            kmer = genome[x:x+k]
            print("K "+kmer)
            if kmer in start[k]:
                for j in range(1, round(10000/k)):
                    print(j)
                    nkmer = genome[x+(j*k):x+((j+1)*k)]
                    print("N "+nkmer)
                    if nkmer in other[k]:
                        print("O "+str(genes))
                        continue
                    if nkmer in endkmers:
                        print("E "+str(genes))
                        genes += 1
                        genetable[k][kmer+" "+nkmer] += 1
                        break
                    if nkmer in startkmers:
                        print("Restart")
                        break
                    if nkmer not in other:
                        print("Not valid")
                        break
        genecounts[k] = genes

geneparser(filename)


print("Gene count table")
print(genetable)
print("Total gene counts")
print(genecounts)
