import sys

from itertools import izip

from Bio import Entrez

#TODO return sequences in the form ( name, fold change, sequence )

def parseGeneFile( f ):
    # call fetchSeq on each gene id in f, which has is in tab-delimited format
    # (gene_id, score)
    for line in f:
        geneId, score = line.strip().split()
        yield (geneId, float(score), fetchSeq(geneId))

def setEmail( email ):
    # set email
    Entrez.email = email

def fetchSeq( geneId ):
    handle = Entrez.efetch(db="nucleotide", id=geneId, retmode="fasta", rettype="text")
    return handle.read().partition('\n')[-1].replace('\n', '')

def parseGEO2R( f ):
    # parse GEO2R output
    for line in f:
        logFC, geneIds = line.strip().split()[5:7]
        logFC = float(logFC[1:-1])
        geneId = geneIds[1:-1].split(',')[0]
        if not geneId:
            continue
        yield (geneId, float(logFC), fetchSeq(geneId))

def parse( f ):
    # parses the given file
    line = f.readline()
    if len(line.strip().split()) == 2:
        f.seek(0)
        return parseGeneFile(f)
    elif line.strip().split()[0] == '"ID"':
        return parseGEO2R(f)
    else:
        return None # FIXME unknown file type

def main():
   
    setEmail("yhtgrace@gmail.com") # set email first

    test_2col = "../MESS_test/GSE29780_de.geneId-fc.txt"
    with open(test_2col) as f:
        it = parse(f)
        record = it.next()
        print >> sys.stderr, record

    test_geo2r = "../MESS_test/GSE29780_de.txt"
    with open(test_geo2r) as f:
        it = parse(f)
        record = it.next()
        print >> sys.stderr, record

if __name__ == '__main__':
    main()
