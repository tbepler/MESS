import sys

from Bio import Entrez

#TODO return sequences in the form ( name, fold change, sequence )

def parseGeneFile( f ):
	# call parseGeneIds on each gene id in f 
    return parseGeneIds((line.strip() for line in f))

def setEmail(email):
    # set email
    Entrez.email = email

def parseGeneIds( geneIds ):
	# return generator that provides gene sequences in fasta format 
    setEmail("yhtgrace@gmail.com")
    for geneId in geneIds:
        handle = Entrez.efetch(db="nucleotide", id=geneId, retmode="fasta", rettype="text")
        # FIXME handle exceptions if id not found 
        yield (geneId, "".join(handle.read().split('\n')[1:]))

def parseFastaFile( f ):
	# return generator that provides each sequence
    # from http://stackoverflow.com/questions/7654971/
    name, seq = None, []
    for line in f:
        line = line.rstrip()
        if line[0] == '>':
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield(name, ''.join(seq))

def parse( f ):
	# parses the given file
    line = f.readline()
    f.seek(0)
    if line[0] == '>':
        return parseFastaFile(f)
    else:
        return parseGeneFile(f)

def main():
    test_gene = "../MESS_test/GSE29780_p0-05_down.gene_ids.txt"
    with open(test_gene) as f:
        it = parse(f)
        record = it.next()
        print >> sys.stderr, record
    test_fa = "../MESS_test/test.fa"
    with open(test_fa) as f:
        it = parse(f)
        record = it.next()
        print >> sys.stderr, record

if __name__ == '__main__':
    main()
