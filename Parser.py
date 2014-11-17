
# fasta sequence = ( name, sequence )

def parseGeneFile( f ):
	#TODO call parseGeneIds on each gene id in f
	raise NotImplementedError

def parseGeneIds( geneIds ):
	#TODO return generator that provides gene sequences in fasta format
	raise NotImplementedError

def parseFastaFile( f ):
	#TODO return generator that provides each sequence
	raise NotImplementedError

def getParser( f ):
	#TODO returns a parsing function appropriate for parsing the given file format
	raise NotImplementedError
