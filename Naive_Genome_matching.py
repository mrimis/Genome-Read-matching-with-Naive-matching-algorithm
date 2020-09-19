

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
    

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences
    
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def naive_with_rc(r,genome):
	matches_unique=[]
	matches=naive(r,genome)
	matches.extend(naive(reverseComplement(r), genome))
	for x in matches:
		if x not in matches_unique:
			matches_unique.append(x)
	return matches_unique

        	
p = 'AGGAGGTT'
filename= raw_input('Enter filename with path for genome fasta file")
t = readGenome(filename)
occurrences = naive(p,t)
occurrences = naive_with_rc(p, t)
print(occurrences)

def naive_with_mismatches(p, t, mismatches=2):#naive matching function for matching pattern to text with mismatches
    occurrences = []
    for i in range(0, len(t) - len(p) + 1): # for all alignments
        nmm = 0
        for j in range(0, len(p)):          # for all characters
            if t[i+j] != p[j]:               # does it match?
                nmm += 1                     # mismatch
                if nmm > mismatches:
                    break                    # exceeded maximum distance
        if nmm <= mismatches:
            # approximate match; return pair where first element is the
            # offset of the match and second is the Hamming distance
            occurrences.append(i)
    return occurrences
    
occurrences = naive_with_mismatches(p,t, mismatches=2)

def readFastq(filename):#Function for extracting Sequence and base qaulity from fastq file
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

seqs, quals = readFastq('ERR037900_1_first1000.fastq')

def phred33ToQ(qual):#function for converting phred33 scores to Quality scores
	return ord(qual) - 33
	
def createHist(qualities):#creating table for histogram with frequency of quality scores across reads position
	hist = [0]*50
	for qual in qualities:
		for phred in qual:
			q=phred33ToQ(phred)
			hist[q] +=1
	return hist

h=createHist(quals)



def average_quality(qualities):#function for finding average quality for each base across the read position
	hist=[0]*100
	totals=[0]*100
	
	for qual in qualities:
		for i in range(len(qual)):
			hist[i]+=phred33ToQ(qual[i])
			totals[i]+=1
	
	for i in range(len(hist)):
		hist[i] /=float(totals[i])
		
	return hist

hist = average_quality(quals)

    
    
    

	
			
	
