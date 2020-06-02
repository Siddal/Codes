# From 3' to 5'
DNANucleotides = ["A", "C", "G", "T"]
DNA_Reverse = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
RNANucleotides = ["A", "C", "G", "U"]
DNAtoRNA = {}

# Validate DNA sequence
# valiDNASeq()
def valiDNASeq(dna_seq):
    tmpDNAseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in DNANucleotides:
            return print("Errors in DNA sequence")
    return tmpDNAseq

# Validate RNA sequence
# valiRNASeq()
def valiRNASeq(rna_seq):
    tmpRNAseq = rna_seq.upper()
    for nuc in tmpseq:
        if nuc not in RNANucleotides:
            return print("Errors in RNA sequence")
    return tmpRNAseq

# Count DNA nucleotides
# countDNANucFrequency()
def countDNANucFrequency(dna_seq):
    tmpDNAFreqDict = {"A":0, "C":0, "G":0, "T":0}
    for nuc in dna_seq:
        tmpDNAFreqDict[nuc] += 1
    return tmpDNAFreqDict

# Count RNA nucleotides
# countRNANucFrequency()
def countRNANucFrequency(rna_seq):
    tmpRNAFreqDict = {"A":0, "C":0, "G":0, "U":0}
    for nuc in rna_seq:
        tmpRNAFreqDict[nuc] += 1
    return tmpRNAFreqDict

def DNA_ReserveComplement(dna_seq):
    return ''.join([DNA_Reverse[nuc] for nuc in dna_seq])[::-1]

def DNAtoRNA_Transcription(dna_seq):
    return dna_seq.replace("T", "U")

dna ="ATTAATTTTGGG"
DNA_ReserveComplement(dna)
DNAtoRNA_Transcription(dna)