DNANucleotides = ["A", "C", "G", "U"]
RNANucleotides = ["A", "C", "G", "T"]

# Validate DNA sequence
# valiDNASeq()
def valiDNASeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in DNANucleotides:
            return False
    return tmpseq

# Validate RNA sequence
# valiRNASeq()
def valiRNASeq(rna_seq):
    tmpseq = rna_seq.upper()
    for nuc in tmpseq:
        if nuc not in RNANucleotides:
            Print("Errors in RNA sequence")
    return tmpseq
