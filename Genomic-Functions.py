
# September 15, 2023
# by Olivia Radcliiffe


# 1.takes a DNA sequence S as input and returns the
# percentage of the DNA which is “G” or “C”.
def GCPercentage(S):
    if len(S) > 0:
        percent = (S.count('G') + S.count('C'))/len(S) * 100
    else:
        percent = 0

    return percent

# 2. takes a DNA sequence S as input and returns 
# its reverse compliment.
def ReverseCompliment(S):
    
    base_pairs_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    revcompliment = ""
    
    for base in S:
        revcompliment += base_pairs_dict.get(base, 'X')

    # checks for unknown bases
    if 'X' in revcompliment or not S:
        return "Incorrect input sequence"
    
    # reverse sequence
    return revcompliment[::-1]

# 3. takes a DNA sequence as input and returns both
#  transcribed mRNA, and translated protein sequence.
def RNAtoProtein(DNA):

    base_pairs_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    compliment = ""
    
    for base in DNA:
        compliment += base_pairs_dict.get(base, 'X')

    # transcribe to mRNA
    mRNA = compliment.replace('T', 'U')

    codon_dict = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
    
    # find start codon
    start = mRNA.index("AUG")
    proSeq = ""

    # translate mRNA to protein
    for i in range(start,len(mRNA), 3):
        code = codon_dict.get(mRNA[i:i+3], 'X')

        # check for stopping
        if code == "*":
            return
        proSeq += code

        # check for inccorect bases
        if 'X' in proSeq or not DNA:
            return "Incorrect input sequence"
        
    return mRNA, proSeq

#4. receives two DNA sequences S1 and S2 with equal
# lengths as input and returns the Hamming distance
# (the number of mismatches) between them.
def HammingDistance(S1, S2):

    # checks the lengths
    if len(S1) != len(S2):
        return "Sequences are not equal length!"
    
    # counts mismatches
    mismatches = 0
    if S1 != S2:
        for i in range(len(S1)):
            if S1[i] != S2[i]:
                mismatches += 1
    
    return mismatches

# 5. takes a DNA sequence and k-mer Pattern and 
# returns the number of times that the k-mer Pattern
# appears as a substring of DNA
def Count(DNA, Pattern):
    k = len(Pattern)
    freq = 0

    # counts apperances of Pattern 
    for i in range(len(DNA)-k+1):
        if DNA[i:i+k] == Pattern:
            freq += 1
    
    return freq


# 6. takes a DNA string and an integer k, and 
# returns all the most frequent k-mers in the 
# sequence, and the indices that the k-mer is 
# appears in the sequence
def mostFrequentKmer(DNA, k):
    max = 0
    mostFreq = []

    # finds all kmers and counts occurences
    for start in range(len(DNA)-k+1):
        kmer = DNA[start:start+k]
        num = Count(DNA, kmer)
        
        if kmer not in mostFreq:
            if num > max:
                max = num
                mostFreq = [(kmer)]
            elif num == max:
                mostFreq.append(kmer)

    # finds all indicies
    index = []
    for seq in mostFreq:
        for i in range(len(DNA)):
            if DNA.startswith(seq, i):
                index.append(i)

    return mostFreq, index

# 7. takes a DNA string Pattern and an integer d, 
# and returns the list of all the neighbors of the
# Pattern.
def dNeighborhood(Pattern, d):

    # special cases
    if d == 0:
        return [Pattern]

    if len(Pattern) == 1:
        return ["A", "C", "G", "T"]

    neighbors = []

    # generate neighbors
    def recursive_neighbors(cur_pattern, hamming):
        if hamming == 0:
            if cur_pattern not in neighbors:
                neighbors.append(cur_pattern)
            return

        for base in "ACGT":
            for i in range(len(cur_pattern)):
                new_pattern = cur_pattern[:i] + base + cur_pattern[i + 1:]
                recursive_neighbors(new_pattern, hamming-1)

    recursive_neighbors(Pattern, d)

    return neighbors


# 8. takes DNA sequence, Pattern, and maximum 
# hamming distance d, and returns approximate 
# occurrences of Pattern and their starting 
# positions in the DNA sequence
def approxOccurence(DNA, Pattern, max_d):
    
    approxOccur = set()

    # finds all kmers and calcs the hamming distance
    for i in range(len(DNA)-len(Pattern)+1):
        kmer = DNA[i:i+len(Pattern)]
        hamming = HammingDistance(Pattern, kmer)
        if hamming <= max_d:
            approxOccur.add((kmer,i))

    return list(approxOccur)

# 9. takes a Genome sequence, and integers k, L, 
# and t, and returns all distinct k-mers forming 
# (L, t)-clumps in Genome.
def clumpsInString(k, L, t, DNA):
    clumps = []
    kmers = []

    #find kmers
    for start in range(len(DNA)-k+1):
        kmers.append(DNA[start:start+k])

    #check if kmer appears in L intervals >= t times
    for i in range(len(DNA)-L):
        for kmer in kmers:
            if kmer not in clumps:
                if Count(DNA[i:i+L], kmer) >= t:
                    clumps.append(kmer)

    return clumps


# 10. takes an integer k and two sequences and 
# returns all k-mers shared by these strings, and 
# the ordered pairs (x, y) corresponding to starting
# positions of these k-mers in the respective 
# strings.
def sharedKmers(k, S1, S2):

    shared = []
    starts = []


    for i in range(len(S1)-k+1):
        kmer = S1[i:i+k]

        # check if kmer is shared
        if kmer in S2:
            shared.append(kmer)
            starts.append((i,S2.index(kmer)))

        # find kmer compliment
        base_pairs_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        compkmer = ""
        for base in kmer:
            compkmer += base_pairs_dict.get(base, 'X')

        # check if compliment is in s2
        if compkmer in S2:
            shared.append((kmer, compkmer))
            starts.append((i,S2.index(compkmer)))

    return shared, starts

# 11. takes a  a sequence and return the skew array, 
# plot the skew diagram, and the position where the
# skew is minimum.
import matplotlib.pyplot as plt

def Skew(Genome):

    skewArray = [0]
    indicies = [0]
    for i in range(len(Genome)):
        # determine next skew
        if Genome[i] == "A" or Genome[i] == "T":
            skewArray.append(skewArray[i])
        elif Genome[i] == "C":
            skewArray.append(skewArray[i]-1)
        elif Genome[i] == "G":
            skewArray.append(skewArray[i]+1)

        # add index
        indicies.append(i+1)  

    # plot skew diagram
    plt.title("Skew Diagram")
    plt.xlabel("Position")
    plt.ylabel("Skew")
    plt.plot(indicies, skewArray)
    plt.plot(skewArray.index(min(skewArray)), min(skewArray), marker='*')
    plt.show()
    
    return skewArray

# 12. takes a sequence S and an integer k, and 
# returns the frequency array of k-mers, along with
#  the k-mers and indices array in lexicographic 
# order.
import numpy as np
import itertools 

def frequenceKmer(S, k):
    kmers = [''.join(p) for p in itertools.product("ACGT", repeat=k)]
    
    indicies = np.arange(0,len(kmers))

    # calc frequency
    freqArray = []
    for kmer in kmers:
        freqArray.append(Count(S, kmer))

    return kmers, freqArray, indicies


# Tests all functions - not sure if it is necessary
def testFunctions():
    print("GCPercentage: ", GCPercentage("AGCTTCAGTTCA"))
    print("GCPercentage: ", GCPercentage(""))

    print("\n")

    print("ReverseCOmpliment: ", ReverseCompliment("AGCTTCAGTTCA"))
    print("ReverseCOmpliment: ", ReverseCompliment("AGCTTCPAGTTCA"))
    print("ReverseCOmpliment: ", ReverseCompliment(""))

    print("\n")

    print("mRNA, proSeq: ", RNAtoProtein("TACAGCTTCAGTTCA"))

    print("\n")

    print("HammingDistance: ", HammingDistance("ATTGACGGAT", "AGTGTCGAAT"))

    print("\n")

    print("Counting K-mers: ", Count("ACAACTATGCATACTATCGGGAACTATCCT","ACTAT"))

    print("\n")

    print("most frequent K-mer: ", mostFrequentKmer("ACAACTATGCATCACTATCGGGAACTATCCT", 5))

    print("\n")

    print("1-neighborhood of CCAGTCAATG: ", dNeighborhood("CCAGTCAATG", 1))
    print("Number of neighbors: ", len(dNeighborhood("CCAGTCAATG", 1)))

    print("\n")

    print("2-neighborhood of CCAGTCAATG: ", dNeighborhood("CCAGTCAATG", 2))
    print("Number of neighbors: ", len(dNeighborhood("CCAGTCAATG", 2)))

    print("\n")

    print("Approx occurences: of CCAGTCAATG", approxOccurence("CCAGTCAATG", "CCT", 2))

    print("\n")

    print("clumps: ", clumpsInString(5, 75, 4, "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC"))
    print("clumps: ", clumpsInString(4, 25, 3, "gatcagcatatgcagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac"))
    
    print("\n")

    shared, starts = sharedKmers(3, "AAACTCATC", "TTTCAAATC")
    print("shared: ", shared, "starts: ", starts)

    print("\n")

    Skew("CATGGGCATCGGCCATACGCC")

    print("\n")
    kmers, freqArray, indicies = frequenceKmer("AAGCAAAGGTGGG", 2)
    print("kmer: ", kmers)
    print("index: ", indicies)
    print("frequency: ", freqArray)

def main():
    testFunctions()

if __name__ == "__main__":
    main()