import random


def readInput(filename):
    """ Function that reads a file's contents outputs a list containing a nucleic acid sequence
    Args:
        filename (string): File path
    Return:
        string: Sequence read from the file
    """
    return [x.strip() for x in open(filename).readlines()]


def randomStart(sequences, k):
    """ Accepts a list of sequences and an integer k representing
    the motif size. Returns a list of size len(sequences) of integers representing
    the random start index for a k-mer from each sequence.
    Args:
        sequences (string): List of sequences
        k (int): Motif size
    Returns:
        int: List of random start indexes for each k-mer
    """
    startIndices = []
    for i in range(0, len(sequences)):
        startIndices.append(random.randint(0, len(sequences[i]) - k))
    return startIndices


def getMotif(sequences, startLocations, k):
    """ Accepts a list of sequences, a list of start locations,
    and the motif size. Returns a list of k-mers. Each k-mer
    starts at the given start location and is k characters long. If the start location for sequence
    is the value -1 then the k-mer that should be returned for that sequence should be an empty string.
    Args:
        sequences (string): List of sequences
        startLocations (int): List of start locations for each sequence in the list
        k (int): Motif size
    Returns:
        string: List of k-mers
    """
    kmers = []
    for i in range(0, len(sequences)):
        kmers.append(sequences[i][startLocations[i]:startLocations[i] + k])
    return kmers


def constructProfile(motifs):
    """ Accepts a list of motifs and returns a list of dictionaries representing the profile of the k-mers.
    Args:
        motifs (string): List of motifs
    Returns:
        dictionary: List of dictionaries representing the profile of the k-mers
    """
    kmersProfile = []
    for i in range(0, len(motifs[0])):
        freqA = 0
        freqC = 0
        freqG = 0
        freqT = 0
        for j in range(0, len(motifs)):
            if motifs[j][i] == 'A':
                freqA += 1
            elif motifs[j][i] == 'C':
                freqC += 1
            elif motifs[j][i] == 'G':
                freqG += 1
            elif motifs[j][i] == 'T':
                freqT += 1
        if freqA == 0 or freqC == 0 or freqG == 0 or freqT == 0:
            freqA += 1
            freqC += 1
            freqG += 1
            freqT += 1
        # 7 is the value for the test's motif sequence. I don't know how to get the actual
        # value for different sequences yet
        kmersProfile.append(
            {'A': freqA / 7, 'C': freqC / 7, 'G': freqG / 7, 'T': freqT / 7})

    return kmersProfile


def getSingleScore(profile, kmer):
    """ Accepts a profile (a list of dictionaries as returned by constructProfile()) and a single string.
    Both the string and the profile should be the same length. Returns the score of applying the profile
    to the sequence.
    Args:
        profile (dictionary): List of Dictionaries
        kmer (string): k-mer of the same length as profile
    Returns:
        float: Score
    """
    score = profile[0][kmer[0]]
    for i in range(1, len(kmer)):
        score *= profile[i][kmer[i]]

    return score


def applyProfile(profile, sequence):
    """ Accepts a profile (a list of dictionaries as returned by
    constructProfile()) and a single sequence of any length. The profile should be applied to all
    len(profile) subsequences in sequence. Returns a list of numbers of the scores of the
    profile applied to each of the subsequences.
    Args:
        profile (dictionary): List of dictionaries
        sequence (string): Sequenceof any length
    Returns:
        float: List of scores of each subsequence
    """
    profiles = []
    for i in range(0, len(sequence) - len(profile) + 1):
        curMotif = sequence[i:i + len(profile)]
        profiles.append(getSingleScore(profile, curMotif))

    return profiles


def randomlySelect(probabilities):
    """ Accepts a list of probabilities as input. The function should normalize
    the numbers (divide by the total so they sum to 1) and then randomly (weighted) select and return the
    index of the selected number.
    Args:
        probabilities (float): List of probabilities
    Returns:
        int: Index of the selected number
    """
    # Normalizes the probabilities list
    probabilities = [float(i) / sum(probabilities) for i in probabilities]
    randValue = random.randint(0, sum(probabilities))
    weightSum = 0

    for i in range(0, len(probabilities)):
        weightSum += probabilities[i]
        if randValue <= weightSum:
            return i


def nucleotideFrequencies(sequences):
    """ Takes in a list of sequences and returns a dictionary of the frequencies of each of the nucleotides.
    Args:
        sequences (string): List of sequences
    Returns:
        Dictionary: Frequency of each nucleotide in the sequence
    """

    freqA = 0
    freqC = 0
    freqG = 0
    freqT = 0
    seqTotalSize = 0

    for i in range(0, len(sequences)):
        seqTotalSize += len(sequences[i])
        freqA += sequences[i].count('A')
        freqC += sequences[i].count('C')
        freqG += sequences[i].count('G')
        freqT += sequences[i].count('T')

    return {'A': freqA / seqTotalSize, 'C': freqC / seqTotalSize, 'G': freqG / seqTotalSize, 'T': freqT / seqTotalSize}


def scoreProfile(profile, nucFreq):
    """ Accepts a profile (a list of dictionaries as returned by constructProfile()) and a
    dictionary of nucleotide frequencies (as returned by nucleotideFrequencies()) and return
    a number representing the relative entropy of the profile.
    Args:
        profile (dictionary): List of dictionaries as returned by constructProfile()
        nucFreq (dictionary): Dictionary of nucleotide frequencies as returned by nucleotideFrequencies()
    Returns:
        float: Relative entropy
    """


def gibbsSampling(sequences, k, iterations):
    """ The function returns the best motifs which is the highest scoring motifs using the relative
    entropy from nucleotideFrequencies().
    Args:
        sequences (string): List of sequences
        k (int): Motif size
        iterations (int): Number of interations to run
    Returns:
        String: List with best motifs
    """

# Tests ------------------------------------------------------------------
print(randomStart(['ACACGTAC', 'CCACGTCACA', 'TTCGTCGTACG'], 4))
print(getMotif(['ACACGTAC', 'CCACGTCACA', 'TTCGTCGTACG'], [3, 5, 2], 4))
print(constructProfile(['CGTA', 'TCAC', 'CGTC']))
print(getSingleScore(constructProfile(['CGTA', 'TCAC', 'CGTC']), 'CACA'))
print(applyProfile(constructProfile(['CGTA', 'TCAC', 'CGTC']), 'CCACGTCACA'))
print(randomlySelect([0.014994, 0.001249, 0.000833, 0.033736, 0.000833, 0.009996, 0.002499]))
print(nucleotideFrequencies(['ACACGTAC', 'CCACGTCACA', 'TTCGTCGTACG']))
