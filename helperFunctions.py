""" Library containing the extra functions we need for the assignment
    @authors: Lucas, Josh, Amy, Daniele
"""
import motifFinder as mf


def getSecondStrand(sequences):
    """ Returns the complement DNA strand and converts it to 3' - 5' given a list of sequences
    Args:
        sequences (string): List of sequences
    Returns:
        string: List of the complement sequences
    """
    compDNA = []
    for dna in sequences:
        compDNAAux = dna.replace('A', 't')
        compDNAAux = compDNAAux.replace('T', 'a')
        compDNAAux = compDNAAux.replace('C', 'g')
        compDNAAux = compDNAAux.replace('G', 'c')
        compDNA.append(compDNAAux.upper())

    for i in range(0, len(compDNA)):
        compDNA[i] = compDNA[i][::-1]

    return compDNA


def gatherScores(profile, sequences, worstScoringMotif, strand):
    """ Finds all locations in genome (sequences) where the given
    profile achieves a score as high as the worst scoring motif
    Args:
        profile (dictionary): Profile of the motifs returned by gibbsSampling
        sequences (string): List of all DNA sequences from TraR.txt
        worstScoringMotif (dictionary): Contains important info about the worst scoring motif (from gibbsSampling)
        such as Motif and its Score
        strand (int): Current strand's number (1 for the original DNA strand and 2 for the complement)
    Returns:
        dictionary: A list of dictionaries containing the locations where the given
        profile achieves a score as high as the worst scoring motif
    """
    for sequenceNumber in range(0, len(sequences)):
        dnaScores = []
        # Applies a profile (gets single score for each subsequence) for each kmer in the current sequence
        applyProfileScores = mf.applyProfile(profile, sequences[sequenceNumber])
        for i in range(0, len(applyProfileScores)):
            if applyProfileScores[i] >= worstScoringMotif['Score']:
                dnaScores.append({'Sequence': sequences[sequenceNumber][i: i + len(profile)], 'Position': i, 'Score': applyProfileScores[i], 'DNASequence#': sequenceNumber + 1, 'Strand #': strand})

        return dnaScores
