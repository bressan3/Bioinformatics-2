""" Library containing the extra functions we need for the main part of the assignment
    @authors: Lucas, Josh, Amy, Daniele
"""
import motifFinder as mf
import json
from time import localtime, strftime


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


def gatherFinalResults(profile, sequences, worstScoringMotif, strand):
    """ Find all locations in the Agrobacterium tumefaciens C58 genome where our
        profile achieves a score as high as the worst scoring motif (returned by gibbsSampling())
    Args:
        profile (dictionary): Profile of the motifs returned by gibbsSampling()
        sequences (string): List of DNA sequences
        worstScoringMotif (dictionary): Dictionary containing important info about the worst scoring motif
            (from the best scoring ones got from gibbsSampling())
        strand (int): 1 for the sequence from the file / 2 from the reverse complement
    Returns:
        dictionary: A list of dictionaries of the motifs that score as high as worstScoringMotif
    """
    dnaScores = []
    for sequenceNumber in range(0, len(sequences)):
        # Applies a profile (gets single score for each subsequence) for each kmer in the current sequence
        applyProfileScores = mf.applyProfile(profile, sequences[sequenceNumber])
        for i in range(0, len(applyProfileScores)):
            # Find motifs that has a score as high as the worst scoring motif's score
            if applyProfileScores[i] >= worstScoringMotif['Score']:
                # Searches for the closest protein start codon
                str1 = sequences[sequenceNumber][0: i + 1]
                str2 = sequences[sequenceNumber][i + len(profile): len(sequences[sequenceNumber])]
                str1Start = str1.rfind('ATG')
                str2Start = str2.find('ATG')
                if str1Start == -1:
                    str1Start = 1000
                if str2Start == -1:
                    str2Start = 1000
                if str1Start < str2Start:
                    start = str1Start
                elif str1Start > str2Start:
                    start = str2Start + i + len(profile)
                else:
                    start = 'None'
                # Adds the found motif and its info to a list of dictionaries that is gonna be returned by the function
                dnaScores.append({'50mer-Sequence': sequences[sequenceNumber][i: i + len(profile) + (50 - len(profile))], 'Position': i, 'Score': applyProfileScores[i], 'Closest-Protein-Coding-Gene': start, 'DNASequence#': sequenceNumber + 1, 'Strand #': strand})

    return dnaScores


def writeFile(fileName, profile, singleScores, bestMotifs, dnaScores, bestMotif):
    """ Writes the gathered data to the Results.json file
    Args:
        fileName (string): Location of the file we want to write on
        profile (dictionary): Profile of the motifs returned by gibbsSampling()
        singleScores (float): List of the single scores of each motif returned by gibbsSampling()
        bestMotifs (string): List containing the best scoring motifs from gibbsSampling()
        dnaScores (dictionary): List of fictionaries containing the info of the motifs returned
            by gatherFinalResults()
        bestMotif (dictionary): Dictionary returned by gibbsSampling() where we have the best scoring motif
    Returns:
        None: No return
    """
    with open(fileName, 'w+') as f:
        f.write(strftime("Created on: %Y-%m-%d %H:%M:%S\n", localtime()))
        f.write('Best Motifs: ')
        f.write('\n')
        json.dump(bestMotif, f)
        f.write('\n')
        f.write('Motifs Profile: ')
        f.write('\n')
        json.dump(profile, f)
        f.write('\n')
        f.write('Single Scores: ')
        f.write('\n')
        for i in range(0, len(singleScores)):
            json.dump(bestMotifs[i], f)
            f.write(': ')
            json.dump(singleScores[i], f)
            f.write('\n')
        f.write('Motifs that have a better score than the worst scoring one: ')
        f.write('\n')
        for scores in dnaScores:
            json.dump(scores, f)
            f.write('\n')
