""" File containing the main() function to run our program. Since it tents to take
too long to get and analize the sequences it has been implemented using multiprocessing
    @authors: Lucas, Josh, Amy, Daniele
"""
import motifFinder as mf
import helperFunctions as hf
from time import time, sleep, localtime, strftime
from multiprocessing import Pool
from functools import partial
import sys
import json


# Running the program:
print('Running...')
print('1/3 - Applying Gibbs Sampling...')
# List for the multiprogramming pool
bestMotifsRes = []
# List of best motifs returned by gibbsSampling()
bestMotifsDict = []
fileToRead = mf.readInput('TraR.txt')
# Only used to print time values on the screen for reference
startTime = time()

# Number of times we'll run gibbsSampling()
iterable = range(2000)
# Max and min motifs sizes
minMotifSize = 12
maxMotifSize = 20
for kmerSize in range(minMotifSize, maxMotifSize + 1):
    pool = Pool()
    function = partial(mf.gibbsSampling, fileToRead, kmerSize, 200)
    bestMotifsRes = pool.map_async(function, iterable)
    # Updates the percentage on the console screen
    while not bestMotifsRes.ready():
        remaining = 100 - (bestMotifsRes._number_left * bestMotifsRes._chunksize / (len(iterable) / 100))
        sys.stderr.write('\r\033[2KK-mer: %d of %d, Progress: %d%%' % (kmerSize, maxMotifSize, remaining))
        sys.stderr.flush()
        sleep(.1)
    pool.close()
    pool.join()
    for dicts in bestMotifsRes.get():
        bestMotifsDict.append(dicts)
print('\n')
# print(bestMotifsDict)
print('Gibbs Sampling Done!')


# --------------------------------------------------------------------------------------------------------------


# Copy the best motifs that we got from gibbsSampling to a list so we can analize them later
print('2/3 - Gathering the Best Motifs...')
# Finds the index in the list of dictionaries returned by gibbsSampling() where maxScore is
maxScoreIndex = next(index for (index, d) in enumerate(bestMotifsDict) if d['highestScore'] == max(item['highestScore'] for item in bestMotifsDict))
# Adds the best found motifs into a list
bestMotifs = bestMotifsDict[maxScoreIndex]['motifs']
# print('Max score ========', maxScore, 'Index ===========', maxScoreIndex, 'Motifs ========', bestMotifs)
# Creates a file that reports the best scoring motifs, k and the scoreProfile()
print('Gathering the Best Motifs Done!')


# --------------------------------------------------------------------------------------------------------------


print('3/3 - Applying Profile to Genome...')
profile = mf.constructProfile(bestMotifs)
singleScores = []
# Get each motif's single score and puts the worstScoring motif into a dictionary
for motif in bestMotifs:
    singleScores.append(mf.getSingleScore(profile, motif))
worstScoringMotif = {'Motif': bestMotifs[singleScores.index(min(singleScores))], 'Score': min(singleScores)}
# List of dictionaries containing the positions where the profile achieves a score as high as the
# worst scoring motif is located in each given DNA sequence
dnaScores = []
for sequenceNumber in range(0, len(fileToRead)):
    # Applies a profile (gets single score for each subsequence) for each kmer in the current sequence
    applyProfileScores = mf.applyProfile(profile, fileToRead[sequenceNumber])
    for i in range(0, len(applyProfileScores)):
        if applyProfileScores[i] >= worstScoringMotif['Score']:
            str1 = fileToRead[sequenceNumber][0: i + 1]
            str2 = fileToRead[sequenceNumber][i + len(profile): len(fileToRead[sequenceNumber])]
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
            dnaScores.append({'50mer-Sequence': fileToRead[sequenceNumber][i: i + len(profile) + (50 - len(profile))], 'Position': i, 'Score': applyProfileScores[i], 'Closest-Protein-Coding-Gene': start, 'DNASequence#': sequenceNumber + 1, 'Strand #': 1})

# Getting the reverse complement
reverseComplement = hf.getSecondStrand(fileToRead)
for sequenceNumber in range(0, len(reverseComplement)):
    # Applies a profile (gets single score for each subsequence) for each kmer in the current sequence
    applyProfileScores = mf.applyProfile(profile, reverseComplement[sequenceNumber])
    for i in range(0, len(applyProfileScores)):
        if applyProfileScores[i] >= worstScoringMotif['Score']:
            str1 = fileToRead[sequenceNumber][0: i + 1]
            str2 = fileToRead[sequenceNumber][i + len(profile): len(fileToRead[sequenceNumber])]
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
            dnaScores.append({'50mer-Sequence': reverseComplement[sequenceNumber][i: i + len(profile) + (50 - len(profile))], 'Position': i, 'Score': applyProfileScores[i], 'Closest-Protein-Coding-Gene': start, 'DNASequence#': sequenceNumber + 1, 'Strand #': 2})

print('Applying Profile to Genome Done!')

# Reporting the results
with open('Results/Profile.json', 'w+') as f:
    f.write(strftime("Created on: %Y-%m-%d %H:%M:%S\n", localtime()))
    f.write('Best Motifs: ')
    f.write('\n')
    json.dump(bestMotifsDict[maxScoreIndex], f)
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
print('Done!')


# --------------------------------------------------------------------------------------------------------------


# Calulates and converts total running time to hh:mm:ss
finalTime = time() - startTime
minutes, seconds = divmod(finalTime, 60)
hours, minutes = divmod(minutes, 60)
print('Total running time: %d:%d:%d' % (hours, minutes, seconds))
