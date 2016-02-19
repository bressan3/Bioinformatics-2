""" File containing the main() function to run our program. Since it tents to take
too long to get and analize the sequences it has been implemented using multiprocessing
    @authors: Lucas, Josh, Amy, Daniele
"""
import motifFinder as mf
import helperFunctions as hf
from time import time, sleep
from multiprocessing import Pool
from functools import partial
import sys


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

# Getting the reverse complement
reverseComplement = hf.getSecondStrand(fileToRead)
dnaScores = hf.gatherFinalResults(profile, fileToRead, worstScoringMotif, 1) + hf.gatherFinalResults(profile, reverseComplement, worstScoringMotif, 2)

print('Applying Profile to Genome Done!')

# Reporting the results to a file
hf.writeFile('Results/Profile.json', profile, singleScores, bestMotifs, dnaScores, bestMotifsDict[maxScoreIndex])

print('Done!')


# --------------------------------------------------------------------------------------------------------------


# Calulates and converts total running time to hh:mm:ss
finalTime = time() - startTime
minutes, seconds = divmod(finalTime, 60)
hours, minutes = divmod(minutes, 60)
print('Total running time: %d:%d:%d' % (hours, minutes, seconds))
