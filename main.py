""" File containing the main() function to run our program. Since it tents to take
too long to get and analize the sequences it has been implemented using multiprocessing
    @authors: Lucas, Josh, Amy, Daniele
"""
import motifFinder as mf
import time
from multiprocessing import Pool
from functools import partial
import sys


# Running the program:
print('Running...')
print('1/3 - Applying Gibbs Sampling...')
# List for the multiprogramming pool
bestMotifsRes = []
# List for the dictionaries returned by gibbsSampling()
bestMotifsDict = []
# List for the motifs in the dictionary returned by gibbsSampling()
bestMotifs = []
fileToRead = mf.readInput('TraR.txt')
# Only used to print values on the screen for reference
startTime = time.time()

# Number of times we'll run gibbsSampling()
iterable = range(2000)
# Max and min motifs sizes
minMotifSize = 12
maxMotifSize = 20
for kmerSize in range(minMotifSize, maxMotifSize + 1):
    pool = Pool()
    function = partial(mf.gibbsSampling, fileToRead, kmerSize, 200)
    bestMotifsRes = pool.map_async(function, iterable)
    # Updates the percentage on the console
    while not bestMotifsRes.ready():
        remaining = 100 - (bestMotifsRes._number_left * bestMotifsRes._chunksize / (len(iterable) / 100))
        sys.stderr.write('\r\033[2KK-mer: %d of %d, Progress: %d%%' % (kmerSize, maxMotifSize, remaining))
        sys.stderr.flush()
        time.sleep(.1)
    pool.close()
    pool.join()
    for dicts in bestMotifsRes.get():
        bestMotifsDict.append(dicts)
# print(bestMotifsDict)
print('Gibbs Sampling Done!')

# Copy the motifs that we got from gibbsSampling to a list so we can analize them later
print('2/3 - Creating the Best Motifs List...')
for dictionary in bestMotifsDict:
    for motifs in dictionary['motifs']:
        bestMotifs.append(motifs)
# print(bestMotifs)

# Non multi-threaded part. Only used for comparison purposes
"""mot = []
percentage = 0

for i in range(0, 2000):
    if int(100 * (i / 2000)) != percentage:
        print(percentage, '%')
    for kmerSize in range(12, 21):
        curBestMotif = mf.gibbsSampling(mf.readInput('TraR.txt'), kmerSize, 200, 0)
        for motif in curBestMotif:
            bestMotifs.append(motif)
    percentage = int(100 * (i / 2000))

print('Gibbs Sampling Done!')
print('2/2 - Applying Profile to Genome...')
print('Applying Profile to Genome Done!')"""

print('Done!')

# Calulates and converts total running time to hh:mm:ss
finalTime = time.time() - startTime
minutes, seconds = divmod(finalTime, 60)
hours, minutes = divmod(minutes, 60)
print('Total running time: %d:%d:%d' % (hours, minutes, seconds))
