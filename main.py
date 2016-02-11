"""
@authors: Lucas, Josh, Amy, Daniele
"""
import motifFinder as mf
import time


# Running the program:
print('Running...')
print('Applying Gibbs Sampling...')
# Only used to print values on the screen for reference
startTime = time.time()
bestMotifs = []
# Only used to print values on the screen for reference
percentage = 0

for i in range(0, 2000):
    if int(100 * (i / 2000)) != percentage:
        print(percentage, '%')
    for kmerSize in range(4, 30):
        curBestMotif = mf.gibbsSampling(mf.readInput('TraR.txt'), kmerSize, 200)
        for motif in curBestMotif:
            bestMotifs.append(motif)
    percentage = int(100 * (i / 2000))

print('Gibbs Sampling Done!')
print('Applying Profile to Genome...')
profile = mf.constructProfile(bestMotifs)
singleScores = []
for motif in bestMotifs:
    singleScores.append(mf.getSingleScore(profile, motif))
print('Applying Profile to Genome Done!')
print('Done!')
print('Total running time: ', time.time() - startTime, 'seconds')
