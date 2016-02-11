""" File containing the functions test for motifFinder.py """
import motifFinder


# Tests ------------------------------------------------------------------

print(motifFinder.readInput('TraR.txt'))
print(motifFinder.randomStart(['ACACGTAC', 'CCACGTCACA', 'TTCGTCGTACG'], 4))
print(motifFinder.getMotif(['ACACGTAC', 'CCACGTCACA', 'TTCGTCGTACG'], [3, 5, 2], 4))
print(motifFinder.constructProfile(['CGTA', 'TCAC', 'CGTC']))
print(motifFinder.getSingleScore(motifFinder.constructProfile(['CGTA', 'TCAC', 'CGTC']), 'CACA'))
print(motifFinder.applyProfile(motifFinder.constructProfile(['CGTA', 'TCAC', 'CGTC']), 'CCACGTCACA'))
print(motifFinder.randomlySelect([0.014994, 0.001249, 0.000833, 0.033736, 0.000833, 0.009996, 0.002499]))
print(motifFinder.nucleotideFrequencies(['ACACGTAC', 'CCACGTCACA', 'TTCGTCGTACG']))
print(motifFinder.scoreProfile([{'A': 0.142857, 'C': 0.428571, 'G': 0.142857, 'T': 0.285714}, {'A': 0.142857, 'C': 0.285714, 'G': 0.428571, 'T': 0.142857}, {'A': 0.285714, 'C': 0.142857, 'G': 0.142857, 'T': 0.428571}, {'A': 0.285714, 'C': 0.428571, 'G': 0.142857, 'T': 0.142857}], {'A': 0.241379, 'C': 0.379310, 'G': 0.172413, 'T': 0.206897}))
