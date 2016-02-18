""" Library containing the extra functions we need for the assignment
    @authors: Lucas, Josh, Amy, Daniele
"""


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
