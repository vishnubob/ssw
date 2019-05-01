import six
from six.moves import range

__all__ = (
    "NucleotideTable",
    "NucleotideAlphabet",
    "NucleotideComplementTable",
    "nucleotide_reverse_complement",
)

def _build_compliment_table():
    _ctable = list(map(chr, range(0xff + 1)))
    for symbol in NucleotideAlphabet:
        complement = NucleotideTable[symbol]["complement"]
        _ctable[ord(symbol)] = complement
        _ctable[ord(symbol.lower())] = complement.lower()
    return str.join('', _ctable)

def _iupac_info(complement, matches):
    return {
        "complement": complement.upper(),
        "matches": tuple(matches.upper())
    }

NucleotideTable = {
    'A': _iupac_info('T', 'A'),
    'G': _iupac_info('C', 'G'),
    'C': _iupac_info('G', 'C'),
    'T': _iupac_info('A', 'T'),
    'U': _iupac_info('A', 'T'),
    'M': _iupac_info('K', 'AC'),
    'R': _iupac_info('Y', 'AG'),
    'Y': _iupac_info('R', 'CT'),
    'S': _iupac_info('S', 'CG'),
    'W': _iupac_info('W', 'AT'),
    'K': _iupac_info('M', 'GT'),
    'B': _iupac_info('V', 'CGT'),
    'D': _iupac_info('H', 'AGT'),
    'H': _iupac_info('D', 'ACT'),
    'V': _iupac_info('B', 'ACG'),
    'N': _iupac_info('N', 'AGCT')
}

NucleotideAlphabet = str.join('', tuple(NucleotideTable.keys()))
NucleotideComplimentTable = _build_compliment_table()

def nucleotide_reverse_complement(sequence):
    return sequence[::-1].translate(NucleotideComplimentTable)
