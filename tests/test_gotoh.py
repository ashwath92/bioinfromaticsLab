import sys,os
import pytest
sys.path.append(os.path.abspath('.'))
from prakt.gt import GotohBase
from gotoh import Gotoh

def test_instance():
    """Check inheritance."""
    assert issubclass(Gotoh, GotohBase)
    assert isinstance(Gotoh(), GotohBase)


def test_example_success():
    """This calls the run method of the Needleman-Wunsch program
    and tests if it works as expected (positive test)"""
    gt = Gotoh()
    seq_fasta_1 = os.path.join('data','sequences','seq1.fasta')
    seq_fasta_2 = os.path.join('data','sequences','seq2.fasta')
    seq_fasta_3 = os.path.join('data','sequences','seq3.fasta')
    result = gt.run(seq_fasta_1,
                    seq_fasta_2,
                    'pam250',
                    -11,
                    -1,
                    True)
    (id_seq1, seq1, id_seq2, seq2, score, alignments,num_alignments) = result

    assert id_seq1 == "ID1"
    assert id_seq2 == "ID2"
    assert seq1 == "ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN"
    assert seq2 == "RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT"
    assert score == 33
    assert alignments== [['ILDMDVVEGSAARFDCKVEG-YPDPEVMWFKDDNP---VKESRHFQIDYDEEGN', 
    'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFV----SQTT', ':::::::**::::::*:::: **:::::*:::::*   ::::***:    ::::'], ['ILDMDVVEGSAARFDCKVEG-YPDPEVMWFKDDNP---VKESRHFQIDYDEEGN', 'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHF----VSQTT', ':::::::**::::::*:::: **:::::*:::::*   ::::***    :::::']]
    assert num_alignments == 2


def test_example_invalid_format_fail():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for the failure is invalid file format: the first line does not start with >"""

    gt = Gotoh()
    seq_fasta_1 = os.path.join('data','sequences','seq1.fasta')
    seq_fasta_2 = os.path.join('data','sequences','Invalid_format.fasta')
    with pytest.raises(SystemExit) as InvalidFileException:
        result = gt.run(seq_fasta_2,
                    seq_fasta_1,
                    'pam250',
                    -11,
                    -1,
                    False)
        (id_seq1, seq1, id_seq2, seq2, score, alignments, num_alignments) = result

        assert InvalidFileException.type == SystemExit
        assert InvalidFileException.code == 1

def test_example_invalid_characters_fail():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for failure is non-amino acid characters in file 2 (error code 12)"""

    gt = Gotoh()
    seq_fasta_1 = os.path.join('data','sequences','seq1.fasta')
    seq_fasta_2 = os.path.join('data','sequences','Invalid_characters.fasta')
    with pytest.raises(SystemExit) as InvalidCharactersException:
        result = gt.run(seq_fasta_1,
                    seq_fasta_2,
                    'pam250',
                    -11,
                    -1,
                    False)
        (id_seq1, seq1, id_seq2, seq2, score, alignments, num_alignments) = result

        assert InvalidCharactersException.type == SystemExit
        assert InvalidCharactersException.code == 12


def test_too_few_arguments():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for failure is non-amino acid characters in file 2 (error code 12)"""

    gt = Gotoh()
    seq_fasta_1 = os.path.join('data','sequences','seq1.fasta')
    seq_fasta_2 = os.path.join('data','sequences','seq2.fasta')
    # test is a variable which becomes True when there are too few arguments
    test = False
    try:
        with pytest.raises(SystemExit) as TooFewArguments:
            result = gt.run(seq_fasta_1,
                        seq_fasta_2,
                        'pam250',
                        -1,
                        False)
            (id_seq1, seq1, id_seq2, seq2, score, alignments, num_alignments) = result
    # A TypeError is thrown when there are too few arguments (we are missing 1 argument)
    except TypeError:
        test = True
    assert test == True

if __name__ == '__main__':
    test_instance()
    test_example_success()
    test_example_invalid_format_fail()
    test_example_invalid_characters_fail()
    test_too_few_arguments()