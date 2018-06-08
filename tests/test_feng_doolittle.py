import sys,os
import pytest
sys.path.append(os.path.abspath('.'))
from prakt.fd import FengDoolittleBase
from feng_doolittle import FengDoolittle

def test_instance():
    """Check inheritance."""
    assert issubclass(FengDoolittle, FengDoolittleBase)
    assert isinstance(FengDoolittle('blosum62','-1'), FengDoolittleBase)

def test_example_successs():
    """This calls the run method of the Feng-Doolittle program
    and tests if it works as expected."""

    fd = FengDoolittle('blosum62','-1')
    seq_fasta = os.path.join('data','sequences','mult.fasta')
    result = fd.run(seq_fasta,
                    'pam250',
                    -1,
                    'UPGMA')
    SOP, newick, newickNoDist = result

    assert newickNoDist == "(((C0,C1),C3),C2);"
    # Sum of pairs changes based on the random alignment obtained from Needleman-Wunsch, not possible to test.
    #assert SOP ==


def test_example_invalid_format_fail():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for the failure is invalid file format: the first line does not start with >"""

    fd = FengDoolittle('blosum62','-1')
    seq_fasta = os.path.join('data','sequences','Invalid_formatMult.fasta')
    with pytest.raises(SystemExit) as InvalidFileException:
        result = fd.run(seq_fasta,
                    'blosum62',
                    -1,
                    'UPGMA')
        SOP, newick = result
        assert InvalidFileException.type == SystemExit
        assert InvalidFileException.code == 1

def test_example_invalid_characters_fail():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for failure is non-amino acid characters in file 2 (error code 12)"""

    fd = FengDoolittle('blosum62','-1')
    seq_fasta = os.path.join('data','sequences','Invalid_charactersMult.fasta')
    with pytest.raises(SystemExit) as InvalidCharactersException:
        fd.run(seq_fasta,
                    'blosum62',
                    -1,
                    'UPGMA')
        SOP, newick = result

        assert InvalidCharactersException.type == SystemExit
        assert InvalidCharactersException.code == 12

def test_too_few_arguments():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for failure is non-amino acid characters in file 2 (error code 12)"""

    fd = FengDoolittle('blosum62','-1')
    seq_fasta = os.path.join('data','sequences','mult.fasta')
    # test is a variable which becomes True when there are too few arguments
    test = False
    try:
        with pytest.raises(SystemExit) as TooFewArguments:
            result = fd.run(seq_fasta,
                    'blosum62',
                    'UPGMA')
        SOP, newick = result
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
