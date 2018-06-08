import sys,os
import pytest
sys.path.append(os.path.abspath('.'))
from prakt.xpgma import XpgmaBase
from uwpgma import Xpgma

def test_instance():
    """Check inheritance."""
    assert issubclass(Xpgma, XpgmaBase)
    assert isinstance(Xpgma(), XpgmaBase)

def test_example_success_UPGMA():
    """This calls the run method of the UWPGMA (UPGMA + WPGMA) program
    and tests if it works as expected for UPGMA with PAM(positive test)"""

    gma = Xpgma()
    seq_fasta = os.path.join('data','sequences','mult.fasta')
    result = gma.run(seq_fasta,
                    'pam250',
                    -1,
                    'UPGMA')
    newick,newickNoDistance,newickIds = result

    assert newickNoDistance == "(((C0,C1),C3),C2);"

def test_example_success_WPGMA():
    """This calls the run method of the UWPGMA (UPGMA + WPGMA) program
    and tests if it works as expected for WPGMA with BLOSUM(positive test)"""

    gma = Xpgma()
    seq_fasta = os.path.join('data','sequences','mult.fasta')
    result = gma.run(seq_fasta,
                    'blosum62',
                    -1,
                    'WPGMA')
    newick, newickNoDistance, newickIds = result

    assert newickNoDistance == "(((C0,C1),C3),C2);"

def test_example_invalid_format_fail():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for the failure is invalid file format: the first line does not start with >"""

    gma = Xpgma()
    seq_fasta = os.path.join('data','sequences','Invalid_formatMult.fasta')
    with pytest.raises(SystemExit) as InvalidFileException:
        result = gma.run(seq_fasta,
                    'blosum62',
                    -1,
                    'WPGMA')
        newick, newickNoDistance = result
        assert InvalidFileException.type == SystemExit
        assert InvalidFileException.code == 1

def test_example_invalid_characters_fail():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for failure is non-amino acid characters in file 2 (error code 12)"""

    gma = Xpgma()
    seq_fasta = os.path.join('data','sequences','Invalid_charactersMult.fasta')
    with pytest.raises(SystemExit) as InvalidCharactersException:
        result = gma.run(seq_fasta,
                    'blosum62',
                    -1,
                    'WPGMA')
        newick, newickNoDistance = result

        assert InvalidCharactersException.type == SystemExit
        assert InvalidCharactersException.code == 12

def test_too_few_arguments():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for failure is non-amino acid characters in file 2 (error code 12)"""

    gma = Xpgma()
    seq_fasta = os.path.join('data','sequences','mult.fasta')
    # test is a variable which becomes True when there are too few arguments
    test = False
    try:
        with pytest.raises(SystemExit) as TooFewArguments:
            result = gma.run(seq_fasta,
                    'blosum62',
                    'UPGMA')
        newick, newickNoDistance = result
    # A TypeError is thrown when there are too few arguments (we are missing 1 argument)
    except TypeError:
        test = True
    assert test == True

if __name__ == '__main__':
    test_instance()
    test_example_success_UPGMA()
    test_example_success_WPGMA()
    test_example_invalid_format_fail()
    test_example_invalid_characters_fail()
    test_too_few_arguments()
