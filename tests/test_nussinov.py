import sys,os
import pytest
sys.path.append(os.path.abspath('.'))
from prakt.nussinov import NussinovBase
from nussinov_alg import Nussinov

def test_instance():
    """Check inheritance."""
    assert issubclass(Nussinov, NussinovBase)
    assert isinstance(Nussinov(), NussinovBase)


def test_example_success():
    """This calls the run method of the Nussinov program
    and tests if it works as expected (positive test)"""

    nu = Nussinov()
    seq_fasta = os.path.join('data','sequences','nuss.fasta')
    result = nu.run(seq_fasta)
    (N,dot_bracket) = result

    assert dot_bracket == "()((((()())....((((())().)))).)())......."


def test_example_invalid_format_fail():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for the failure is invalid file format: the first line does not start with >"""

    nu = Nussinov()
    seq_fasta = os.path.join('data','sequences','Invalid_formatRNA.fasta')
    with pytest.raises(SystemExit) as InvalidFileException:
        result = nu.run(seq_fasta)
        (N,dot_bracket) = result

        assert InvalidFileException.type == SystemExit
        assert InvalidFileException.code == 1

def test_example_invalid_characters_fail():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for failure is non-amino acid characters in file 2 (error code 12)"""

    nu = Nussinov()
    seq_fasta = os.path.join('data','sequences','Invalid_charactersRNA.fasta')
    with pytest.raises(SystemExit) as InvalidCharactersException:
        result = nu.run(seq_fasta)
        (N,dot_bracket) = result

        assert InvalidCharactersException.type == SystemExit
        assert InvalidCharactersException.code == 12

def test_too_few_arguments():
    """This function does a negative test: it checks if it fails when it is supposed to.
    The reason for failure is non-amino acid characters in file 2 (error code 12)"""

    nu = Nussinov()
    seq_fasta_1 = os.path.join('data','sequences','nuss.fasta')
    # test is a variable which becomes True when there are too few arguments
    test = False
    try:
        with pytest.raises(SystemExit) as TooFewArguments:
            result = nu.run()
            (N,dot_bracket) = result
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
