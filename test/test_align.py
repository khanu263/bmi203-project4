# Importing Dependencies
import pytest
import numpy as np
from align import NeedlemanWunsch, read_fasta

def test_nw_alignment():
    """

    Write your unit test using test_seq1.fa and test_seq2.fa by asserting
    that you have correctly filled out the three alignment matrices. Use
    the BLOSUM62 matrix with a gap open penalty of -10 and a gap extension
    penalty of -1.

    """

    # Read in test sequences
    seq1, seq1_header = read_fasta("./data/test_seq1.fa")
    seq2, seq2_header = read_fasta("./data/test_seq2.fa")
    assert seq1_header == ">test seq 1", "seq1 header is different from expected."
    assert seq2_header == ">test seq 2", "seq2 header is different from expected."

    # Define expected matrices (from doing alignment by hand)
    expected_align = np.array([[      0, -np.inf, -np.inf, -np.inf],
                               [-np.inf,       5,     -11,     -13],
                               [-np.inf,     -12,       4,      -8],
                               [-np.inf,     -12,      -1,       5],
                               [-np.inf,     -14,      -6,       4]])
    expected_gapA = np.array([[    -10, -11, -12, -13],
                              [-np.inf, -22,  -6,  -7],
                              [-np.inf, -23, -17,  -7],
                              [-np.inf, -24, -18, -12],
                              [-np.inf, -25, -19, -17]])
    expected_gapB = np.array([[-10, -np.inf, -np.inf, -np.inf],
                              [-11,     -22,     -23,     -24],
                              [-12,      -6,     -17,     -18],
                              [-13,      -7,      -7,     -18],
                              [-14,      -8,      -8,      -6]])

    # Run alignment and make sure matrices are equal
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    _ = nw.align(seq1, seq2)
    assert np.array_equal(nw._align_matrix, expected_align), "align matrix does not match expected."
    assert np.array_equal(nw._gapA_matrix, expected_gapA), "gapA matrix does not match expected."
    assert np.array_equal(nw._gapB_matrix, expected_gapB), "gapB matrix does not match expected."

def test_nw_backtrace():
    """

    Write your unit test using test_seq3.fa and test_seq4.fa by asserting
    that the backtrace is correct. Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.

    """

    # Read in test sequences
    seq3, seq3_header = read_fasta("./data/test_seq3.fa")
    seq4, seq4_header = read_fasta("./data/test_seq4.fa")
    assert seq3_header == ">test seq 3", "seq3 header is different from expected."
    assert seq4_header == ">test seq 4", "seq4 header is different from expected."

    # Run alignment and check backtrace
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, seq3_align, seq4_align = nw.align(seq3, seq4)
    assert score == 17, "Alignment score is incorrect."
    assert seq3_align == "MAVHQLIRRP", "seq3 alignment is incorrect."
    assert seq4_align == "M---QLIRHP", "seq4 alignment is incorrect."
    