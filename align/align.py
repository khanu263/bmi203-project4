# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """

        This function aligns the two given sequences using the Needleman-Wunsch
        algorithm, as per the scoring matrix and gap penalties specified when
        creating the aligner.

        Parameters:
            seqA: str
                One of the two sequences to align.
            seqB: str
                The other sequence to align.

        Returns:
            alignment_score: float
                The final alignment score.
            seqA_align: str
                Aligned representation of seqA, with '-' representing gaps.
            seqB_align: str
                Aligned representation of seqB, with '-' representing gaps.

        """
        
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        ### Implement the global sequence alignment here ###
        
        # Initialize first element / row / column of score matrices
        self._align_matrix[0, 0] = 0
        self._gapA_matrix[0, :] = self.gap_open + self.gap_extend * np.array(range(len(seqB) + 1))
        self._gapB_matrix[:, 0] = self.gap_open + self.gap_extend * np.array(range(len(seqA) + 1))

        # Initalize backtrace matrices (highroad alignment)
        # 0: we came from gapA
        # 1: we came from a match
        # 2: we came from gapB
        self._back_A[0, 1:] = 0
        self._back_B[1:, 0] = 2

        # Go through remaining matrix elements (since seqB defines columns, j comes first in the loop)
        for j in range(1, len(seqB) + 1):
            for i in range(1, len(seqA) + 1):

                # Compute scores for the three possibilities for all three matrices
                match = self.sub_dict[(seqA[i - 1], seqB[j - 1])] + np.array([self._gapA_matrix[i - 1, j - 1], self._align_matrix[i - 1, j - 1], self._gapB_matrix[i - 1, j - 1]])
                gapA = self.gap_extend + np.array([self._gapA_matrix[i, j - 1], self.gap_open + self._align_matrix[i, j - 1], self.gap_open + self._gapB_matrix[i, j - 1]])
                gapB = self.gap_extend + np.array([self.gap_open + self._gapA_matrix[i - 1, j], self.gap_open + self._align_matrix[i - 1, j], self._gapB_matrix[i - 1, j]])

                # Pick the maximums to fill the current elements
                self._align_matrix[i, j] = np.max(match)
                self._gapA_matrix[i, j] = np.max(gapA)
                self._gapB_matrix[i, j] = np.max(gapB)

                # Update backtrace matrices with argmax
                self._back[i, j] = np.argmax(match)
                self._back_A[i, j] = np.argmax(gapA)
                self._back_B[i, j] = np.argmax(gapB)

        # Return results of backtrace
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """

        Backtrace the alignment using the highroad heuristic. In the backtrace
        matrices, a 0 means the score came from a gap in A, a 1 means the 
        score came from a match, and a 2 means the score came from a gap in B.

        """

        # Initialize tracking variables (we start at the end and work backwards)
        i, j = len(self._seqA), len(self._seqB)
        end = (self._gapA_matrix[i, j], self._align_matrix[i, j], self._gapB_matrix[i, j])
        self.alignment_score = np.max(end)
        trace = np.argmax(end)

        # Keep going until we hit the beginning of both sequences
        while not (i == 0 and j == 0):

            # There was a gap in A
            if trace == 0:
                self.seqA_align = "-" + self.seqA_align
                self.seqB_align = self._seqB[j - 1] + self.seqB_align
                trace = self._back_A[i, j]
                j -= 1

            # There was a match
            elif trace == 1:
                self.seqA_align = self._seqA[i - 1] + self.seqA_align
                self.seqB_align = self._seqB[j - 1] + self.seqB_align
                trace = self._back[i, j]
                i -= 1
                j -= 1

            # There was a gap in B
            if trace == 2:
                self.seqA_align = self._seqA[i - 1] + self.seqA_align
                self.seqB_align = "-" + self.seqB_align
                trace = self._back_B[i, j]
                i -= 1

        # Return score and strings
        return self.alignment_score, self.seqA_align, self.seqB_align


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
