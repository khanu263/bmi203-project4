# Project 4

Needleman-Wunsch algorithm.

# Assignment Overview

The purpose of this assigment is to have you implement the Needleman-Wunsch global pairwise sequence alignment algorithm (dynamic programming). See this [video](https://www.youtube.com/watch?v=NqYY0PJbD3s) for a walk through of the algorithm implementation. NOTE: This video misses an important step. In the video the X and Y gap matrices assume only two possible transitions when in reality there are three. The missing possibility is the chance that immediately go from a gap in one sequence to a gap in a second sequence. Please follow Lee's Lecture in terms of this approach. Please ask the TAs if you have any questions.

# Assignment Tasks

## Coding Assessment

**Note: All modules you need have already been imported.**

* [TODO] Complete the `NeedlemanWunsch.align` method found in the align/align.py 
* [TODO] Complete the `NeedlemanWunsch._backtrace` method found in align/align.py
* [TODO] Complete the `main` function in main.py to 
    1. align all provided species BRD2 sequences to the human BRD2 sequence and print the species in order of most similar to least similar with respect to human BRD2.
    2. print the alignment scores corresponding to each species alignment to the human BRD2 sequence.

## Software Development Assessment

### Unit Tests

* [TODO] Complete the `test_nw_alignment` function in test/test_align.py to test for proper matrix filling in your `NeedlemanWunsch.align` method.
* [TODO] Complete the `test_nw_backtrace` function in test/test_align.py to test for proper backtracing in your `NeedlemanWunsch._backtrace` method.

Note: test_seq3.fa and test_seq4.fa should have an alignment score of **17** and an alignment of:

```
MAVHQLIRRP	
M---QLIRHP
```

### Packaging
* [TODO] make .toml file with flit and ensure that your package can be installed with pip

### Automate Testing with [Github Actions](https://docs.github.com/en/actions)
See blogposts below on helping set up Github actions with pytest:
* [TODO] install and run pytest via GitHub Actions
* [TODO] add build status badge to README

# Getting Started
To get started you will need to fork this repository onto your own Github account. Work on the codebase from your own repo and commit changes. 

The following packages will be needed:
* numpy
* pytest

# Completing the assignment
Make sure to push all your code to Github, ensure that your unit tests are correct, and submit a link to your Github through the Google classroom assignment.

# Grading
## Code (6 points)
* Pairwise global alignment works properly (6)
    * Correct implementation of Needleman-Wunsch algorithm (4)
    * Produces correct order of species in main.py (1) 
    * Produces correct NW alignment scores in main.py (1)

## Unit tests (3 points)
* `test_nw_alignment` function properly checks that matrices are filled in correctly for alignment of test_seq1.fa and test_seq2.fa (1.5)
* `test_nw_backtrace` function properly checks that backtrace works correctly (1.5)

## Style (1 points)
* Readable code with clear comments and method descriptions (1)

