# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """

    # Read in sequences
    hs_seq, _ = read_fasta("./data/Homo_sapiens_BRD2.fa")
    species = ["Balaeniceps rex", "Gallus gallus", "Mus musculus", "Tursiops truncatus"]
    seqs = [read_fasta(f"./data/{s.replace(' ', '_')}_BRD2.fa")[0] for s in species]

    # Run alignments
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    scores = [nw.align(seq, hs_seq)[0] for seq in seqs]

    # Sort and print results
    scores, species = zip(*sorted(zip(scores, species), reverse = True))
    print("From most similar to least similar:")
    print("\n".join(f"{species[i]} (score: {scores[i]})" for i in range(4)))

if __name__ == "__main__":
    main()
