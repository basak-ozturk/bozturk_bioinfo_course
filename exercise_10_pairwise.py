from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import pandas as pd
# Parse the FASTA file and store sequences in a list
sequences = list(SeqIO.parse("Sequences_MSA_01.fasta", "fasta"))
aligner = Align.PairwiseAligner()
#print(sequences)
print(f"Number of sequences read: {len(sequences)}")
for seq_record in sequences:
    print(seq_record.id)
score=[]
print(f"The total number of sequences is: {len(sequences)}.")
for i in range(len(sequences)):
    for j in range(i+1, len(sequences)):
        seq1=sequences[i]
        seq2=sequences[j]
#         print(seq1)
#         print(seq2)
        
        alignments = aligner.align(seq1, seq2)
#        print(len(alignments))
        score.append(aligner.score(seq1.seq, seq2.seq))
        print(score)

# Example: Perform pairwise alignment between each pair of sequences
# for i in range(len(sequences)):
#     for j in range(i + 1, len(sequences)):
#         seq1 = sequences[i].seq
#         seq2 = sequences[j].seq
# 
#         # Perform global pairwise alignment using pairwise2
#         alignments = pairwise2.align.globalxx(seq1, seq2)
#         
#         # Access the best alignment (optional)
#         best_alignment = alignments[0]
# 
#         # Print or process alignment
#         print(f"Alignment between sequence {i+1} and sequence {j+1}:")
#         print(pairwise2.format_alignment(*best_alignment))
