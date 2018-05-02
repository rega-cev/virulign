// This may look like C code, but it's really -*- C++ -*-
#ifndef ALIGNMENT_ALGORITHM_H_
#define ALIGNMENT_ALGORITHM_H_

#include <NTSequence.h>
#include <AASequence.h>

/**
 * libseq namespace
 */
namespace seq {

class AlignmentAlgorithm {
  public:
    /**
     * Pair-wise align two nucleotide sequences.
     *
     * The two sequences seq1 and seq2 are aligned in-place: gaps are inserted
     * according to a global alignment, and they will have equal length.
     */
    virtual double align(NTSequence& seq1, NTSequence& seq2) = 0;

    /**
     * Pair-wise align two amino acid sequences.
     *
     * The two sequences seq1 and seq2 are aligned in-place: gaps are inserted
     * according to a global alignment, and they will have equal length.
     */
    virtual double align(AASequence& seq1, AASequence& seq2) = 0;

    virtual double computeAlignScore(const NTSequence& seq1, 
				     const NTSequence& seq2) = 0;

    /**
     * Similarity weights matrix for nucleotides.
     *
     * Compares also IUB ambiuguity codes, and is the matrix used by BLAST.
     *
     * Taken from: ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4
     */
    static double** IUB();

    /**
     * Similarity weights matrix for amino acids.
     *
     * This is from the famous BLOSUM series of weight matrices, the one
     * that is the default use by ClustalX.
     *
     * From: ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM30
     */
    static double** BLOSUM30();
  };

}

#endif // ALIGNMENT_ALGORITHM_H_
