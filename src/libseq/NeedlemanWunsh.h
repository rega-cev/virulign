// This may look like C code, but it's really -*- C++ -*-
#ifndef NEEDLEMAN_WUNSH_H_
#define NEEDLEMAN_WUNSH_H_

#include <AlignmentAlgorithm.h>

/**
 * libseq namespace
 */
namespace seq {

class NeedlemanWunsh : public AlignmentAlgorithm 
{
  public:
    NeedlemanWunsh(double gapOpenScore = -10,
		   double gapExtensionScore = -3.3,
		   double **ntWeightMatrix = 
		   AlignmentAlgorithm::IUB(),
		   double **aaWeightMatrix = 
		   AlignmentAlgorithm::BLOSUM30());
  /**
   * Pair-wise align two nucleotide sequences, using a modified
   * NeedleMan-Wunsh algorithm.
   *
   * The two sequences seq1 and seq2 are aligned in-place: gaps are inserted
   * according to a global alignment, and they will have equal length.
   *
   * The algorithm is NeedleMan-Wunsh, with two popular modifications:
   *  - there is a different cost for opening a gap or for extending a gap.
   *  - there is no gap open cost for a gap at the beginning or the end.
   */
  virtual double align(NTSequence& seq1, NTSequence& seq2);

  /**
   * Pair-wise align two amino acid sequences, using a modified
   * NeedleMan-Wunsh algorithm.
   *
   * The two sequences seq1 and seq2 are aligned in-place: gaps are inserted
   * according to a global alignment, and they will have equal length.
   *
   * The algorithm is NeedleMan-Wunsh, with two popular modifications:
   *  - there is a different cost for opening a gap or for extending a gap.
   *  - there is no gap open cost for a gap at the beginning or the end.
   */
  virtual double align(AASequence& seq1, AASequence& seq2);

  virtual double computeAlignScore(const NTSequence& seq1, 
				   const NTSequence& seq2);

private:
  double gapOpenScore_;
  double gapExtensionScore_;
  double **ntWeightMatrix_;
  double **aaWeightMatrix_;

  template <typename Symbol>
  double needlemanWunshAlign(std::vector<Symbol>& seq1,
			     std::vector<Symbol>& seq2,
			     double** weigthMatrix);
};

}

#endif // NEEDLEMAN_WUNSH_H_
