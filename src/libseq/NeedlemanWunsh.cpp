#include "NeedlemanWunsh.h"

#include <algorithm>

namespace seq {

NeedlemanWunsh::NeedlemanWunsh(double gapOpenScore,
			       double gapExtensionScore,
			       double **ntWeightMatrix,
			       double **aaWeightMatrix)
{
  gapOpenScore_ = gapOpenScore;
  gapExtensionScore_ = gapExtensionScore;
  ntWeightMatrix_ = ntWeightMatrix;
  aaWeightMatrix_ = aaWeightMatrix;
}

/*
 * A straight-forward implementation of Neeldeman-Wunsh algorithm
 * for a pairwise global alignment, with the difference that a
 * gapOpenScore is not added at the beginning or end of the sequence
 * (like ClustalW does).
 */
template <typename Symbol>
double NeedlemanWunsh::needlemanWunshAlign(std::vector<Symbol>& seq1,
					   std::vector<Symbol>& seq2,
					   double** weightMatrix)
{
  /*
   * Remove gaps, and warn that we did.
   */
  bool foundGaps = false;
  for (unsigned i = 0; i < seq1.size(); ++i) {
    if (seq1[i] == Symbol::GAP) {
      if (!foundGaps) {
	std::cerr << "Warning: NeedlemanWunsh: sequence contained gaps? "
	             "Removed them." << std::endl;
	foundGaps = true;
      }
      seq1.erase(seq1.begin() + i);
      --i;
    }
  }

  for (unsigned i = 0; i < seq2.size(); ++i) {
    if (seq2[i] == Symbol::GAP) {
      if (!foundGaps) {
	std::cerr << "Warning: NeedlemanWunsh: sequence contained gaps? "
	             "Removed them." << std::endl;
	foundGaps = true;
      }
      seq2.erase(seq2.begin() + i);
      --i;
    }
  }

  const int seq1Size = seq1.size();
  const int seq2Size = seq2.size();

  double **dnTable = new double* [seq1Size+1];
  for (unsigned i = 0; i < seq1Size+1; ++i)
    dnTable[i] = new double[seq2Size+1];
  int    **gapsLengthTable = new int *[seq1Size+1];
  for (unsigned i = 0; i < seq1Size+1; ++i)
    gapsLengthTable[i] = new int[seq2Size+1]; // >0: horiz, <0: vert

  double edgeGapExtensionScore = 0;

  /*
   * compute table
   */
  dnTable[0][0] = 0;
  gapsLengthTable[0][0] = 0;
  for (unsigned i = 1; i < seq1Size+1; ++i) {
    dnTable[i][0] = dnTable[i-1][0] + edgeGapExtensionScore;
    gapsLengthTable[i][0] = gapsLengthTable[i-1][0] + 1;
  }
  for (unsigned j = 1; j < seq2Size+1; ++j) {
    dnTable[0][j] = dnTable[0][j-1] + edgeGapExtensionScore;
    gapsLengthTable[0][j] = gapsLengthTable[0][j-1] - 1;
  }

  for (unsigned i = 1; i < seq1Size+1; ++i) {
    for (unsigned j = 1; j < seq2Size+1; ++j) {

      double sextend
	= dnTable[i-1][j-1]
	+ weightMatrix[seq1[i-1].intRep()][seq2[j-1].intRep()];

      double ges = (j == seq2Size) ? edgeGapExtensionScore : gapExtensionScore_;

      double horizGapScore = ((gapsLengthTable[i-1][j] > 0) || (j == seq2Size)
			      ? ges : gapOpenScore_ + ges);
      double sgaphoriz
	= dnTable[i-1][j] + horizGapScore;

      ges = (i == seq1Size) ? edgeGapExtensionScore : gapExtensionScore_;

      double vertGapScore = (gapsLengthTable[i][j-1] < 0 || (i == seq1Size)
			     ? ges : gapOpenScore_ + ges);
      double sgapvert
	= dnTable[i][j-1] + vertGapScore;

      if ((sextend >= sgaphoriz) && (sextend >= sgapvert)) {
	dnTable[i][j] = sextend;
	gapsLengthTable[i][j] = 0;
      } else {
	if (sgaphoriz > sgapvert) {
	  dnTable[i][j] = sgaphoriz;
	  gapsLengthTable[i][j] = std::max(0, gapsLengthTable[i-1][j]) + 1;
	} else {
	  dnTable[i][j] = sgapvert;
	  gapsLengthTable[i][j] = std::min(0, gapsLengthTable[i][j-1]) - 1;
	}
      }
    }
  }

  /*
   * reconstruct best solution alignment.
   */
  int i = seq1Size+1, j = seq2Size+1;
  do {
    if (gapsLengthTable[i-1][j-1] == 0) {
      --i; --j;
    } else if (gapsLengthTable[i-1][j-1] > 0) {
      --i;
      seq2.insert(seq2.begin() + (j-1), Symbol::GAP);
    } else {
      --j;
      seq1.insert(seq1.begin() + (i-1), Symbol::GAP);
    }
  } while (i > 1 || j > 1);

  double score = dnTable[seq1Size][seq2Size];

  for (unsigned i = 0; i < seq1Size+1; ++i) {
    delete[] dnTable[i];
    delete[] gapsLengthTable[i];
  }
  delete[] dnTable;
  delete[] gapsLengthTable;

  return score;
}
  
double NeedlemanWunsh::align(NTSequence& seq1, NTSequence& seq2)
{
  return needlemanWunshAlign(seq1, seq2, ntWeightMatrix_);
}

double NeedlemanWunsh::align(AASequence& seq1, AASequence& seq2)
{
  return needlemanWunshAlign(seq1, seq2, aaWeightMatrix_);
}

double NeedlemanWunsh::computeAlignScore(const NTSequence& seq1, 
					 const NTSequence& seq2)
{
  double score = 0;
  int seq1GapLength = 0;
  int seq2GapLength = 0;

  bool seq1LeadingGap = true;
  bool seq2LeadingGap = true;

  double edgeGapExtensionScore = 0;
  
  for (unsigned i = 0; i < seq1.size(); ++i) {
    if (seq1[i] == Nucleotide::GAP) {
      ++seq1GapLength;
    } else {
      if (seq1GapLength) {
	    if (seq1LeadingGap)
	        score += seq1GapLength * edgeGapExtensionScore;
	    else
	    score += gapOpenScore_ + seq1GapLength * gapExtensionScore_;
      }
      seq1GapLength = 0;

      if (seq2[i] == Nucleotide::GAP) {
	++seq2GapLength;
      } else {
	if (seq2GapLength) {
	  if (seq2LeadingGap)
	    score += seq2GapLength * edgeGapExtensionScore;
	  else
	    score += gapOpenScore_ + seq2GapLength * gapExtensionScore_;
    }
	seq2GapLength = 0;

	score += ntWeightMatrix_[seq1[i].intRep()][seq2[i].intRep()];
      }
    }
  }

  score += seq1GapLength * edgeGapExtensionScore;
  score += seq2GapLength * edgeGapExtensionScore;

  return score;
}

}
