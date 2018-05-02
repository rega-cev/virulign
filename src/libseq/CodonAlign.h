// This may look like C code, but it's really -*- C++ -*-
#ifndef CODON_ALIGN_H_
#define CODON_ALIGN_H_

#include <AlignmentAlgorithm.h>

/**
 * libseq namespace
 */
namespace seq {
/**
  * Thrown when alignment failed.
  */
class AlignmentError : public std::exception
{
public:
  AlignmentError(double ntScore, double codonScore,
		  const NTSequence& ntRef, const NTSequence& ntTarget,
		  const std::string& message = std::string("Alignment error."));
  virtual ~AlignmentError() throw();

  /** %Nucleotide alignment score.
   */
  double nucleotideAlignmentScore() const { return ntScore_; }

  /** Codon-based alignemnt score.
   */
  double codonAlignmentScore() const { return codonScore_; }

  /** %Nucleotide aligned reference sequence
   */
  const NTSequence& nucleotideAlignedRef() const { return ntRef_; }

  /** %Nucleotide aligned target sequence
   */
  const NTSequence& nucleotideAlignedTarget() const { return ntTarget_; }

  /** Error message
   */
  const std::string& message() const{ return message_; }

private:
  std::string message_;
  double ntScore_, codonScore_;
  NTSequence ntRef_, ntTarget_;
};

/**
 * Error thrown by CodonAlign when apparent frame shifts cannot be corrected.
 *
 * Details in CodonAlign.
 */
class FrameShiftError : public AlignmentError
{
public:
  FrameShiftError(double ntScore, double codonScore,
		  const NTSequence& ntRef, const NTSequence& ntTarget);
  ~FrameShiftError() throw();

  const char *what() const throw() { return "Frameshift error"; }
};


class CodonAlign {
public:
  /**
   * Constructor
   */
  CodonAlign(AlignmentAlgorithm* algorithm);

 /**
 * Perform codon-based alignment of nucleotide sequences.
 *
 * Two nucleotide sequences are pair-wise aligned, but so that gaps are
 * at codon boundaries. Optionally, frameshifts may be detected and corrected.
 *
 * The reference sequence must be of length a multiple of 3, and is assumed
 * to represent an Open Reading Frame (ORF).
 *
 * The procedure translates the target sequence in the 3 ORFs,
 * and for each ORF performs an amino-acid alignment against the translated
 * reference sequence. The best alignment is used to create the nucleotide
 * alignment.
 *
 * Then, the score of the codon aligned nucleotide alignment is computed, and
 * compared with a direct nucleotide alignment of both nucleotide sequences.
 * The codon alignment is accepted only if the difference is smaller than 100.
 * Otherwise, if maxFrameShifts > 0, the frameshift is searched, corrected
 * by inserting 1 or 2 'N' symbols in the target sequence, and repeating the
 * codon alignment. This is repeated for up to maxFrameShifts of times.
 *
 * The result is the nucleotide alignment score of the codon alignment, and
 * the number of frameshifts that have been corrected.
 *
 * @throws FrameShiftError when frameshifts could not be corrected, or
 *         the number of detected frameshifts exceeds maxFrameShifts.
 */
 std::pair<double, int>
 align(NTSequence& ref, NTSequence& target, int maxFrameShifts = 1);

private:
  bool haveGaps(const NTSequence& seq, int from, int to);
  double alignLikeAA(NTSequence& seq1, NTSequence& seq2, 
		     int ORF, 
		     const AASequence& seqAA1, const AASequence& seqAA2);
  bool noGapAt(const NTSequence& seq, unsigned int i) const;

  AlignmentAlgorithm* algorithm_;
};
}

#endif // CODON_ALIGN_H_
