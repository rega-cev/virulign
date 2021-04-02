// This may look like C code, but it's really -*- C++ -*-
#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "ReferenceSequence.h"
#include <AlignmentAlgorithm.h>

class IsolateMutation;

class Alignment
{
public:
  bool   success;
  bool   tooShort;
  bool   failure;
  int    correctedFrameshifts;
  double score;

  ReferenceSequence ref;
  seq::NTSequence   target;

  std::string 
  mutations(const ReferenceSequence::Region& region) const;
  std::string 
  codonMutations(const ReferenceSequence::Region& region,
		 int& start,
		 int& end) const;
  void isolateMutations(const ReferenceSequence::Region& regioni, std::vector<IsolateMutation>& mutations) const;

  /*! \brief Return the amino acid position of the given mutation, if there
   *         is information on that mutation in the alignment
   *
   * The first result (bool) indicates if the target sequence contains
   * the mutation.
   *
   * The second result is the amino acid position. If the mutation is
   * an insertion which is not contained in the sequence, this value is
   * -1.
   */
  std::pair<bool, int> findAminoAcid(const ReferenceSequence::Region& region,
				     int positionInRegion, int insertion)
    const;

  static Alignment compute(const ReferenceSequence& ref,
			   const seq::NTSequence& target,
			   seq::AlignmentAlgorithm* algorithm,
			   int maxFrameShifts = 5);

  static Alignment given(const ReferenceSequence& ref,
			 const seq::NTSequence& target);

  void revert(const IsolateMutation& mutation);

  Alignment();

private:
  Alignment(const ReferenceSequence& aref,
	    const seq::NTSequence&   atarget);

  void     computeAlignedRanges(int referenceSequenceLength);
  int      alignedPos(int refPos) const;
  int      firstPos(int begin, int end) const;
  int      lastPos(int begin, int end) const;
};

#endif // ALIGNMENT_H_
