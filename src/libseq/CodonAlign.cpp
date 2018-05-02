#include "CodonAlign.h"

#include <algorithm>

namespace seq {

CodonAlign::CodonAlign(AlignmentAlgorithm* algorithm)
{ 
  algorithm_ = algorithm;
}

double CodonAlign::alignLikeAA(NTSequence& seq1, 
			       NTSequence& seq2, 
			       int ORF, 
			       const AASequence& seqAA1,
			       const AASequence& seqAA2)
{
  NTSequence seq2ORFLead(seq2.begin(), seq2.begin() + ORF);
  seq2.erase(seq2.begin(), seq2.begin() + ORF);
  int aaLength = seq2.size() / 3;
  NTSequence seq2ORFEnd(seq2.begin() + aaLength*3, seq2.end());
  seq2.erase(seq2.begin() + aaLength*3, seq2.end());

  int firstNonGap = -1;
  int lastNonGap = -1;

  for (unsigned i = 0; i < seqAA1.size(); ++i) {
    if (seqAA1[i] == AminoAcid::GAP && noGapAt(seq1, i)) {
      if (i*3 < seq1.size())
	seq1.insert(seq1.begin() + (i*3), 3, Nucleotide::GAP);
      else
	seq1.insert(seq1.end(), 3, Nucleotide::GAP);
    }

    if (seqAA2[i] == AminoAcid::GAP && noGapAt(seq2, i)) {
      if (i*3 < seq2.size())
	seq2.insert(seq2.begin() + (i*3), 3, Nucleotide::GAP);
      else
	seq2.insert(seq2.end(), 3, Nucleotide::GAP);
    } else {
      if (firstNonGap == -1)
	firstNonGap = i*3;
      lastNonGap = i*3 + 3;
    }
  }

  for (int i = 0; i < (int)seq2ORFLead.size(); ++i)
    if ((firstNonGap - (int)seq2ORFLead.size() + i) >= 0)
      seq2[firstNonGap - (int)seq2ORFLead.size() + i] = seq2ORFLead[i];

  for (unsigned i = 0; i < seq2ORFEnd.size(); ++i)
    if (lastNonGap + i < seq2.size())
      seq2[lastNonGap + i] = seq2ORFEnd[i];

  return algorithm_->computeAlignScore(seq1, seq2);
}

bool CodonAlign::noGapAt(const NTSequence& seq, unsigned int i) const
{
  if((i * 3) == seq.size())
    return true;
  else
    return seq[i * 3] != Nucleotide::GAP
        && seq[(i * 3) + 1] != Nucleotide::GAP
        && seq[(i * 3) + 2] != Nucleotide::GAP;
}

bool CodonAlign::haveGaps(const NTSequence& seq, int from, int to)
{
  for (unsigned i = std::max(from, 0); i < std::min((int)seq.size(), to); ++i)
    if (seq[i] == Nucleotide::GAP)
      return true;

  return false;
}

std::pair<double, int>
CodonAlign::align(NTSequence& ref, NTSequence& target, int maxFrameShifts)
{
  /*
   * 1. translate the reference sequence
   * 2. for every open reading frame:
   *   - translate the target sequence
   *   - perform the alignment
   * 3. take the alignment with best score and align nucleotide
   *    sequence like amino acid sequence
   * 4. compute nucleotide alignment score
   * 5. make nucleotide sequence alignment, compare score, if difference
   *    too big then correct the frame shift and repeat.
   */
  AASequence refAA = AASequence::translate(ref);

  NTSequence refNTAligned = ref;
  NTSequence targetNTAligned = target;
  double ntScore = algorithm_->align(refNTAligned, targetNTAligned);

  if(ntScore < 200)
    throw AlignmentError(ntScore,0,refNTAligned,targetNTAligned);

  int bestFrameShift = -1;
  double bestScore = -1E10;
  AASequence bestRefAA;
  AASequence bestTargetAA;

  for (unsigned i = 0; i < 3; ++i) {
    int last = i + ((target.size() - i) / 3) * 3;
    AASequence targetAA
      = AASequence::translate(target.begin() + i, target.begin() + last);

    AASequence refCopyAA = refAA;
    double score = algorithm_->align(refCopyAA, targetAA);

    if (score > bestScore) {
      bestFrameShift = i;
      bestScore = score;
      bestRefAA = refCopyAA;
      bestTargetAA = targetAA;
    }
  }

  NTSequence refCodonAligned = ref;
  NTSequence targetCodonAligned = target;

  double ntCodonScore = alignLikeAA(refCodonAligned,
				    targetCodonAligned,
				    bestFrameShift,
				    bestRefAA,
				    bestTargetAA);


  if (ntScore - ntCodonScore > 100) {
    /*
     * a possible frameshift
     */
    if (maxFrameShifts) {
      /*
       * try to fix: walk through the nucleotide alignment, and find
       * an "isolated" gap that is not of size multiple of 3.
       */
      const int BOUNDARY=10;
      int seq2pos = 0;
      int refGapStart = 0;
      int targetGapStart = 0;
      bool fixed = false;

      for (unsigned i = 0; i < refNTAligned.size(); ++i) {
	if (refNTAligned[i] == Nucleotide::GAP) {
	  if (refGapStart == -1)
	    refGapStart = i;
	} else { 
	  if (refGapStart > 0) {
	    int refGapStop = i;

	    if ((refGapStop - refGapStart) % 3) {
	      /*
	       * check it is isolated: no gaps in either sequence around
	       * this gap
	       */
	      if (haveGaps(refNTAligned,
			   refGapStart - BOUNDARY, refGapStart)
		  || haveGaps(refNTAligned,
			      refGapStop, refGapStop + BOUNDARY)
		  || haveGaps(targetNTAligned,
			      refGapStart - BOUNDARY, refGapStart)
		  || haveGaps(targetNTAligned,
			      refGapStop, refGapStop + BOUNDARY)) {
		/*
		 * not isolated: skip this gap.
		 */
	      } else {
		/*
		 * fix it !
		 */
		target.insert(target.begin() + seq2pos,
			      3 - (refGapStop - refGapStart) % 3,
			      Nucleotide::N);
		fixed = true;
		break;		
	      }
	    }
	  }

	  refGapStart = -1;
	}

	if (targetNTAligned[i] == Nucleotide::GAP) {
	  if (targetGapStart == -1)
	    targetGapStart = i;
	} else {
	  if (targetGapStart > 0) {
	    int targetGapStop = i;

	    if ((targetGapStop - targetGapStart) % 3) {
	      /*
	       * check it is isolated: no gaps in either sequence around
	       * this gap
	       */
	      if (haveGaps(refNTAligned,
			   targetGapStart - BOUNDARY, targetGapStart)
		  || haveGaps(refNTAligned, targetGapStop,
			      targetGapStop + BOUNDARY)
		  || haveGaps(targetNTAligned,
			      targetGapStart - BOUNDARY, targetGapStart)
		  || haveGaps(targetNTAligned,
			      targetGapStop, targetGapStop + BOUNDARY)) {
		/*
		 * not isolated: skip this gap.
		 */
	      } else {
		/*
		 * fix it !
		 */
		target.insert(target.begin() + seq2pos,
			      (targetGapStop - targetGapStart) % 3,
			      Nucleotide::N);
		fixed = true;
		break;
	      }
	    }
	  }

	  targetGapStart = -1;
	  ++seq2pos;
	}
      }

      if (!fixed)
	throw FrameShiftError(ntScore, ntCodonScore,
			      refNTAligned, targetNTAligned);
      else {
	std::pair<double, int> result
	  = align(ref, target, maxFrameShifts - 1);
	++result.second;
	return result;
      }
    } else {
      throw FrameShiftError(ntScore, ntCodonScore,
			    refNTAligned, targetNTAligned);
    }
  } else {
    ref = refCodonAligned;
    target = targetCodonAligned;
/*
    std::cerr << "Scores: " << ntScore << " " << ntCodonScore << " " << bestScore << std::endl;
    std::cerr << refNTAligned.asString() << std::endl;
    std::cerr << targetNTAligned.asString() << std::endl;
    std::cerr << refCodonAligned.asString() << std::endl;
    std::cerr << targetCodonAligned.asString() << std::endl;
    std::cerr << bestRefAA.asString() << std::endl;
    std::cerr << bestTargetAA.asString() << std::endl;
*/
    return std::make_pair(ntCodonScore, 0);
  }
}

AlignmentError::AlignmentError(double ntScore, double codonScore,
				 const NTSequence& ntRef,
				 const NTSequence& ntTarget,
				 const std::string& message)
  :ntScore_(ntScore),codonScore_(codonScore),
   ntRef_(ntRef),ntTarget_(ntTarget),
   message_(message)
{ }

AlignmentError::~AlignmentError() throw()
{ }


FrameShiftError::FrameShiftError(double ntScore, double codonScore,
				 const NTSequence& ntRef,
				 const NTSequence& ntTarget)
  :AlignmentError(ntScore,codonScore,ntRef,ntTarget,std::string("Frameshift error"))
{ }

FrameShiftError::~FrameShiftError() throw()
{ }

};

