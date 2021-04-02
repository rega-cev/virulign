#include "Utils.h"

#include <Codon.h>
#include <CodonAlign.h>
#include <NeedlemanWunsh.h>

#include "Alignment.h"

#include <algorithm>

Alignment Alignment::compute(const ReferenceSequence& ref,
			     const seq::NTSequence&   target,
			     seq::AlignmentAlgorithm* algorithm,
			     int maxFrameShifts)
{
  seq::CodonAlign codonAlign(algorithm);
  Alignment result(ref, target);

  for (unsigned j = 0; j < result.target.size(); ++j)
    if (result.target[j] == seq::Nucleotide::GAP) {
      result.target.erase(result.target.begin() + j);
      --j;
    }

  try {
    if (result.target.size() > 6) {
      std::pair<double, int> res
	= codonAlign.align(result.ref, result.target, maxFrameShifts);

      result.score = res.first;
      result.correctedFrameshifts = res.second;
      result.success = true;
    } else
      result.tooShort = true;
  } catch (seq::AlignmentError e) {
    result.failure = true;
    std::cerr << e.nucleotideAlignedTarget().name() << ": " << e.message()
      << " (scores nt: " << e.nucleotideAlignmentScore() << "; codon: "
      << e.codonAlignmentScore() << ")" << std::endl;
  }

  result.computeAlignedRanges(ref.size()/3);

  return result;
}

Alignment Alignment::given(const ReferenceSequence& ref,
			   const seq::NTSequence&   target)
{
  if (ref.size() != target.size()) {
    std::cerr << ref.name() << ".length: " << ref.size()
	      << ", " << target.name() << ".length: " << target.size()
	      << std::endl;

    assert(ref.size() == target.size());
  }

  Alignment result(ref, target);

  result.success = true;

  result.computeAlignedRanges(ref.size()/3);

  return result;
}

Alignment::Alignment(const ReferenceSequence& aref,
		     const seq::NTSequence&   atarget)
  : success(false),
    tooShort(false),
    failure(false),
    correctedFrameshifts(0),
    ref(aref),
    target(atarget)    
{ }

void Alignment::computeAlignedRanges(int referenceSequenceLength)
{
  for (unsigned r = 0; r < ref.regions().size(); ++r) {
    ReferenceSequence::Region& region = ref.regions()[r];

    int regionEnd = std::min(region.end(), referenceSequenceLength);

    if (success) {
      region.alignedBegin  = alignedPos(region.begin());
      region.alignedEnd    = alignedPos(regionEnd);
      region.targetBegin = firstPos(region.begin(), regionEnd);
      region.targetEnd   = lastPos(region.begin(), regionEnd);
    } else {
      region.alignedBegin = region.begin();
      region.alignedEnd   = regionEnd;
      region.targetBegin = ref.size();
      region.targetEnd   = -1;
    }
  }
}

int Alignment::alignedPos(int refPos) const
{
  int j = -1;
  for (unsigned i = 0; i < ref.size(); i += 3) {
    if (ref[i] != seq::Nucleotide::GAP)
      ++j;

    if (j == refPos)
      return i/3;
  }
  if (j == refPos - 1)
    return ref.size() / 3;
  else {
    std::cerr << refPos << " " << ref.size() << " " << j << std::endl;
    assert(false);
  }
}

int Alignment::firstPos(int begin, int end) const
{
  int refPos = -1;

  for (unsigned i = 0; i < ref.size(); i += 3) {
    if (ref[i] != seq::Nucleotide::GAP)
      ++refPos;

    if (refPos >= begin) {
      if (refPos >= end)
	return end;

      if (target[i] != seq::Nucleotide::GAP)
	return refPos;
    }
  }

  return end;
}

int Alignment::lastPos(int begin, int end) const
{
  int refPos = -1;
  int lastPos = -1;

  for (unsigned i = 0; i < ref.size(); i += 3) {
    if (ref[i] != seq::Nucleotide::GAP)
      ++refPos;

    if (refPos >= begin) {
      if (refPos >= end)
	return lastPos;

      if (target[i+2] != seq::Nucleotide::GAP)
	lastPos = refPos;
    }
  }

  return lastPos;
}

std::pair<bool, int>
Alignment::findAminoAcid(const ReferenceSequence::Region& region,
			 int posInRegion, int insertion) const
{
  bool withinTarget
    = ((region.targetBegin < region.targetEnd)
       && posInRegion >= region.targetBegin - region.begin() + 1
       && posInRegion <= region.targetEnd   - region.begin() + 1);

  int pos = 0;
  int gap = 0;

  for (int i = region.alignedBegin; i < region.alignedEnd; ++i) {
    if (ref[i*3] != seq::Nucleotide::GAP) {
      ++pos;
      gap = 0;
    } else
      ++gap;

    if (pos == posInRegion
	&& gap == insertion
	&& (!withinTarget || (target[i*3] != seq::Nucleotide::GAP))) {
      return std::make_pair(withinTarget, i);
    } else if (pos > posInRegion) {
      return std::make_pair(withinTarget, -1);
    }
  }

  assert(false);
  return std::make_pair(false, 0);
}

std::string Alignment::mutations(const ReferenceSequence::Region& region) const
{
  std::string result;
  int fp    = region.targetBegin;
  int lp    = region.targetEnd;

  if (fp >= lp)
    return result;

  int refPos = -1;

  for (unsigned i = 0; i < ref.size(); i += 3) {
    if (ref[i] != seq::Nucleotide::GAP)
      ++refPos;

    if (refPos >= fp) {
      if (refPos > lp)
	return result;

      seq::AminoAcid refAA = seq::Codon::translate(ref.begin() + i);
      std::set<seq::AminoAcid> targetAAs
	= seq::Codon::translateAll(target.begin() + i);

      if (((targetAAs.size() > 1)
	   || (*targetAAs.begin() != refAA))
	  && (*targetAAs.begin() != seq::AminoAcid::GAP)) {

	if (!result.empty())
	  result += ' ';

	result += refAA.toChar()
	  + to_string(refPos - region.begin() + 1);

	for (std::set<seq::AminoAcid>::const_iterator k = targetAAs.begin();
	     k != targetAAs.end(); ++k)
	  result += k->toChar();

      }
    }
  }

  return result;
}

std::string Alignment::
codonMutations(const ReferenceSequence::Region& region,
	       int& start,
	       int& end) const
{
  std::string result;
  int fp    = region.begin();
  int lp    = region.end() - 1;

  start = -1;
  end = -1;

  if (fp >= lp)
    return result;

  int refPos = -1;

  for (unsigned i = 0; i < ref.size(); i += 3) {
    if (ref[i] != seq::Nucleotide::GAP)
      ++refPos;

    int pos = refPos - region.begin() + 1;

    if (refPos >= fp) {
      if (refPos > lp)
	return result;

      if (target[i] == seq::Nucleotide::GAP &&
	  target[i + 1] == seq::Nucleotide::GAP &&
	  target[i + 2] == seq::Nucleotide::GAP &&
	  (refPos < region.targetBegin || refPos > region.targetEnd))
	continue;

      if (refPos == region.targetEnd && 
          ref[i] == seq::Nucleotide::GAP && 
          ref[i + 1] == seq::Nucleotide::GAP && 
          ref[i + 2] == seq::Nucleotide::GAP)
        continue;

      //skip incomplete begin codon
      if(refPos == region.targetBegin-1 &&
          target[i] == seq::Nucleotide::GAP)
        continue;

      //skip incomplete end codon
      if(refPos == region.targetEnd+1 &&
          target[i + 2] == seq::Nucleotide::GAP)
        continue;

      if (start == -1)
	start = pos;
      end = pos;

      bool mutation;
      mutation = ref[i] != target[i] ||
                 ref[i + 1] != target[i + 1] ||
                 ref[i + 2] != target[i + 2];

      if(mutation) {
	if (!result.empty())
	  result += ' ';

        seq::AminoAcid refAA = seq::Codon::translate(ref.begin() + i);
        std::set<seq::AminoAcid> targetAAs = seq::Codon::translateAll(target.begin() + i);

        result += refAA.toChar()
               + to_string(refPos - region.begin() + 1);

        for (std::set<seq::AminoAcid>::const_iterator k = targetAAs.begin(); k != targetAAs.end(); ++k)
          result += k->toChar();
        result += ';';

	result += ref[i].toChar();
	result += ref[i+1].toChar(); 
	result += ref[i+2].toChar();
	result += to_string(pos);

	result += target[i].toChar();
	result += target[i + 1].toChar();
	result += target[i + 2].toChar();
      }
    }
  }

  return result;
}

Alignment::Alignment() {
}
