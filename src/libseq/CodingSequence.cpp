#include <ctype.h>

#include "Codon.h"
#include "CodingSequence.h"

namespace seq {

CodingSequence::CodingSequence()
  : ntSequence_(),
    aaSequence_(),
    dirty_(D_COMPLETE)
{ }

CodingSequence::CodingSequence(const NTSequence& aNtSequence)
  : ntSequence_(aNtSequence),
    aaSequence_(aNtSequence.size() / 3),
    dirty_(D_COMPLETE)
{ }

const AASequence& CodingSequence::aaSequence() const
{
  if (isDirty())
    updateAASequence();

  return aaSequence_;
}

void CodingSequence::changeNucleotide(int pos, Nucleotide value)
{
  // a small effort to avoid to avoid retranslation of the whole AA sequence.
  if (isDirty() && (dirty_ != D_COMPLETE))
    updateAASequence();

  ntSequence_[pos] = value;

  if (isDirty())
    dirty_ = D_COMPLETE;
  else
    dirty_ = pos;
}

int CodingSequence::whatIfMutation(int pos, Nucleotide value,
				   AminoAcid& oldAA,
				   AminoAcid& newAA) const
{
  if (isDirty())
    updateAASequence();  

  const int aaPos = pos / 3;
  const int codonPos = pos % 3;

  NTSequence newcodon(ntSequence_.begin() + aaPos * 3,
		      ntSequence_.begin() + (aaPos * 3 + 3));
  newcodon[codonPos] = value;

  oldAA = aaSequence_[aaPos];
  newAA = Codon::translate(newcodon.begin());

  return aaPos;
}

bool CodingSequence::isSynonymousMutation(int pos, Nucleotide value) const
{
  AminoAcid oldAA, newAA;

  whatIfMutation(pos, value, oldAA, newAA);

  return (oldAA == newAA);
}

void CodingSequence::updateAASequence() const
{
  if (dirty_ == D_COMPLETE) {
    aaSequence_ = AASequence::translate(ntSequence_);
  } else {
    dirty_ /= 3;
    aaSequence_[dirty_]
      = Codon::translate(ntSequence_.begin() + (dirty_ * 3)); 
  }

  dirty_ = D_CLEAN;
}

void CodingSequence::allAASequences(std::vector<std::set<AminoAcid> >& result)
  const
{
  for (unsigned i = 0; i < ntSequence_.size(); i += 3) {
    result.push_back(Codon::translateAll(ntSequence_.begin() + i));
  }
};

extern void printAmbiguousAASequence(std::ostream& out,
				     const CodingSequence& cs)
{
  std::vector<std::set<AminoAcid> > aas;
  cs.allAASequences(aas);

  for (unsigned i = 0; i < aas.size(); ++i) {
    if (aas[i].size() > 1)
      out << "{";
    for (std::set<AminoAcid>::const_iterator j = aas[i].begin();
	 j != aas[i].end(); ++j)
      out << *j;
    if (aas[i].size() > 1)
      out << "}";
  }
}

}
