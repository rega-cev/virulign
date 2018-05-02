// This may look like C code, but it's really -*- C++ -*-
#ifndef CODING_SEQUENCE_H_
#define CODING_SEQUENCE_H_

#include <set>
#include <iostream>

#include "NTSequence.h"
#include "AASequence.h"

namespace seq {

/**
 * A coding sequence represents a nucleotide sequence that codes for
 * an amino acid sequence (an oligo- or polypeptide).
 *
 * It is useful when one wants to track the effect of changes in the
 * nucleotide sequence for the amino acid sequence, and to investigate
 * properties of nucleotide mutations.
 */
class CodingSequence
{
 public:
  /**
   * Construct a coding sequence with empty nucleotide sequence.
   */
  CodingSequence();

  /**
   * Construct a coding sequence based on the given nucleotide
   * sequence. The sequence must be translatable as per
   * AASequence::translate(const NTSequence&).
   */
  CodingSequence(const NTSequence& aNtSequence);

  /**
   * Get the nucleotide sequence.
   */
  const NTSequence& ntSequence() const { return ntSequence_; }

  /**
   * Get the amino acid sequence.
   *
   * If needed, the amino acid sequence is updated to reflect changes
   * in the nucleotide sequence.
   */
  const AASequence& aaSequence() const;

  /**
   * Change a nucleotide at a given position in the nucleotide sequence to
   * a new value.
   */
  void changeNucleotide(int pos, Nucleotide value);

  /**
   * Investigate the effect of a nucleotide mutation on the amino acid
   * sequence. This returns both the old (oldAA) and new amino acid (newAA)
   * encoded by the mutation, as well as the position (return value).
   */
  int whatIfMutation(int pos, Nucleotide value,
		     AminoAcid& oldAA, AminoAcid& newAA) const;

  /**
   * Investigate whether a give nucleotide mutation is synonymous or
   * non-synonymous with respect to the amino acid sequence.
   */
  bool isSynonymousMutation(int pos, Nucleotide value) const;

  /**
   * Get the amino acid sequence possibilities, taking into account
   * all ambiguities
   */
  void allAASequences(std::vector<std::set<AminoAcid> >& result) const;

 protected:
  void updateAASequence() const;

 private:
  NTSequence         ntSequence_;
  mutable AASequence aaSequence_;

  bool               isDirty() const { return dirty_ != D_CLEAN; }

  mutable int        dirty_;

  static const int   D_CLEAN = -1;
  static const int   D_COMPLETE = -2;
};

  /**
   * Write an amino acid sequence with all possible ambiguities
   * to the stream.
   *
   * The format is e.g. TW{LM}YS
   */
extern void printAmbiguousAASequence(std::ostream& out,
				     const CodingSequence& cs);

};

#endif // CODING_SEQUENCE_H_
