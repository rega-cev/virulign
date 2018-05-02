// This may look like C code, but it's really -*- C++ -*-
#ifndef AASEQUENCE_H_
#define AASEQUENCE_H_

#include <vector>
#include <string>
#include <iostream>

#include "NTSequence.h"
#include "AminoAcid.h"

namespace seq {

/**
 * An amino acid sequence.
 *
 * The sequence may have a name and a description.
 *
 * The sequence data is stored by publicly inheriting
 * std::vector<AminoAcid>, so you can use all
 * <a HREF="http://wwwasd.web.cern.ch/wwwasd/lhc++/RW/stdlibcr/vec_0251.htm">
 * std::vector</a> manipulations to access the amino acid data.
 */
class AASequence : public std::vector<AminoAcid>
{
public:
  /**
   * Create an empty amino acid sequence with emtpy name and empty
   * description.
   */
  AASequence();

  /**
   * Create an amino acid sequence of length size, filled with
   * AminoAcid::X, with empty name and emtpy description.
   */
  AASequence(unsigned size);

  /**
   * Create an amino acid sequence with given name and description, and
   * with the given sequence string. Each character in the sequence
   * string will be interpreted as an AminoAcid using the
   * AminoAcid::AminoAcid(char) constructor.
   */
  AASequence(const std::string name,
	     const std::string description,
	     const std::string aSeqString);

  /**
   * Create a nucleotide sequence with empty name and emtpy
   * description, and copy the sequence data from the range [first, last[.
   */
  AASequence(const const_iterator first, const const_iterator last);

  /**
   * Represent the sequence data as a string.
   */
  std::string asString() const;

  /**
   * Get the name.
   */
  std::string name() const { return name_; }

  /**
   * Get the description.
   */
  std::string description() const { return description_; }

  /**
   * Set the name.
   */
  void setName(std::string name) { name_ = name; }

  /**
   * Set the description.
   */
  void setDescription(std::string description) { description_ = description; }

  /**
   * Translate a nucleotide sequence to an amino acid sequence. The
   * nucleotide sequence must have a length that is a multiple of
   * three.
   *
   * The resulting amino acid sequence will contain an amino acid for
   * every triplet of nucleotides in the nucleotide sequence.  The
   * amino acid sequence will have the same name and description as
   * the nucleotide sequence.
   *
   * \sa translate(const NTSequence::const_iterator, const NTSequence::const_iterator), Codon::translate(const NTSequence::const_iterator)
   */
  static AASequence translate(const NTSequence& ntSequence);

  /**
   * Translate a nucleotide sequence, defined by the range begin to
   * end, to an amino acid sequence. The nucleotide sequence must have a
   * length that is a multiple of three.
   *
   * The resulting amino acid sequence will contain an amino acid for
   * every triplet of nucleotides in the nucleotide sequence, and will
   * have an empty name and empty description.
   *
   * \sa translate(const NTSequence&), Codon::translate(const NTSequence::const_iterator)
   */
  static AASequence translate(const NTSequence::const_iterator begin,
			      const NTSequence::const_iterator end);
private:
  std::string name_;
  std::string description_;
};

/**
 * Read an amino acid sequence in FASTA format from the given stream.
 */
extern std::istream& operator>>(std::istream& i, AASequence& sequence);

/**
 * Write an amino acid sequence to the given stream in FASTA format.
 */
extern std::ostream& operator<<(std::ostream& o, const AASequence& sequence);

};

#endif // AASEQUENCE_H_
