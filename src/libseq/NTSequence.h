// This may look like C code, but it's really -*- C++ -*-
#ifndef NTSEQUENCE_H_
#define NTSEQUENCE_H_

#include <vector>
#include <string>
#include <iostream>
#include <set>

#include "ParseException.h"
#include "Nucleotide.h"

namespace seq {

/**
 * A nucleotide sequence.
 *
 * The sequence may have a name and a description.
 *
 * The sequence data is stored by publicly inheriting
 * std::vector<Nucleotide>, so you can use all
 * <a HREF="http://wwwasd.web.cern.ch/wwwasd/lhc++/RW/stdlibcr/vec_0251.htm">
 * std::vector</a> manipulations to access the nucleotide data.
 */
class NTSequence : public std::vector<Nucleotide>
{
public:
  /**
   * Create an empty nucleotide sequence with emtpy name
   * and empty description.
   */
  NTSequence();
 
  /**
   * Create a nucleotide sequence of length size, filled with
   * Nucleotide::N, with empty name and emtpy description.
   */
  NTSequence(unsigned size);

  /**
   * Create a nucleotide sequence with given name and description, and
   * with the given sequence string. Each character in the sequence
   * string will be interpreted as a Nucleotide using the
   * Nucleotide::Nucleotide(char) constructor.
   *
   * If sampleAmbiguities = true, then sampleAmbiguities() is
   * performed during construction.
   *
   * \sa sampleAmbiguities()
   */
  NTSequence(const std::string name,
	     const std::string description,
	     const std::string aSeqString,
	     bool sampleAmbiguities = false);

  /**
   * Create a nucleotide sequence with empty name and emtpy
   * description, and copy the sequence data from the range [first, last[.
   */
  NTSequence(const const_iterator first,
	     const const_iterator last);

  /**
   * Remove ambiguity nucleotide symbols by replacing them by sampling
   * a random non-ambiguous nucleotide that is represented by the
   * ambiguity symbol.
   *
   * \sa Nucleotide::sampleAmbiguity()
   */
  void sampleAmbiguities();

  NTSequence reverseComplement() const;

  /**
   * Add all the possible non-ambiguous sequences possibly represented by
   * this sequence to result.
   */
  void nonAmbiguousSequences(std::vector<NTSequence>& result) const;

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

private:
  std::string name_;
  std::string description_;

  void iterateNonAmbiguous(const NTSequence& head,
			   std::vector<NTSequence>& result) const;
};

/**
 * Write a set of sequences to Stockholm format
 */
extern void writeStockholm(std::ostream& o,
                           const std::vector<NTSequence>& sequences,
                           int length=10000, int labelsize=0,
                           int seqsize=0, int pos=0);

/**
 * Read a nucleotide sequence in FASTA format from the given stream.
 */
extern std::istream& operator>>(std::istream& i, NTSequence& sequence);

/**
 * Write a nucleotide sequence to the given stream in FASTA format.
 */
extern std::ostream& operator<<(std::ostream& o, const NTSequence& sequence);

};

#endif // NTSEQUENCE_H_
