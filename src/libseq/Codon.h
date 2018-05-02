// This may look like C code, but it's really -*- C++ -*-
#ifndef CODON_H_
#define CODON_H_

#include <string>
#include <set>

#include "NTSequence.h"
#include "AminoAcid.h"

namespace seq {

/**
 * Utility class that defines the genetic code.
 */
class Codon
{
public:
  /**
   * Translate a nucleotide triplet (given by the range starting and
   * the indicated start point in a NTSequence) into an AminoAcid.
   *
   * If the triplet is three gaps, then the result is AminoAcid::GAP.
   * If the triplet contains ambiguity codes or gaps, then the result
   * is AminoAcid::X. Otherwise, the result is the translated amino
   * acid.
   */
  static AminoAcid translate(const NTSequence::const_iterator triplet);

  static std::set<AminoAcid>
     translateAll(const NTSequence::const_iterator triplet);

  static std::set<NTSequence> codonsFor(AminoAcid a);
};

};

#endif // CODON_H_
