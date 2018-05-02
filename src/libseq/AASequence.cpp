#include "AASequence.h"
#include "Codon.h"

namespace seq {

AASequence::AASequence()
  : std::vector<AminoAcid>()
{ }

AASequence::AASequence(unsigned size)
  : std::vector<AminoAcid>(size)
{ }

AASequence::AASequence(const const_iterator first,
		       const const_iterator last)
  : std::vector<AminoAcid>(first, last)
{ }

AASequence::AASequence(const std::string name,
		       const std::string description,
		       const std::string aSeqString)
  : std::vector<AminoAcid>(aSeqString.length()),
    name_(name),
    description_(description)
{
  for (unsigned i = 0; i < aSeqString.length(); ++i) {
    (*this)[i] = AminoAcid(aSeqString[i]);
  }
}

std::string AASequence::asString() const
{
  std::string result(size(), '-');

  for (unsigned i = 0; i < size(); ++i) {
    result[i] = (*this)[i].toChar();
  }

  return result;
}

inline bool contains(const std::set<AminoAcid>& possibilities, const AminoAcid& aa)
{
  return possibilities.find(aa) != possibilities.end();
}

AASequence AASequence::translate(const NTSequence::const_iterator begin,
				 const NTSequence::const_iterator end)
{
  const int size = end - begin;
  assert(size % 3 == 0);

  AASequence result(size / 3);

  for (NTSequence::const_iterator i = begin; i < end; i += 3) {
    std::set<AminoAcid> possibilities = Codon::translateAll(i);
    if (possibilities.size() > 2) {
      result[(i - begin)/3] = AminoAcid::X;
    } else if (possibilities.size() == 2) {
      if (contains(possibilities,AminoAcid::D) && contains(possibilities,AminoAcid::N))
        result[(i - begin)/3] = AminoAcid::B;
      else if (contains(possibilities,AminoAcid::E) && contains(possibilities,AminoAcid::Q))
        result[(i - begin)/3] = AminoAcid::Z;
      else if (contains(possibilities,AminoAcid::L) && contains(possibilities,AminoAcid::I))
        result[(i - begin)/3] = AminoAcid::J;
      else
        result[(i - begin)/3] = AminoAcid::X;
    } else {
      result[(i - begin)/3] = *possibilities.begin();
    }
  }

  return result;
}

AASequence AASequence::translate(const NTSequence& ntSequence)
{
  return translate(ntSequence.begin(), ntSequence.end());
}

// defined in NTSequence.C:
extern void readFastaEntry(std::istream& i,
			   std::string& name,
			   std::string& description,
			   std::string& sequence);
extern void writeFastaEntry(std::ostream& o,
			    const std::string& name,
			    const std::string& description,
			    const std::string& sequence);

std::istream& operator>>(std::istream& i, AASequence& sequence)
{
  std::string name, description, seqString;

  readFastaEntry(i, name, description, seqString);
  sequence = AASequence(name, description, seqString);

  return i;
}

std::ostream& operator<<(std::ostream& o, const AASequence& sequence)
{
  writeFastaEntry(o, sequence.name(), sequence.description(),
		  sequence.asString());

  return o;
}

};
