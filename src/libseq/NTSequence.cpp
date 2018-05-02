#include <fstream>
#include <algorithm>

#include "NTSequence.h"
#include "ParseException.h"

namespace seq {

NTSequence::NTSequence()
  : std::vector<Nucleotide>()
{ }

NTSequence::NTSequence(unsigned size)
  : std::vector<Nucleotide>(size)
{ }

NTSequence::NTSequence(const std::string name, const std::string description,
		       const std::string aSeqString,
		       bool sampleAmbiguities)
  : std::vector<Nucleotide>(aSeqString.length()),
    name_(name),
    description_(description)
{
  for (unsigned i = 0; i < aSeqString.length(); ++i) {
    try {
      Nucleotide nt(aSeqString[i]);
      if (sampleAmbiguities)
	nt.sampleAmbiguity();

      (*this)[i] = nt;
    } catch (ParseException& e) {
      throw ParseException(name, e.message(), e.recovered());
    }
  }
}

NTSequence::NTSequence(const const_iterator first,
		       const const_iterator last)
  : std::vector<Nucleotide>(first, last)
{ }

void NTSequence::sampleAmbiguities()
{
  for (unsigned i = 0; i < size(); ++i) {
    (*this)[i].sampleAmbiguity();
  }
}

NTSequence NTSequence::reverseComplement() const
{
  NTSequence result(size());
  result.name_ = name_;
  result.description_ = description_;

  for (unsigned i = 0; i < size(); ++i)
    result[size() - i - 1] = (*this)[i].reverseComplement();

  return result;
}

void NTSequence::nonAmbiguousSequences(std::vector<NTSequence>& result) const
{
  iterateNonAmbiguous(NTSequence(), result);
}

void NTSequence::iterateNonAmbiguous(const NTSequence& head,
				     std::vector<NTSequence>& result) const
{
  /*
   * find the next ambigous codon (if any)
   */
  NTSequence s = head;
  unsigned i = head.size();
  for (; i < size(); ++i)
    if ((*this)[i].isAmbiguity())
      break;
    else
      s.push_back((*this)[i]);

  if (i == size())
    result.push_back(s);
  else {
    std::vector<Nucleotide> ambiguities;
    (*this)[i].nonAmbiguousNucleotides(ambiguities);
    for (unsigned i = 0; i < ambiguities.size(); ++i) {
      s.push_back(ambiguities[i]);
      iterateNonAmbiguous(s, result);
      s.pop_back();
    }
  }
}

std::string NTSequence::asString() const
{
  std::string result(size(), '-');

  for (unsigned i = 0; i < size(); ++i) {
    result[i] = (*this)[i].toChar();
  }

  return result;
}

/// \cond

void readFastaEntry(std::istream& i,
		    std::string& name,
		    std::string& description,
		    std::string& sequence)
{
    char ch;
    char c[512];

    i.getline(c, 511);
    if (i) {
      if (c[0] != '>') {
	throw ParseException(std::string(),
			     std::string("FASTA file expected '>', got: '")
			     + c[0] + "'", false);
      }

      std::string nameDesc = c + 1;
      std::string::size_type spacepos = nameDesc.find(" ");
      name = nameDesc.substr(0, spacepos);
      description = (spacepos == std::string::npos
		     ? ""
		     : nameDesc.substr(spacepos));

      for (ch = i.get(); (ch != EOF) && (ch != '>'); ch = i.get()) {
	if ((ch != '\n') && (ch != '\r') && (ch != ' ')) {
	  if (((ch >= 'a') && (ch <= 'z'))
	      || ((ch >= 'A') && (ch <= 'Z'))
	      || (ch == '-') || (ch == '*')) {
	    sequence += ch;
	  } else {
	    char failedCh = ch;
	    /*
	     * Wind further to the next possible sequence.
	     */
	    for (ch = i.get(); (ch != EOF) && (ch != '>'); ch = i.get())
	      ;

	    if (ch == '>')
	      i.putback(ch);

	    throw ParseException
	      (name, std::string("Illegal character in FASTA: '")
	       + (char)failedCh + "'", true);
	  }
	}

	if (i.peek() == EOF)
	  break;
      }

      if (ch == '>')
	i.putback(ch);
    }
}

void writeFastaEntry(std::ostream& o,
		     const std::string& name,
		     const std::string& description,
		     const std::string& sequence)
{
  o << ">" << name << " " << description << std::endl;
  if (sequence.size() == 0)
    o << std::endl;
  else {
    for (unsigned i = 0; i <= (sequence.size() - 1) / 60; ++i) {
      int s = i * 60;
      o << sequence.substr(s, 60) << std::endl;
    }
  }
}

void writeStockholm(std::ostream& o, const std::vector<NTSequence>& sequences, int length, int labelsize, int seqsize, int pos)
{
  if(labelsize < 1 && seqsize < 1){
    for(std::vector<NTSequence>::const_iterator i = sequences.begin();
        i < sequences.end(); ++i){
      labelsize = std::max(labelsize, (int)i->name().length());
      seqsize = std::max(seqsize, (int)i->size());
    }

    o << "# STOCKHOLM 1.0" << std::endl;
  }

  int epos = pos+length - (labelsize + 1);
  for(std::vector<NTSequence>::const_iterator i = sequences.begin();
      i < sequences.end(); ++i){
    o << i->name();
    for(int j = 0; j < labelsize - i->name().length() + 1; ++j)
      o << ' ';
    
    int n = std::min(epos, (int)i->size());
    for(int spos=pos; spos<n; ++spos)
      o << (*i)[spos];
  
    o << std::endl;
  }

  if(epos >= seqsize){
    o << "//";
  }
  else{
    writeStockholm(o, sequences, length, labelsize, seqsize, epos);
  }
}

/// \endcond

std::istream& operator>>(std::istream& i, NTSequence& sequence)
{
  std::string name, description, seqString;

  readFastaEntry(i, name, description, seqString);
  sequence = NTSequence(name, description, seqString);

  return i;
}

std::ostream& operator<<(std::ostream& o, const NTSequence& sequence)
{
  writeFastaEntry(o, sequence.name(), sequence.description(),
		  sequence.asString());
  return o;
}

};
