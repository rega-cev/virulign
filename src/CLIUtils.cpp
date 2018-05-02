// This may look like C code, but it's really -*- C++ -*-

#include <fstream> 
#include <string.h>
#include <stdexcept>

#include "ReferenceSequence.h" 
#include "CLIUtils.h" 

ReferenceSequence loadRefSeqFromFile(const char* refSeqFileName) {
  std::ifstream f_ref(refSeqFileName);
  if (!f_ref) {
    throw std::runtime_error(std::string("Could not open ") + refSeqFileName);
  }

  ReferenceSequence* ref;

  try {
    seq::NTSequence refNt;
    f_ref >> refNt;

    if (!f_ref) {
      throw std::runtime_error(std::string("RefSeq loading:: File ") + refSeqFileName + " does not contain a FASTA sequence ?");
    }

    ref = new ReferenceSequence(refNt);
    ref->addRegion(ReferenceSequence::Region(0, refNt.size()/3, "P"));

    f_ref >> refNt;

  if (f_ref) {
      throw std::runtime_error(std::string("RefSeq loading:: File ") + refSeqFileName + " contains multiple sequences ?");
    }
	return *ref;
  } catch (seq::ParseException& e) {
      throw std::runtime_error(std::string("RefSeq loading:: Fatal error: ") + e.message());
  }
}

bool equalsS(char* str1, char* str2) {
  return strcmp(str1, str2) == 0;
}

bool equalsString(std::string str1, std::string str2){
  return str1.compare(str2) == 0;
}
