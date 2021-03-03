// This may look like C code, but it's really -*- C++ -*-
#ifndef REFERENCE_SEQUENCE_H_
#define REFERENCE_SEQUENCE_H_

#include <NTSequence.h>

#include <map>
#include <vector>

class ReferenceSequence : public seq::NTSequence
{
public:
  class Region {
  public:

    Region(int begin, int end, std::string prefix)
      : begin_(begin),
	end_(end),
	prefix_(prefix)
    { }

    int         begin()   const { return begin_; }   // AA position [0 -- N[
    int         end()     const { return end_; }     // AA position [0 -- N[
    std::string prefix()  const { return prefix_; }
    
    // aligned positions of begin, end
    int         alignedBegin, alignedEnd; // AA position [0 -- N[
    // reference position of first/last non-gap in target within region
    int         targetBegin, targetEnd;   // AA position [0 -- N[

  private:
    int                       begin_, end_;
    std::string               prefix_;

    friend class ReferenceSequence;
  };

  ReferenceSequence(const seq::NTSequence& seq);
  
  const std::vector<Region>& regions() const { return regions_; }
  std::vector<Region>&       regions() { return regions_; }
  void addRegion(Region r) { regions_.push_back(r); }

  static std::map<std::string, std::vector<ReferenceSequence> > 
  parseProteinReferences(std::string genomesXmlFile);
  static ReferenceSequence 
  parseOrfReferenceFile(const std::string& fileName);
  ReferenceSequence();

private:
  std::vector<Region> regions_;
};

#endif // REFERENCE_SEQUENCE_H_
