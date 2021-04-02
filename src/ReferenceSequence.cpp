#include "ReferenceSequence.h"

#include "mxml-utils/MXMLUtils.h"
#include "mxml/mxml.h"
#include <fstream>
#include <stdexcept>

#include "Utils.h"

ReferenceSequence::ReferenceSequence(const seq::NTSequence& seq)
  : seq::NTSequence(seq)
{

}

ReferenceSequence parseOrfReference(mxml_node_t* node)
{
  std::string orf;
  std::string refSeq;
  attributeValue(node, "name", orf);
  attributeValue(node, "referenceSequence", refSeq);
	  
  std::vector<mxml_node_t *> proteins = childElements(node, "protein");
	  
  ReferenceSequence seq 
    (seq::NTSequence(orf, orf, refSeq));
  for (int k = 0; k < proteins.size(); k++) {
    std::string protein;
    std::string start;
    std::string end;
    
    attributeValue(proteins[k], "abbreviation", protein);
    attributeValue(proteins[k], "startPosition", start);
    attributeValue(proteins[k], "stopPosition", end);
	    
    if (protein == "") 
        throw std::runtime_error("protein abbreviation is invalid");
    if (start == "") 
        throw std::runtime_error("protein start is invalid");
    if (end == "") 
        throw std::runtime_error("protein end is invalid");

    int startPos = (atoi(start.c_str()) - 1) / 3;
    int endPos = (atoi(end.c_str()) - 1) / 3;
	    
    ReferenceSequence::Region region(startPos, 
		  endPos, 
		  protein);
    seq.addRegion(region);
  }

  return seq;
}

ReferenceSequence ReferenceSequence::parseOrfReferenceFile(const std::string& fileName)
{
  FILE *fp = fopen(fileName.c_str(), "r");
  if (fp) {
    mxml_node_t *top = mxmlNewElement(MXML_NO_PARENT, "top");
  
    mxml_node_t *first = mxmlLoadFile(top, fp, MXML_NO_CALLBACK);
    
    if (first) {
      mxml_node_t *root = singleChildElement(top, "orf");
      return parseOrfReference(root);
    }
  }
    
  throw std::runtime_error("Error parsing ORF reference file");
}

std::map<std::string, std::vector<ReferenceSequence> >
ReferenceSequence::parseProteinReferences(std::string genomesXmlFile) 
{
  std::map<std::string, std::vector<ReferenceSequence> > genomesMap;

  FILE *fp = fopen(genomesXmlFile.c_str(), "r");
  if (fp) {
    mxml_node_t *top = mxmlNewElement(MXML_NO_PARENT, "top");

    mxml_node_t *first = mxmlLoadFile(top, fp, MXML_NO_CALLBACK);

    if (first) {
      mxml_node_t *root = singleChildElement(top, "genomes");

      std::vector<mxml_node_t *> genomes
        = childElements(root, "genome");

      for (int i = 0; i < genomes.size(); i++) {
	std::string organism;
	std::vector<ReferenceSequence> refs;
	attributeValue(genomes[i], "organismName", organism);
    if (organism == "") 
        throw std::runtime_error("organism name is invalid");

	std::vector<mxml_node_t *> orfs
	  = childElements(genomes[i], "openReadingFrame");
	for (int j = 0; j < orfs.size(); j++) {
	  ReferenceSequence seq = parseOrfReference(orfs[j]);
	  refs.push_back(seq);
	}
	genomesMap[organism] = refs;
      }
    }
  }

  return genomesMap;
}

ReferenceSequence::ReferenceSequence() {
}
