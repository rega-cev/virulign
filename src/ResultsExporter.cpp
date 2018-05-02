#include "Utils.h"

#include <set>
#include <algorithm>

#include <AASequence.h>
#include <Codon.h>

#include "Alignment.h"
#include "ResultsExporter.h"
#include "ReferenceSequence.h"

ResultsExporter::ResultsExporter(const std::vector<Alignment>& results,
				 ExportKind kind,
				 ExportAlphabet alphabet,
				 bool withInsertions)
  : results_(results),
    kind_(kind),
    alphabet_(alphabet),
    withInsertions_(withInsertions)
{ }

void ResultsExporter::streamData(std::ostream& stream)
{
  switch (kind_) {
  case Mutations:
    streamMutationsCsv(stream);
    break;
  case PairwiseAlignments:
    streamPairwiseAlignments(stream);
    break;
  case GlobalAlignment:
    streamGlobalAlignment(stream);
    break;
  case PositionTable:
    streamPositionTable(stream);
    break;
  case MutationTable:
    streamMutationTable(stream);
  }
}

namespace {
  seq::AASequence translate(seq::NTSequence seq)
  {
    seq::AASequence result = seq::AASequence::translate(seq);
    result.setName(seq.name());
    result.setDescription(seq.description());

    return result;
  }
}

void ResultsExporter::streamPairwiseAlignments(std::ostream& s)
{
  for (unsigned i = 0; i < results_.size(); ++i) {
    const Alignment& result = results_[i];

    if (alphabet_ == Nucleotides) {
      seq::NTSequence seq = result.ref;
      seq.setDescription(seq.description() + " aligned for "
			 + result.target.name());
      s << seq;
      s << result.target;
    } else {
      seq::AASequence seq = ::translate(result.ref);
      seq.setDescription(seq.description() + " aligned for "
			 + result.target.name());
      s << seq;
      s << ::translate(result.target);
    }
  }
}

namespace {

void alignToGlobalAlignment(seq::NTSequence& globalRef,
			    std::vector<seq::NTSequence >& globalTargets,
			    seq::NTSequence& ref,
			    seq::NTSequence& target,
			    bool withInsertions)
{
  /*
   * Copy gaps from globalRef to ref and vice-versa, adjusting the
   * respective targets.
   */
  for (unsigned i = 0; i < std::max(globalRef.size(), ref.size()); ++i) {
    if (i >= globalRef.size() || i >= ref.size() || (globalRef[i] != ref[i])) {
      if (i < globalRef.size() && globalRef[i] == seq::Nucleotide::GAP) {
	if (!withInsertions) {
	  std::cerr << globalRef
		    << ref
		    << target
		    << "i = " << i << std::endl;
	}
	assert(withInsertions);

	ref.insert(ref.begin() + i, seq::Nucleotide::GAP);
	target.insert(target.begin() + i, seq::Nucleotide::GAP);
      } else {
	assert(i > ref.size() || ref[i] == seq::Nucleotide::GAP);

	if (withInsertions) {
	  globalRef.insert(globalRef.begin() + i, seq::Nucleotide::GAP);
	  for (unsigned j = 0; j < globalTargets.size(); ++j) {
	    seq::NTSequence& t = globalTargets[j];
	    t.insert(t.begin() + i, seq::Nucleotide::GAP);
	  }
	} else {
	  ref.erase(ref.begin() + i);
	  target.erase(target.begin() + i);
	  --i;
	}
      }
    }
  }

  assert(ref == globalRef);
}

}

void ResultsExporter
::computeGlobalAlignment(seq::NTSequence& globalRef,
			 std::vector<seq::NTSequence>& globalAlignment)
{
  std::cerr << "Computing global alignment...";
  globalRef = results_[0].ref;

  if (!withInsertions_) {
    for (;;) {
      seq::NTSequence::iterator i
	= std::find(globalRef.begin(), globalRef.end(), seq::Nucleotide::GAP);
      if (i != globalRef.end())
	globalRef.erase(i);
      else
	break;
    }
  } else {
    /*
     * No insertions before first and after last nucleotide
     */
    while (globalRef[0] == seq::Nucleotide::GAP)
      globalRef.erase(globalRef.begin());
    while (globalRef[globalRef.size() - 1] == seq::Nucleotide::GAP)
      globalRef.erase(globalRef.begin() + globalRef.size() - 1);
  }

  for (unsigned j = 0; j < results_.size(); ++j) {
    if (results_[j].success) {
      seq::NTSequence ref = results_[j].ref;
      seq::NTSequence target = results_[j].target;

      /*
       * No insertions before first and after last nucleotide
       */
      while (ref[0] == seq::Nucleotide::GAP) {
	ref.erase(ref.begin());
	target.erase(target.begin());
      }
      while (ref[ref.size() - 1] == seq::Nucleotide::GAP) {
	ref.erase(ref.begin() + ref.size() - 1);
	target.erase(target.begin() + target.size() - 1);
      }

      alignToGlobalAlignment(globalRef, globalAlignment, ref, target,
			     withInsertions_);

      globalAlignment.push_back(target);
    }
  }

  std::cerr << " done." << std::endl;
}

void ResultsExporter::streamGlobalAlignment(std::ostream& s)
{
  if (results_.empty())
    return;

  seq::NTSequence globalRef;
  std::vector<seq::NTSequence> globalAlignment;
  computeGlobalAlignment(globalRef, globalAlignment);

  if (alphabet_ == Nucleotides) {
    for (unsigned i = 0; i < globalAlignment.size(); ++i)
      s << globalAlignment[i];
  } else {
    for (unsigned i = 0; i < globalAlignment.size(); ++i)
      s << ::translate(globalAlignment[i]);
  }
}

void ResultsExporter::streamConsensusSequence(std::ostream& s)
{
	if (results_.empty())
		return;

	seq::NTSequence globalRef;
	std::vector<seq::NTSequence> globalAlignment;
	computeGlobalAlignment(globalRef, globalAlignment);
	
	s << ">consensus" << std::endl;
	for (unsigned pos = 0; pos < globalRef.size(); ++pos) {
		std::set<seq::Nucleotide> all;
		for (unsigned i = 0; i < globalAlignment.size(); ++i) {
			all.insert(globalAlignment[i][pos]);
		}
		s << seq::Nucleotide::singleNucleotide(all);
	}
}

namespace {

int alignedAAPos(const seq::NTSequence& seq, int aapos)
{
  int j = 0;
  int pos = 0;

  // -XLFM--
  // aapos = 3 -> result = 4
  // ---XLFM
  // aapos = 1000 -> result = 6

  while ((pos < aapos)
	 && (((j+1)*3) < seq.size())) {
    if (seq[j*3] != seq::Nucleotide::GAP)
      ++pos;
    ++j;
  }

  return j;
}

}

void ResultsExporter::streamPositionTable(std::ostream& s)
{
  if (results_.empty())
    return;

  const ReferenceSequence& ref = results_[0].ref;

  seq::NTSequence globalRef;
  std::vector<seq::NTSequence> globalAlignment;
  computeGlobalAlignment(globalRef, globalAlignment);

  s << "seqid";

  for (unsigned r = 0; r < ref.regions().size(); ++r) {
    const ReferenceSequence::Region& region = ref.regions()[r];

    int first = alignedAAPos(globalRef, region.begin());
    int last = alignedAAPos(globalRef, region.end() - 1);

    int pos = 0;
    int insert = 0;

    for (int j = first; j <= last; ++j) {
      if (globalRef[j*3] != seq::Nucleotide::GAP) {
	++pos;
	if (alphabet_ == Nucleotides)
	  s << "," << region.prefix() << "_" << pos << "_1"
	    << "," << region.prefix() << "_" << pos << "_2"
	    << "," << region.prefix() << "_" << pos << "_3";
	else
	  s << "," << region.prefix() << "_" << pos;

	insert = 0;
      } else {
	++insert;

	if (alphabet_ == Nucleotides)
	  s << "," << region.prefix() << "_" << pos << "ins" << insert << "_1"
	    << "," << region.prefix() << "_" << pos << "ins" << insert << "_2"
	    << "," << region.prefix() << "_" << pos << "ins" << insert << "_3";
	else
	  s << "," << region.prefix() << "_" << pos << "ins" << insert;
      }
    }
  }
  s << std::endl;

  for (unsigned i = 0; i < globalAlignment.size(); ++i) {
    const seq::NTSequence& seq = globalAlignment[i];

    s << seq.name();

    for (unsigned r = 0; r < ref.regions().size(); ++r) {
      const ReferenceSequence::Region& region = ref.regions()[r];

      int first = alignedAAPos(globalRef, region.begin());
      int last = alignedAAPos(globalRef, region.end() - 1);

      int seqLast = last;
      while ((seqLast >= first)
	     && (seq[seqLast * 3] == seq::Nucleotide::GAP))
	--seqLast;

      bool beforeFirst = true;

      for (int j = first; j <= seqLast; ++j) {
	if (seq[j*3] == seq::Nucleotide::GAP && beforeFirst) {
	  if (alphabet_ == Nucleotides)
	    s << ",,,";
	  else
	    s << ",";
	} else {
	  beforeFirst = false;
	  if (alphabet_ == Nucleotides)
	    s << "," << seq[j*3]
	      << "," << seq[j*3 + 1]
	      << "," << seq[j*3 + 2];
	  else {
	    std::set<seq::AminoAcid>
	      aas = seq::Codon::translateAll(seq.begin() + j*3);

	    s << ",";
	    for (std::set<seq::AminoAcid>::const_iterator k = aas.begin();
		 k != aas.end(); ++k)
	      s << *k;
	  }
	}
      }

      for (int j = seqLast+1; j <= last; ++j) {
	if (alphabet_ == Nucleotides)
	  s << ",,,";
	else
	  s << ",";
      }
    }

    s << std::endl;
  }
}

void ResultsExporter::streamMutationTable(std::ostream& s)
{
  if (results_.empty())
    return;

  const ReferenceSequence& ref = results_[0].ref;

  seq::NTSequence globalRef;
  std::vector<seq::NTSequence> globalAlignment;
  computeGlobalAlignment(globalRef, globalAlignment);

  std::vector<std::set<seq::AminoAcid> > aminoAcids;
  for (unsigned j = 0; j < globalRef.size(); ++j)
    aminoAcids.push_back(std::set<seq::AminoAcid>());

  for (unsigned i = 0; i < globalAlignment.size(); ++i) {
    const seq::NTSequence& seq = globalAlignment[i];

    for (unsigned j = 0; j < seq.size(); j += 3) {
      std::set<seq::AminoAcid> 
	aas = seq::Codon::translateAll(seq.begin() + j);

      for (std::set<seq::AminoAcid>::const_iterator k = aas.begin();
	   k != aas.end(); ++k)
	if (*k != seq::AminoAcid::GAP)
	  aminoAcids[j/3].insert(*k);
    }
  }

  s << "seqid";

  for (unsigned r = 0; r < ref.regions().size(); ++r) {
    const ReferenceSequence::Region& region = ref.regions()[r];

    int first = alignedAAPos(globalRef, region.begin());
    int last = alignedAAPos(globalRef, region.end() - 1);

    int pos = 0;
    int insert = 0;

    for (int j = first; j <= last; ++j) {
      std::string varName;
      if (globalRef[j*3] != seq::Nucleotide::GAP) {
	++pos;
	varName = region.prefix() + "_" + 
	  to_string(pos);
	insert = 0;
      } else {
	++insert;
	varName = region.prefix() + "_" + to_string(pos)
	  + "ins" + to_string(insert);
      }
      
      for (std::set<seq::AminoAcid>::const_iterator k = aminoAcids[j].begin();
	   k != aminoAcids[j].end(); ++k)
	s << "," << varName << *k;
    }
  }

  s << std::endl;

  for (unsigned i = 0; i < globalAlignment.size(); ++i) {
    const seq::NTSequence& seq = globalAlignment[i];

    s << seq.name();

    for (unsigned r = 0; r < ref.regions().size(); ++r) {
      const ReferenceSequence::Region& region = ref.regions()[r];

      int first = alignedAAPos(globalRef, region.begin());
      int last = alignedAAPos(globalRef, region.end() - 1);

      int seqLast = last;
      while ((seqLast >= first)
	     && (seq[seqLast * 3] == seq::Nucleotide::GAP))
	--seqLast;

      bool beforeFirst = true;

      for (int j = first; j <= last; ++j) {
	if (seq[j*3] != seq::Nucleotide::GAP)
	  beforeFirst = false;

	std::set<seq::AminoAcid>
	  aas = seq::Codon::translateAll(seq.begin() + j*3);

	for (std::set<seq::AminoAcid>::const_iterator 
	       k = aminoAcids[j].begin();
	     k != aminoAcids[j].end(); ++k)
	  if (aas.find(*k) != aas.end())
	    s << ",y";
	  else
	    if (beforeFirst || j > seqLast)
	      s << ",";
	    else
	      s << ",n";
      }
    }

    s << std::endl;
  }  
}

void ResultsExporter::streamMutationsCsv(std::ostream& s)
{
  if (results_.empty())
    return;

  const ReferenceSequence& ref = results_[0].ref;

  s << "seqid,status,score,frameshifts";

  for (unsigned i = 0; i < ref.regions().size(); ++i) {
    std::string prefix;
    if (ref.regions().size() > 1)
      prefix = " " + ref.regions()[i].prefix();

    s << ",begin" << prefix
      << ",end" << prefix
      << ",mutations" << prefix;
  }
  s << std::endl;

  for (unsigned i = 0; i < results_.size(); ++i) {
    const Alignment& result = results_[i];

    s << result.target.name();

    if (result.success)
      s << ",Success";
    else if (result.tooShort)
      s << ",FailTooShort";
    else if (result.failure)
      s << ",Failure";
    else
      s << ",InternalError";

    if (result.success) {
      s << "," << result.score
	<< "," << result.correctedFrameshifts;

      for (unsigned i = 0; i < result.ref.regions().size(); ++i) {
	const ReferenceSequence::Region& region = result.ref.regions()[i];
	int begin = region.targetBegin;
	int end = region.targetEnd;

	if (begin < end)
	  s << "," << begin - region.begin() + 1
	    << "," << end - region.begin() + 1;
	else
	  s << ",,";

	s << "," << result.mutations(region);
      }
    } else {
      s << ",,";
      for (unsigned i = 0; i < result.ref.regions().size(); ++i)
	s << ",";
    }

    s << std::endl;
  }
}
