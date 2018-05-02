// This may look like C code, but it's really -*- C++ -*-
#ifndef RESULTS_EXPORTER_H_
#define RESULTS_EXPORTER_H_

#include <iostream>
#include <vector>

class Alignment;

enum ExportKind { Mutations, PairwiseAlignments, GlobalAlignment,
		  PositionTable, MutationTable };
enum ExportAlphabet { Nucleotides, AminoAcids };

class ResultsExporter
{
public:
  ResultsExporter(const std::vector<Alignment>& results, ExportKind kind,
		  ExportAlphabet alphabet, bool withInsertions = false);

  ExportKind     kind()     const { return kind_; }
  ExportAlphabet alphabet() const { return alphabet_; }

  void streamData(std::ostream& stream);
	void streamConsensusSequence(std::ostream& stream);

private:
  const std::vector<Alignment>& results_;
  const ExportKind      kind_;
  const ExportAlphabet  alphabet_;
  const bool            withInsertions_;

  void streamMutationsCsv(std::ostream& stream);
  void streamPairwiseAlignments(std::ostream& stream);
  void streamPositionTable(std::ostream& stream);
  void streamMutationTable(std::ostream& stream);

  void computeGlobalAlignment(seq::NTSequence& globalRef,
			      std::vector<seq::NTSequence>& globalAlignment);
  void streamGlobalAlignment(std::ostream& stream);
};

#endif // RESULTS_EXPORTER_H_
