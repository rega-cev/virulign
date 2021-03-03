#include <limits>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <omp.h>

#include <NeedlemanWunsh.h>

#include "ReferenceSequence.h"
#include "Alignment.h"
#include "ResultsExporter.h"
#include "CLIUtils.h"
#include "Utils.h"

ReferenceSequence loadRefSeq(const std::string& fn) {
  if (ends_with(fn, ".fasta")) {
    return loadRefSeqFromFile(fn.c_str());
  } else if (ends_with(fn, ".xml")) {
    return ReferenceSequence::parseOrfReferenceFile(fn);
  }
  throw std::runtime_error("Unsupported reference sequence format");
}

int main(int argc, char **argv) {
  unsigned int i;
	
  int obligatoryParams = 2;
  if(argc < obligatoryParams+1) {
    std::cerr << "Usage: virulign [reference.fasta orf-description.xml] sequences.fasta" << std::endl 
	      << "Optional parameters (first option will be the default):" << std::endl
	      << "  --exportKind [Mutations PairwiseAlignments GlobalAlignment PositionTable MutationTable]" << std::endl  
	      << "  --exportAlphabet [AminoAcids Nucleotides]" << std::endl
	      << "  --exportWithInsertions [yes no]" << std::endl
	      << "  --exportReferenceSequence [no yes]" << std::endl
	      << "  --gapExtensionPenalty doubleValue=>3.3" << std::endl
	      << "  --gapOpenPenalty doubleValue=>10.0" << std::endl
	      << "  --maxFrameShifts intValue=>3" << std::endl
              << "  --progress [no yes]" << std::endl
              << "  --nt-debug directory" << std::endl
	      << "Output: The alignment will be printed to standard out and any progress or error messages will be printed to the standard error. This output can be redirected to files, e.g.:" << std::endl
              << "   virulign ref.xml sequence.fasta > alignment.mutations 2> alignment.err" << std::endl;
    exit(0);
  }
	
  int amountOfParameters = argc - obligatoryParams - 1;
  if (amountOfParameters%2 == 1) {
    std::cerr << "Please provide parameters as: --parameterName parameterValue" << std::endl;	
    exit(0);
  } 

  std::string refSeqFileName = argv[1];
  if (!ends_with(refSeqFileName, ".fasta") && !ends_with(refSeqFileName, ".xml")) {
    std::cerr << 
      "Unknown reference sequence: "
      "expected a FASTA file or an XML file that describes the ORF" << std::endl;
    exit(1);
  }
  ReferenceSequence refSeq = loadRefSeq(refSeqFileName);

  std::ifstream f_seqs(argv[2]);
  std::vector<seq::NTSequence> targets; 

  try {
    while (f_seqs) {
      seq::NTSequence s;

      f_seqs >> s;

      if (f_seqs) {
	targets.push_back(s);
      }
    }
  } catch (seq::ParseException& e) {
    std::cerr << "Fatal error: " << e.message() << std::endl;
    exit(1);
  }

  ExportKind exportKind = Mutations;
  ExportAlphabet exportAlphabet = AminoAcids;
  bool exportWithInsertions = true;

  double gapExtensionPenalty = 3.3;
  double gapOpenPenalty = 10.0;
  int maxFrameShifts = 3;

  bool progress = false;

  std::string ntDebugDir;
	
  char* parameterName;
  char* parameterValue;
  for(i = obligatoryParams+1; i < amountOfParameters+obligatoryParams; i=i+2) {
    parameterName = argv[i];
    parameterValue = argv[i+1];
    if(equalsString(parameterName,"--exportKind")) {
      if(equalsString(parameterValue, "Mutations")) {
	exportKind = Mutations;
      } else if(equalsString(parameterValue, "PairwiseAlignments")) {
	exportKind = PairwiseAlignments;
      } else if(equalsString(parameterValue, "GlobalAlignment")) {
	exportKind = GlobalAlignment;
      } else if(equalsString(parameterValue, "PositionTable")) {
	exportKind = PositionTable;
      } else if(equalsString(parameterValue, "MutationTable")) {
	exportKind = MutationTable;
      } else {
	std::cerr << "Unkown value " << parameterValue << " for parameter : " << parameterName << std::endl; 
	exit(0);
      }	
    } else if(equalsString(parameterName,"--exportAlphabet")) {
      if(equalsString(parameterValue, "AminoAcids")) {
	exportAlphabet = AminoAcids;
      } else if(equalsString(parameterValue, "Nucleotides")) {
	exportAlphabet = Nucleotides;
      } else {
	std::cerr << "Unkown value " << parameterValue << " for parameter : " << parameterName << std::endl; 
	exit(0);
      }
    } else if(equalsString(parameterName,"--exportReferenceSequence")) {
      if (equalsString(parameterValue,"yes")) {
	seq::NTSequence refNtSeq = refSeq;
	targets.insert(targets.begin(), refNtSeq);
      }
    } else if(equalsString(parameterName,"--exportWithInsertions")) {
      if(equalsString(parameterValue,"yes")) {
	exportWithInsertions = true;
      } else if(equalsString(parameterValue,"no")) {
	exportWithInsertions = false;
      } else {
	std::cerr << "Unkown value " << parameterValue << " for parameter : " << parameterName << std::endl; 
	exit(0);
      } 
    } else if(equalsString(parameterName,"--gapExtensionPenalty")) {
      try {
	gapExtensionPenalty = lexical_cast<double>(parameterValue);
      } catch (std::bad_cast& e) {
	std::cerr << "Unkown value " << parameterValue << " for parameter : " << parameterName << std::endl; 
	exit(0);
      }
    } else if(equalsString(parameterName,"--gapOpenPenalty")) {
      try {
	gapOpenPenalty = lexical_cast<double>(parameterValue);
      } catch (std::bad_cast& e) {
        std::cerr << "Unkown value " << parameterValue << " for parameter : " << parameterName << std::endl;
        exit(0);
      }
    } else if(equalsString(parameterName,"--maxFrameShifts")) {
      try {
        maxFrameShifts = lexical_cast<int>(parameterValue);
      } catch (std::bad_cast& e) {
        std::cerr << "Unkown value " << parameterValue << " for parameter : " << parameterName << std::endl;
        exit(0);
      }
    } else if(equalsString(parameterName,"--progress")) {
      if(equalsString(parameterValue,"yes")) {
	progress = true;
      } else if(equalsString(parameterValue,"no")) {
	progress = false;
      } else {
	std::cerr << "Unkown value " << parameterValue << " for parameter : " << parameterName << std::endl; 
	exit(0);
      } 
    } else if(equalsString(parameterName,"--nt-debug")) {
      ntDebugDir = parameterValue;  
    } else {
      std::cerr << "Unkown parameter name: " << parameterName << std::endl; 
      exit(0);
    }
  }
	
  std::vector<Alignment> results;
 
  seq::NeedlemanWunsh algorithm(-gapOpenPenalty, -gapExtensionPenalty);

  if (!ntDebugDir.empty()) {
	seq::NTSequence r = refSeq;
    for (i = 0; i < targets.size(); ++i) {
      seq::NTSequence t = targets[i];
      double ntScore = algorithm.align(r, t);
      if(ntScore > 200) {
        std::string dbg = ntDebugDir + std::string("/") + t.name() + ".fasta";  
        std::ofstream ofs(dbg.c_str());
        ofs << r;
        ofs << t;
      }
    }
  }

  long int start = current_time_ms();
  
  #pragma omp parallel for
  for (i = 0; i < targets.size(); ++i) {
    std::cerr << "Align target " << i
            << " (" << targets[i].name() << ")" << std::endl;
    results.push_back(Alignment::compute(refSeq, targets[i], &algorithm, maxFrameShifts));
    if (progress) {
      long int end = current_time_ms();
      long int elapsed = end - start;
      double time_per_seq = (double)elapsed / (i + 1);
      double estimated_time_left = time_per_seq * (targets.size() - (i + 1));

      std::cerr << "Progress: " << (i + 1) << "/" << targets.size() << " sequences aligned (" << std::fixed << std::setprecision(2) << (i + 1) / (double)targets.size() * 100 << "%), Estimated time left " <<  format_time(estimated_time_left) << std::endl;
    }
  }

  ResultsExporter exporter(results, exportKind, exportAlphabet, exportWithInsertions);

  exporter.streamData(std::cout);
}
