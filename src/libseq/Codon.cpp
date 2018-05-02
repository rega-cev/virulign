#include "Codon.h"

namespace seq {

AminoAcid Codon::translate(const NTSequence::const_iterator triplet)
{
  const AminoAcid codonTable[4][4][4] = {
  { { AminoAcid::K /* AAA */,
      AminoAcid::N /* AAC */,
      AminoAcid::K /* AAG */,
      AminoAcid::N /* AAT */
    },
    { AminoAcid::T /* ACA */,
      AminoAcid::T /* ACC */,
      AminoAcid::T /* ACG */,
      AminoAcid::T /* ACT */
    },
    { AminoAcid::R /* AGA */,
      AminoAcid::S /* AGC */,
      AminoAcid::R /* AGG */,
      AminoAcid::S /* AGT */
    },
    { AminoAcid::I /* ATA */,
      AminoAcid::I /* ATC */,
      AminoAcid::M /* ATG */,
      AminoAcid::I /* ATT */
    }
  },
  { { AminoAcid::Q /* CAA */,
      AminoAcid::H /* CAC */,
      AminoAcid::Q /* CAG */,
      AminoAcid::H /* CAT */
    },
    { AminoAcid::P /* CCA */,
      AminoAcid::P /* CCC */,
      AminoAcid::P /* CCG */,
      AminoAcid::P /* CCT */
    },
    { AminoAcid::R /* CGA */,
      AminoAcid::R /* CGC */,
      AminoAcid::R /* CGG */,
      AminoAcid::R /* CGT */
    },
    { AminoAcid::L /* CTA */,
      AminoAcid::L /* CTC */,
      AminoAcid::L /* CTG */,
      AminoAcid::L /* CTT */
    }
  },
  { { AminoAcid::E /* GAA */,
      AminoAcid::D /* GAC */,
      AminoAcid::E /* GAG */,
      AminoAcid::D /* GAT */
    },
    { AminoAcid::A /* GCA */,
      AminoAcid::A /* GCC */,
      AminoAcid::A /* GCG */,
      AminoAcid::A /* GCT */
    },
    { AminoAcid::G /* GGA */,
      AminoAcid::G /* GGC */,
      AminoAcid::G /* GGG */,
      AminoAcid::G /* GGT */
    },
    { AminoAcid::V /* GTA */,
      AminoAcid::V /* GTC */,
      AminoAcid::V /* GTG */,
      AminoAcid::V /* GTT */
    }
  },
  { { AminoAcid::STP /* TAA */,
      AminoAcid::Y /* TAC */,
      AminoAcid::STP /* TAG */,
      AminoAcid::Y /* TAT */
    },
    { AminoAcid::S /* TCA */,
      AminoAcid::S /* TCC */,
      AminoAcid::S /* TCG */,
      AminoAcid::S /* TCT */
    },
    { AminoAcid::STP /* TGA */,
      AminoAcid::C /* TGC */,
      AminoAcid::W /* TGG */,
      AminoAcid::C /* TGT */
    },
    { AminoAcid::L /* TTA */,
      AminoAcid::F /* TTC */,
      AminoAcid::L /* TTG */,
      AminoAcid::F /* TTT */
    }
  } };

  if (*triplet == Nucleotide::GAP
      && (*(triplet + 1) == Nucleotide::GAP)
      && (*(triplet + 2) == Nucleotide::GAP))
    return AminoAcid::GAP;

  if (triplet->isAmbiguity()
      || (triplet + 1)->isAmbiguity()
      || (triplet + 2)->isAmbiguity())
    return AminoAcid::X;

  return
    codonTable[triplet->intRep()]
              [(triplet + 1)->intRep()]
              [(triplet + 2)->intRep()];
}

std::set<AminoAcid>
Codon::translateAll(const NTSequence::const_iterator triplet)
{
  std::set<AminoAcid> result;

  NTSequence s(triplet, triplet + 3);

  std::vector<NTSequence> possibilities;
  s.nonAmbiguousSequences(possibilities);

  for (unsigned i = 0; i < possibilities.size(); ++i)
    result.insert(translate(possibilities[i].begin()));

  return result;
}

namespace {
  void addTriplet(std::set<NTSequence>& result,
		  Nucleotide c1, Nucleotide c2, Nucleotide c3)
  {
    NTSequence triplet;
    triplet.push_back(c1);
    triplet.push_back(c2);
    triplet.push_back(c3);

    result.insert(triplet);
  }

}

std::set<NTSequence> Codon::codonsFor(AminoAcid a)
{
  std::set<NTSequence> result;

  switch (a.intRep()) {
  case AminoAcid::AA_A:
    addTriplet(result, Nucleotide::G, Nucleotide::C, Nucleotide::T);
    addTriplet(result, Nucleotide::G, Nucleotide::C, Nucleotide::C);
    addTriplet(result, Nucleotide::G, Nucleotide::C, Nucleotide::A);
    addTriplet(result, Nucleotide::G, Nucleotide::C, Nucleotide::G);
    break;
  case AminoAcid::AA_C:
    addTriplet(result, Nucleotide::T, Nucleotide::G, Nucleotide::T);
    addTriplet(result, Nucleotide::T, Nucleotide::G, Nucleotide::C);
    break;
  case AminoAcid::AA_D:
    addTriplet(result, Nucleotide::G, Nucleotide::A, Nucleotide::T);
    addTriplet(result, Nucleotide::G, Nucleotide::A, Nucleotide::C);
    break;
  case AminoAcid::AA_E:
    addTriplet(result, Nucleotide::G, Nucleotide::A, Nucleotide::A);
    addTriplet(result, Nucleotide::G, Nucleotide::A, Nucleotide::G);
    break;
  case AminoAcid::AA_F:
    addTriplet(result, Nucleotide::T, Nucleotide::T, Nucleotide::T);
    addTriplet(result, Nucleotide::T, Nucleotide::T, Nucleotide::C);
    break;
  case AminoAcid::AA_G:
    addTriplet(result, Nucleotide::G, Nucleotide::G, Nucleotide::T);
    addTriplet(result, Nucleotide::G, Nucleotide::G, Nucleotide::C);
    addTriplet(result, Nucleotide::G, Nucleotide::G, Nucleotide::A);
    addTriplet(result, Nucleotide::G, Nucleotide::G, Nucleotide::G);
    break;   
  case AminoAcid::AA_H:
    addTriplet(result, Nucleotide::C, Nucleotide::A, Nucleotide::T);
    addTriplet(result, Nucleotide::C, Nucleotide::A, Nucleotide::C);
    break;    
  case AminoAcid::AA_I:
    addTriplet(result, Nucleotide::A, Nucleotide::T, Nucleotide::T);
    addTriplet(result, Nucleotide::A, Nucleotide::T, Nucleotide::C);
    addTriplet(result, Nucleotide::A, Nucleotide::T, Nucleotide::A);
    break;    
  case AminoAcid::AA_K:
    addTriplet(result, Nucleotide::A, Nucleotide::A, Nucleotide::A);
    addTriplet(result, Nucleotide::A, Nucleotide::A, Nucleotide::G);
    break;    
  case AminoAcid::AA_L:
    addTriplet(result, Nucleotide::T, Nucleotide::T, Nucleotide::A);
    addTriplet(result, Nucleotide::T, Nucleotide::T, Nucleotide::G);
    addTriplet(result, Nucleotide::C, Nucleotide::T, Nucleotide::T);
    addTriplet(result, Nucleotide::C, Nucleotide::T, Nucleotide::C);
    addTriplet(result, Nucleotide::C, Nucleotide::T, Nucleotide::A);
    addTriplet(result, Nucleotide::C, Nucleotide::T, Nucleotide::G);
    break;    
  case AminoAcid::AA_M:
    addTriplet(result, Nucleotide::A, Nucleotide::T, Nucleotide::G);
    break;
  case AminoAcid::AA_N:
    addTriplet(result, Nucleotide::A, Nucleotide::A, Nucleotide::T);
    addTriplet(result, Nucleotide::A, Nucleotide::A, Nucleotide::C);
    break;    
  case AminoAcid::AA_P:
    addTriplet(result, Nucleotide::C, Nucleotide::C, Nucleotide::T);
    addTriplet(result, Nucleotide::C, Nucleotide::C, Nucleotide::C);
    addTriplet(result, Nucleotide::C, Nucleotide::C, Nucleotide::A);
    addTriplet(result, Nucleotide::C, Nucleotide::C, Nucleotide::G);
    break;
  case AminoAcid::AA_Q:
    addTriplet(result, Nucleotide::C, Nucleotide::A, Nucleotide::A);
    addTriplet(result, Nucleotide::C, Nucleotide::A, Nucleotide::G);
    break;
  case AminoAcid::AA_R:
    addTriplet(result, Nucleotide::C, Nucleotide::G, Nucleotide::T);
    addTriplet(result, Nucleotide::C, Nucleotide::G, Nucleotide::C);
    addTriplet(result, Nucleotide::C, Nucleotide::G, Nucleotide::A);
    addTriplet(result, Nucleotide::C, Nucleotide::G, Nucleotide::G);
    addTriplet(result, Nucleotide::A, Nucleotide::G, Nucleotide::A);
    addTriplet(result, Nucleotide::A, Nucleotide::G, Nucleotide::G);
    break;
  case AminoAcid::AA_S:
    addTriplet(result, Nucleotide::T, Nucleotide::C, Nucleotide::T);
    addTriplet(result, Nucleotide::T, Nucleotide::C, Nucleotide::C);
    addTriplet(result, Nucleotide::T, Nucleotide::C, Nucleotide::A);
    addTriplet(result, Nucleotide::T, Nucleotide::C, Nucleotide::G);
    addTriplet(result, Nucleotide::A, Nucleotide::G, Nucleotide::T);
    addTriplet(result, Nucleotide::A, Nucleotide::G, Nucleotide::C);
    break;
  case AminoAcid::AA_T:
    addTriplet(result, Nucleotide::A, Nucleotide::C, Nucleotide::T);
    addTriplet(result, Nucleotide::A, Nucleotide::C, Nucleotide::C);
    addTriplet(result, Nucleotide::A, Nucleotide::C, Nucleotide::A);
    addTriplet(result, Nucleotide::A, Nucleotide::C, Nucleotide::G);
    break;
  case AminoAcid::AA_V:
    addTriplet(result, Nucleotide::G, Nucleotide::T, Nucleotide::T);
    addTriplet(result, Nucleotide::G, Nucleotide::T, Nucleotide::C);
    addTriplet(result, Nucleotide::G, Nucleotide::T, Nucleotide::A);
    addTriplet(result, Nucleotide::G, Nucleotide::T, Nucleotide::G);
    break;
  case AminoAcid::AA_W:
    addTriplet(result, Nucleotide::T, Nucleotide::G, Nucleotide::G);
    break;
  case AminoAcid::AA_Y:
    addTriplet(result, Nucleotide::T, Nucleotide::A, Nucleotide::T);
    addTriplet(result, Nucleotide::T, Nucleotide::A, Nucleotide::C);
    break;
  case AminoAcid::AA_STP:
    addTriplet(result, Nucleotide::T, Nucleotide::A, Nucleotide::A);
    addTriplet(result, Nucleotide::T, Nucleotide::A, Nucleotide::G);
    addTriplet(result, Nucleotide::T, Nucleotide::G, Nucleotide::A);
    break;
  case AminoAcid::AA_GAP:
  case AminoAcid::AA_Z:
  case AminoAcid::AA_U:
  case AminoAcid::AA_B:
  case AminoAcid::AA_X:
  default:
    break;
  }

  return result;
}

};
