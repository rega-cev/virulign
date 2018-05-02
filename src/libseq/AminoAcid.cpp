#include "AminoAcid.h"
#include "ParseException.h"

namespace seq {

const char AminoAcid::AA_CHAR[]
= { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
    'W', 'Y', '*', '-', 'Z', 'U', 'B', 'X', 'J' };

const char * const AminoAcid::AA_TLA[]
= { "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys",
    "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val",
    "Trp", "Tyr", "STP", "GAP", "Glu/Gln", "Sec", "Asp/Asn", "Any",
    "Leu/Ile" };

const AminoAcid AminoAcid::A(AminoAcid::AA_A);
const AminoAcid AminoAcid::C(AminoAcid::AA_C);
const AminoAcid AminoAcid::D(AminoAcid::AA_D);
const AminoAcid AminoAcid::E(AminoAcid::AA_E);
const AminoAcid AminoAcid::F(AminoAcid::AA_F);
const AminoAcid AminoAcid::G(AminoAcid::AA_G);
const AminoAcid AminoAcid::H(AminoAcid::AA_H);
const AminoAcid AminoAcid::I(AminoAcid::AA_I);
const AminoAcid AminoAcid::K(AminoAcid::AA_K);
const AminoAcid AminoAcid::L(AminoAcid::AA_L);
const AminoAcid AminoAcid::M(AminoAcid::AA_M);
const AminoAcid AminoAcid::N(AminoAcid::AA_N);
const AminoAcid AminoAcid::P(AminoAcid::AA_P);
const AminoAcid AminoAcid::Q(AminoAcid::AA_Q);
const AminoAcid AminoAcid::R(AminoAcid::AA_R);
const AminoAcid AminoAcid::S(AminoAcid::AA_S);
const AminoAcid AminoAcid::T(AminoAcid::AA_T);
const AminoAcid AminoAcid::V(AminoAcid::AA_V);
const AminoAcid AminoAcid::W(AminoAcid::AA_W);
const AminoAcid AminoAcid::Y(AminoAcid::AA_Y);
const AminoAcid AminoAcid::STP(AminoAcid::AA_STP);
const AminoAcid AminoAcid::GAP(AminoAcid::AA_GAP);
const AminoAcid AminoAcid::Z(AminoAcid::AA_Z);
const AminoAcid AminoAcid::U(AminoAcid::AA_U);
const AminoAcid AminoAcid::B(AminoAcid::AA_B);
const AminoAcid AminoAcid::X(AminoAcid::AA_X);
const AminoAcid AminoAcid::J(AminoAcid::AA_J);

AminoAcid::AminoAcid()
  : rep_(AA_Z)
{ }

AminoAcid::AminoAcid(char c)
{
  switch (toupper(c)) {
  case 'A': rep_ = AA_A; break;
  case 'C': rep_ = AA_C; break;
  case 'D': rep_ = AA_D; break;
  case 'E': rep_ = AA_E; break;
  case 'F': rep_ = AA_F; break;
  case 'G': rep_ = AA_G; break;
  case 'H': rep_ = AA_H; break;
  case 'I': rep_ = AA_I; break;
  case 'K': rep_ = AA_K; break;
  case 'L': rep_ = AA_L; break;
  case 'M': rep_ = AA_M; break;
  case 'N': rep_ = AA_N; break;
  case 'P': rep_ = AA_P; break;
  case 'Q': rep_ = AA_Q; break;
  case 'R': rep_ = AA_R; break;
  case 'S': rep_ = AA_S; break;
  case 'T': rep_ = AA_T; break;
  case 'V': rep_ = AA_V; break;
  case 'W': rep_ = AA_W; break;
  case 'Y': rep_ = AA_Y; break;
  case '*': rep_ = AA_STP; break;
  case '-': rep_ = AA_GAP; break;
  case 'Z': rep_ = AA_Z; break;
  case 'U': rep_ = AA_U; break;
  case 'B': rep_ = AA_B; break;
  case 'X': rep_ = AA_X; break;
  case 'J': rep_ = AA_J; break;
  default:
    throw ParseException
      (std::string(),
       std::string("Invalid amino acid character: '") + c + "'", false);
  }
}

std::string AminoAcid::tla() const
{
  return AA_TLA[rep_];
}

std::ostream& operator<< (std::ostream& s, const AminoAcid aa)
{
  return s << aa.toChar();
}

};
