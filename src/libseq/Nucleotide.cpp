#include <ctype.h>
#include <stdlib.h>
#include <set>
#include <stdexcept>

#include "ParseException.h"
#include "Nucleotide.h"

namespace {

#ifdef _WIN32

double drand48()
{
	return (double(rand()) / RAND_MAX);
}

#endif

int sampleUniform(int one, int two)
{
  double d = drand48();
  
  return (d < 0.5 ? one : two);
}

int sampleUniform(int one, int two, int three)
{
  double d = drand48() * 3.;

  return (d < 1. ? one : (d < 2. ? two : three));
}

int sampleUniform(int one, int two, int three, int four)
{
  double d = drand48() * 4.;

  return (d < 1. ? one : (d < 2. ? two : (d < 3. ? three : four)));
}
};

namespace seq {

const char Nucleotide::NT_CHAR[] = {'A', 'C', 'G', 'T',
				    'M', 'R', 'W', 'S',
				    'Y', 'K', 'V', 'H',
				    'D', 'B', 'N', '-' };

const Nucleotide Nucleotide::A(Nucleotide::NT_A);
const Nucleotide Nucleotide::C(Nucleotide::NT_C);
const Nucleotide Nucleotide::G(Nucleotide::NT_G);
const Nucleotide Nucleotide::T(Nucleotide::NT_T);
const Nucleotide Nucleotide::M(Nucleotide::NT_M);
const Nucleotide Nucleotide::R(Nucleotide::NT_R);
const Nucleotide Nucleotide::W(Nucleotide::NT_W);
const Nucleotide Nucleotide::S(Nucleotide::NT_S);
const Nucleotide Nucleotide::Y(Nucleotide::NT_Y);
const Nucleotide Nucleotide::K(Nucleotide::NT_K);
const Nucleotide Nucleotide::V(Nucleotide::NT_V);
const Nucleotide Nucleotide::H(Nucleotide::NT_H);
const Nucleotide Nucleotide::D(Nucleotide::NT_D);
const Nucleotide Nucleotide::B(Nucleotide::NT_B);
const Nucleotide Nucleotide::N(Nucleotide::NT_N);
const Nucleotide Nucleotide::GAP(Nucleotide::NT_GAP);

Nucleotide::Nucleotide()
  : rep_(NT_N)
{ }

void Nucleotide::sampleAmbiguity()
{
  switch (rep_) {
  case NT_A:
  case NT_C:
  case NT_G:
  case NT_T:
  case NT_GAP:
    break;
  case NT_M:
    rep_ = sampleUniform(NT_A, NT_C); break;
  case NT_R:
    rep_ = sampleUniform(NT_A, NT_G); break;  
  case NT_W:
    rep_ = sampleUniform(NT_A, NT_T); break;
  case NT_S:
    rep_ = sampleUniform(NT_C, NT_G); break;
  case NT_Y:
    rep_ = sampleUniform(NT_C, NT_T); break;
  case NT_K:
    rep_ = sampleUniform(NT_G, NT_T); break;
  case NT_V:
    rep_ = sampleUniform(NT_A, NT_C, NT_G); break;
  case NT_H:
    rep_ = sampleUniform(NT_A, NT_C, NT_T); break;
  case NT_D:
    rep_ = sampleUniform(NT_A, NT_G, NT_T); break;
  case NT_B:
    rep_ = sampleUniform(NT_C, NT_G, NT_T); break;
  case NT_N:
    rep_ = sampleUniform(NT_A, NT_C, NT_G, NT_T); break;
  default:
    std::cerr << rep_ << std::endl;
    assert(false);
  }
}

Nucleotide Nucleotide::reverseComplement() const
{
  switch (rep_) {
  case NT_A: return NT_T;
  case NT_C: return NT_G;
  case NT_G: return NT_C;
  case NT_T: return NT_A;
  case NT_GAP: return NT_GAP;
  case NT_M: return /* AC -> TG */ NT_K;
  case NT_R: return /* AG -> TC */ NT_Y;
  case NT_W: return /* AT -> TA */ NT_W;
  case NT_S: return /* CG -> GC */ NT_S;
  case NT_Y: return /* CT -> GA */ NT_R;
  case NT_K: return /* GT -> CA */ NT_M;
  case NT_V: return /* ACG -> TGC */ NT_B;
  case NT_H: return /* ACT -> TGA */ NT_D;
  case NT_D: return /* AGT -> TCA */ NT_H;
  case NT_B: return /* CGT -> GCA */ NT_V;
  case NT_N: return NT_N;
  default:
    std::cerr << rep_ << std::endl;
    assert(false);
  }
}

Nucleotide Nucleotide::singleNucleotide(std::set<Nucleotide>& nucleotides)
{
	std::set<Nucleotide>::iterator itgap = nucleotides.find(GAP);
	if(itgap != nucleotides.end())
		nucleotides.erase(itgap);

	if (nucleotides.size() == 1)
		return *nucleotides.begin();
	
	std::set<Nucleotide> all;
	for(std::set<Nucleotide>::iterator it = nucleotides.begin(); it != nucleotides.end(); ++it) {
		std::vector<Nucleotide> t;
		it->nonAmbiguousNucleotides(t);
		all.insert(t.begin(), t.end());
	}
	bool nta = all.find(A) != all.end();
	bool ntc = all.find(C) != all.end();
	bool ntg = all.find(G) != all.end();
	bool ntt = all.find(T) != all.end();

	if (nta && ntc && ntg && ntt)
		return N;
	if (nta && ntc && ntg)
		return V;
	if (nta && ntc && ntt)
		return H;
	if (nta && ntg && ntt)
		return D;
	if (ntc && ntg && ntt)
		return B;
 	if (nta && ntc)
		return M;
	if (ntg && ntt)
		return K;
	if (nta && ntt)
		return W;
	if (ntg && ntc)
		return S;
	if (ntc && ntt)
		return Y;
	if (nta && ntg)
		return R;		

	throw std::runtime_error
	  ("Internal error in Nucleotide::singleNucleotide()");
}

/**
 * Get all non ambiguous nucleotides represented by this nucleotide.
 */ 
void Nucleotide::nonAmbiguousNucleotides(std::vector<Nucleotide>& result) const
{
  switch (rep_) {
  case NT_A:
  case NT_C:
  case NT_G:
  case NT_T:
  case NT_GAP:
    result.push_back(*this);
    break;
  case NT_M:
    result.push_back(A);
    result.push_back(C);
    break;
  case NT_R:
    result.push_back(A);
    result.push_back(G);
    break;
  case NT_W:
    result.push_back(A);
    result.push_back(T);
    break;
  case NT_S:
    result.push_back(C);
    result.push_back(G);
    break;
  case NT_Y:
    result.push_back(C);
    result.push_back(T);
    break;
  case NT_K:
    result.push_back(G);
    result.push_back(T);
    break;
  case NT_V:
    result.push_back(A);
    result.push_back(C);
    result.push_back(G);
    break;
  case NT_H:
    result.push_back(A);
    result.push_back(C);
    result.push_back(T);
    break;
  case NT_D:
    result.push_back(A);
    result.push_back(G);
    result.push_back(T);
    break;
  case NT_B:
    result.push_back(C);
    result.push_back(G);
    result.push_back(T);
    break;
  case NT_N:
    result.push_back(A);
    result.push_back(C);
    result.push_back(G);
    result.push_back(T);
    break;
  default:
    std::cerr << rep_ << std::endl;
    assert(false);
  }
}


std::ostream& operator<< (std::ostream& s, const Nucleotide nt)
{
  return s << nt.toChar();
}

};
