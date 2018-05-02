// This may look like C code, but it's really -*- C++ -*-
#ifndef PARSE_EXCEPTION_H_
#define PARSE_EXCEPTION_H_

#include <string>

namespace seq {

/**
 * Exception thrown when an error was encountered while parsing the
 * string representation of an nucleotide, nucleotide sequence, amino
 * acid, amino acid sequence, or a FASTA file.
 *
 * \sa Nucleotide::Nucleotide(char), AminoAcid::AminoAcid(char),
 * NTSequence::NTSequence(const std::string, const std::string, const
 * std::string, bool), AASequence::AASequence(const std::string, const
 * std::string, const std::string), operator>> (std::istream&,
 * NTSequence&), operator>> (std::istream&, AASequence&)
 */
class ParseException
{
public:
  ParseException(const std::string name,
		 const std::string message, bool recovered)
    : name_(name), message_(message), recovered_(recovered) { }

  /**
   * The sequence name.
   */
  std::string name() const { return name_; }

  /**
   * The message describing the error.
   */
  std::string message() const { return message_; }

  /**
   * Whether the parser attempted to recover and you could try parsing
   * the next sequence.
   */
  bool recovered() const { return recovered_; }

private:
  std::string name_;
  std::string message_;
  bool recovered_;
};

};

#endif // PARSE_EXCEPTION_H_
