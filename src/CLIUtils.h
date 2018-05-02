// This may look like C code, but it's really -*- C++ -*-
#ifndef CLI_UTILS_H_ 
#define CLI_UTILS_H_ 

#include <string>

ReferenceSequence loadRefSeqFromFile(const char* refSeqFileName); 
bool equalsS(char* str1, char* str2); 
bool equalsString(std::string str1, std::string str2);

#endif // CLI_UTILS_H_ 
