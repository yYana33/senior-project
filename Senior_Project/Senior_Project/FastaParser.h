#pragma once
#include <string>
#include <memory>
#include "DNASequence.h"

class FastaParser {
public:
    static std::unique_ptr<DNASequence> parseFromFile(const std::string& filename);
    static std::unique_ptr<DNASequence> parseFromString(const std::string& data);
};