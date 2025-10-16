#pragma once
#include <string>
#include <iostream>

class DNASequence {
private:
    std::string sequence;
    std::string header;

public:
    DNASequence(const std::string& seq, const std::string& hdr);

    const std::string& getSequence() const { return sequence; }
    const std::string& getHeader() const { return header; }
    size_t length() const { return sequence.length(); }

    bool isValidDNA() const;
    void printSummary() const;
};