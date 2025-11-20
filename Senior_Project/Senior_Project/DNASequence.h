#pragma once
#include <string>
#include <iostream>
#include <vector>
#include "Feature.h"

class DNASequence {
private:
    std::string sequence;
    std::string header;
    std::vector<std::unique_ptr<Feature>> features;

public:
    DNASequence(const std::string& seq, const std::string& hdr);

    const std::string& getSequence() const { return sequence; }
    const std::string& getHeader() const { return header; }
    size_t length() const { return sequence.length(); }

    const std::vector<std::unique_ptr<Feature>>& getFeatures() const { return features; }

  
    void addFeature(std::unique_ptr<Feature> feature);
    void clearFeatures();
    void removeFeaturesOfType(const std::string& type);

    bool isValidDNA() const;
    void printSummary() const;

    static std::string reverseComplement(const std::string& seq);
};



 

