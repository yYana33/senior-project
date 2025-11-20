#include "DNASequence.h"
#include <algorithm>
#include <iostream>

using namespace std;

DNASequence::DNASequence(const std::string& seq, const std::string& hdr) : sequence(seq), header(hdr) {}

bool DNASequence::isValidDNA() const {
    return std::all_of(sequence.begin(), sequence.end(), [](char c) {
        return c == 'A' || c == 'T' || c == 'C' || c == 'G' || c == 'a' || c == 't' || c == 'c' || c == 'g';
    });
}

void DNASequence::printSummary() const {
    cout << "DNA sequence summary" << endl;
    cout << "Header: " << header << endl;
    cout << "Length: " << length() << " bases" << endl;
    cout << "Valid DNA: " << (isValidDNA() ? "Yes" : "No") << endl;

    /*
    string preview = sequence.substr(0, 50);
    if (sequence.length() > 50) {
        preview += "...";
    }
    cout << "Preview: " << preview << endl;
    cout << "\n" << endl;
    */
}

    std::string DNASequence::reverseComplement(const std::string & seq) {
        string result;
        result.reserve(seq.length());

        for (int i = seq.length() - 1; i >= 0; i--) {
            char base = seq[i];
            switch (toupper(base)) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            default: result += 'N';//invalid characters
            }
        }
        return result;
    }

void DNASequence::addFeature(std::unique_ptr<Feature> feature) {
    features.push_back(std::move(feature));
}

void DNASequence::clearFeatures() {
    features.clear();
}

void DNASequence::removeFeaturesOfType(const std::string& type) {
    features.erase(std::remove_if(features.begin(), features.end(), [&](const std::unique_ptr<Feature>& feature) {
        return feature->getType() == type; }),
        features.end());
}