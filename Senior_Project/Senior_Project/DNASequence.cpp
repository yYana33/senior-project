#include "DNASequence.h"
#include <algorithm>
#include <iostream>

using namespace std;

DNASequence::DNASequence(const std::string& seq, const std::string& hdr)
    : sequence(seq), header(hdr) {}

bool DNASequence::isValidDNA() const {
    return std::all_of(sequence.begin(), sequence.end(), [](char c) {
        return c == 'A' || c == 'T' || c == 'C' || c == 'G' || c == 'a' || c == 't' || c == 'c' || c == 'g';
        });
}

void DNASequence::printSummary() const {
    cout << "=== DNA Sequence Summary ===" << endl;
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