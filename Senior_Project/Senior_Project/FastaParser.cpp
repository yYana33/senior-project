#include "FastaParser.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

//opening the file
unique_ptr<DNASequence> FastaParser::parseFromFile(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Could not open the file!" << filename << endl;
        return nullptr;
    }

    stringstream buffer;
    buffer << file.rdbuf();
    string fileContent = buffer.str();
    return parseFromString(fileContent);
}

unique_ptr<DNASequence> FastaParser::parseFromString(const std::string& data) {
    string header;
    string sequence;
    stringstream ss(data);
    string line;
    bool isFirstLine = true;
    bool hasHeader = false;
    bool hasSequence = false;

    while (getline(ss, line)) {
        if (line.empty()) continue;

        //checking if the line is a header line 
        if (line[0] == '>') {
            if (!isFirstLine) {
                cerr << "Error! Multiple sequences found in FASTA file; only single sequences allowed." << endl;
                return nullptr;
            }
            //extracting the header (removing the '>')
            header = line.substr(1);
            //carriage return remove if present 
            if (!header.empty() && header.back() == '\r') {
                header.pop_back();
            }
            hasHeader = true;
            isFirstLine = false;
        }
        else {
            //any line that doesn't start with `>` is part of the DNA sequence.
            string cleanLine = line;
            if (!cleanLine.empty() && cleanLine.back() == '\r') {
                cleanLine.pop_back();
            }

            for (char c : cleanLine) {
                char upperC = toupper(c);
                if (upperC != 'A' && upperC != 'T' && upperC != 'C' && upperC != 'G') { 
                    cerr << "Error! Invalid character '" << c << "' found in sequence." << endl;
                    cerr << "Only A, T, C, and G are allowed." << endl;
                    return nullptr;
                }
            }

            sequence += cleanLine;
            hasSequence = true;
        }
    } 

    if (!hasHeader) {
        cerr << "Error! No header line (starting with '>') found in FASTA file." << endl;
        return nullptr;
    }

    if (!hasSequence) {
        cerr << "Error! No sequence data found in FASTA file." << endl;
        return nullptr;
    }

    if (sequence.empty()) {
        cerr << "Error! Sequence is empty." << endl;
        return nullptr;
    }

    //Final validation
    auto dnaSeq = make_unique<DNASequence>(sequence, header);
    if (!dnaSeq->isValidDNA()) {
        cerr << "Error! Sequence validation failed." << endl;
        return nullptr;
    }

    return dnaSeq;
}
   
    /*
    if (!header.empty() && !sequence.empty()) {
        return make_unique<DNASequence>(sequence, header);
    }

    cerr << "Invalid file format!" << endl;
    return nullptr;
}
*/