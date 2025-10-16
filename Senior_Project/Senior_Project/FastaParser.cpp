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

    while (getline(ss, line)) {
        if (line.empty()) continue;

        //checking if the line is a header line 
        if (line[0] == '>') {
            if (!isFirstLine) {
           
                break;
            }
            //extracting the header (removing the '>')
            header = line.substr(1);
            //carriage return remove if present 
            if (!header.empty() && header.back() == '\r') {
                header.pop_back();
            }
            isFirstLine = false;
        }
        else {
            //any line that doesn’t start with `>` is part of the DNA sequence.
            string cleanLine = line;
            if (!cleanLine.empty() && cleanLine.back() == '\r') {
                cleanLine.pop_back();
            }
            sequence += cleanLine;
        }
    }
    //creating a DNASequence object
    if (!header.empty() && !sequence.empty()) {
        return make_unique<DNASequence>(sequence, header);
    }

    cerr << "Invalid file format!" << endl;
    return nullptr;
}