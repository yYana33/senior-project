#include "TrieIndex.h"
#include <iostream>

using namespace std;

TrieIndex::TrieIndex(int length) : kmerLength(length) {
    root = new TrieNode();
}

TrieIndex::~TrieIndex() {
    delete root; //destructor 
}

void TrieIndex::insertKmer(const std::string& kmer, int position) {
    if (kmer.length() != kmerLength) {
        cerr << "Error: K-mer length must be " << kmerLength << ", got " << kmer.length() << endl;
        return;
    }

    TrieNode* current = root;

    for (char nucleotide : kmer) {//converting to upercase and traversing the trie, more nodes if needed
        char upperNuc = toupper(nucleotide);
        if (!current->hasChild(upperNuc)) {
            current->addChild(upperNuc);
        }

        current = current->getChild(upperNuc); //to the child node
    }

    current->addPosition(position);
}

void TrieIndex::buildIndex(const std::string& sequence) {
    int seqLength = sequence.length();

    for (int i = 0; i <= seqLength - kmerLength; i++) {
        string kmer = sequence.substr(i, kmerLength);
        insertKmer(kmer, i); 
    }

    cout << "Built Trie index with " << (seqLength - kmerLength + 1) << " k-mers" << endl;
}

std::vector<int> TrieIndex::findPositions(const std::string& kmer) const {
    vector<int> results;

    if (kmer.length() != kmerLength) {
        cerr << "Error: Search pattern length must be " << kmerLength << endl;
        return results;
    }

    TrieNode* current = root;

    for (char nucleotide : kmer) { //traversing
        char upperNuc = toupper(nucleotide);
        current = current->getChild(upperNuc);

        if (current == nullptr) {
            return results;
        }
    }

    return current->getPositions();
}

bool TrieIndex::contains(const std::string& kmer) const {
    return !findPositions(kmer).empty();
}

int TrieIndex::getKmerLength() const {
    return kmerLength;
}

std::vector<std::string> TrieIndex::getAllKmers() const {
    vector<string> results;
    string current;
    collectKmers(root, current, results);
    return results;
}

void TrieIndex::collectKmers(TrieNode* node, std::string& current, std::vector<std::string>& results) const {
    //if there is a k-mer of the correct length and it has positions
    if (current.length() == kmerLength && !node->getPositions().empty()) {
        results.push_back(current);
        return;
    }

    if (current.length() < kmerLength) {
        for (const auto& pair : node->getChildren()) {
            current.push_back(pair.first);
            collectKmers(pair.second, current, results);
            current.pop_back();
        }
    }
}