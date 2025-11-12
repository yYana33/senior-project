#pragma once
#include "TrieNode.h"
#include <string>
#include <vector>

class TrieIndex {
private:
    TrieNode* root;
    int kmerLength; //3 for codons

public:
    TrieIndex(int length = 3); 
    ~TrieIndex();

    void insertKmer(const std::string& kmer, int position);

    void buildIndex(const std::string& sequence);

    std::vector<int> findPositions(const std::string& kmer) const;//all positions of a specific k-mer

    bool contains(const std::string& kmer) const;

    int getKmerLength() const;

    std::vector<std::string> getAllKmers() const;//debugging

private:
    void collectKmers(TrieNode* node, std::string& current, std::vector<std::string>& results) const;
};