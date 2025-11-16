#pragma once
#include <string>
#include <vector>
#include "TrieIndex.h"

//a found ORF
struct ORF {
    int start;          
    int end;            
    int length;         
    int frame; //+1, +2, +3 (forward) or -1, -2, -3 (reverse)
    std::string sequence; 
};

class OrfFinder {
private:
    static const std::vector<std::string> STOP_CODONS;

public:
    static std::vector<ORF> findORFs(const std::string& sequence, TrieIndex& trie);

    static void printORFs(const std::vector<ORF>& orfs);

    static std::vector<ORF> filterByLength(const std::vector<ORF>& orfs, int minLength);

    static std::vector<ORF> filterByFrame(const std::vector<ORF>& orfs, int frame);

    static std::vector<ORF> removeOverlaps(const std::vector<ORF>& orfs);

private:
    static void scanFrame(const std::string& frameSequence, const std::string& originalSequence, int frameOffset, int frameNumber, TrieIndex& trie, std::vector<ORF>& results);
};