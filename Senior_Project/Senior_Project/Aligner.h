#pragma once
#include <string>
#include <vector>
#include <utility> 

struct Alignment {
    int score;
    std::string sequence1;  
    std::string sequence2;  
    std::string matchLine;  
    int start1;             
    int end1;               
    int start2;             
    int end2;               
};

class Aligner {
private:
    int matchScore;
    int mismatchScore;
    int gapPenalty;


public:
    Aligner(int match = 2, int mismatch = -1, int gap = -2);//The default scoring parameters
    Alignment align(const std::string& seq1, const std::string& seq2);//performing the alignment
    static void printAlignment(const Alignment& alignment);
    static double calculateSimilarity(const Alignment& alignment);//a similarity percentage

private:
    int score(char a, char b) const;
    void traceback(const std::vector<std::vector<int>>& matrix, const std::string& seq1, const std::string& seq2, int maxI, int maxJ, Alignment& result);
};