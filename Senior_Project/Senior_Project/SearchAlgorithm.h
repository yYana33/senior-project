#pragma once
#include <string>
#include <vector>

class SearchAlgorithm {
public:
    //Boyer-Moore search 
    static std::vector<int> boyerMooreSearch(const std::string& text, const std::string& pattern);

    //displaying result 
    static void printSearchResults(const std::vector<int>& positions, const std::string& pattern);

private:

    static void buildBadCharTable(const std::string& pattern, int badchar[256]);
    static int max(int a, int b) { return (a > b) ? a : b; }
};

