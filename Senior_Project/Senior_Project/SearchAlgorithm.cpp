#include "SearchAlgorithm.h"
#include "SearchMatch.h"
#include <iostream>
#include <cstring>

using namespace std;

vector<int> SearchAlgorithm::boyerMooreSearch(const std::string& text, const std::string& pattern) {
    vector<int> positions;
    int m = pattern.length();
    int n = text.length();

    if (m == 0 || n == 0 || m > n) {
        return positions; //returns empty vector for invalid inputs
    }

    //creating the table
    int badchar[256];
    buildBadCharTable(pattern, badchar);

    int s = 0; //the shift of the pattern
    while (s <= (n - m)) {
        int j = m - 1;

        //reducing index j while characters match
        while (j >= 0 && pattern[j] == text[s + j]) {
            j--;
        }

        if (j < 0) { //if the pattern is present at current shift
            positions.push_back(s);

            //shifting the pattern to find next occurrence
            s += (s + m < n) ? m - badchar[text[s + m]] : 1;
        }
        else { //bad character rule
            s += max(1, j - badchar[text[s + j]]);
        }
    }
    return positions;
}

void SearchAlgorithm::buildBadCharTable(const std::string& pattern, int badchar[256]) {
    int m = pattern.length();

    //initializing all occurrences as -1
    for (int i = 0; i < 256; i++) {
        badchar[i] = -1;
    }
   
    for (int i = 0; i < m; i++) {  //the actual value of last occurrence of a character
        badchar[(int)pattern[i]] = i;
    }
}

void SearchAlgorithm::printSearchResults(const std::vector<int>& positions, const std::string& pattern) {
    if (positions.empty()) {
        cout << "Pattern '" << pattern << "' not found." << endl;
        return;
    }

    cout << "Pattern '" << pattern << "' found at positions: ";
    for (size_t i = 0; i < positions.size(); i++) {
        cout << positions[i];
        if (i < positions.size() - 1) {
            cout << ", ";
        }
    }
    cout << " (Total: " << positions.size() << " occurrences)" << endl;
}

void SearchAlgorithm::searchAndAddMatches(DNASequence& sequence, const std::string& pattern) {
    //clearing existing search matches 
    sequence.removeFeaturesOfType("search_match");

    vector<int> positions = boyerMooreSearch(sequence.getSequence(), pattern);

    //SearchMatch features for each found position
    for (int pos : positions) {
        int endPos = pos + pattern.length() - 1;
        auto match = make_unique<SearchMatch>(pos, endPos, pattern);
        sequence.addFeature(std::move(match));
    }

    cout << "Added " << positions.size() << " search matches for pattern '" << pattern << "'" << endl;
}