#include "SearchAlgorithm.h"
#include "features/SearchMatch.h"
#include <iostream>
#include <cstring>
#include <QDebug>

using namespace std;

vector<int> SearchAlgorithm::boyerMooreSearch(const std::string& text, const std::string& pattern) {
    vector<int> positions;

    //edge cases
    if (pattern.empty() || text.empty() || pattern.length() > text.length()) {
        return positions;
    }

    //both text and pattern to uppercase
    string upperText = text;
    string upperPattern = pattern;

    transform(upperText.begin(), upperText.end(), upperText.begin(), ::toupper);
    transform(upperPattern.begin(), upperPattern.end(), upperPattern.begin(), ::toupper);

    int m = upperPattern.length();
    int n = upperText.length();

    //creating bad character table
    int badchar[256];
    buildBadCharTable(upperPattern, badchar);

    int s = 0; // s is shift of the pattern with respect to text
    while (s <= (n - m)) {
        int j = m - 1;

        //reducing index j of pattern while characters match
        while (j >= 0 && upperPattern[j] == upperText[s + j]) {
            j--;
        }

        //if the pattern is present at current shift
        if (j < 0) {
            positions.push_back(s);

            //shifting the pattern to find next occurrence
            s += (s + m < n) ? m - badchar[upperText[s + m]] : 1;
        } else {
            //shifting the pattern based on bad character rule
            s += max(1, j - badchar[upperText[s + j]]);
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
    sequence.removeFeaturesOfType("search_match");

    vector<int> positions = boyerMooreSearch(sequence.getSequence(), pattern);

    // creating SearchMatch features for each found position
    for (int pos : positions) {
        int endPos = pos + pattern.length() - 1;
        qDebug() << "Creating SearchMatch at position:" << pos << "to" << endPos;
        auto match = std::make_unique<SearchMatch>(pos, endPos, pattern);
        sequence.addFeature(std::move(match));
    }
}
