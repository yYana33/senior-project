#include "Aligner.h"
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

Aligner::Aligner(int match, int mismatch, int gap) : matchScore(match), mismatchScore(mismatch), gapPenalty(gap) {}


Alignment Aligner::align(const std::string& seq1, const std::string& seq2) {
    int m = seq1.length();
    int n = seq2.length();

    //Scoring matrix
    vector<vector<int>> matrix(m + 1, vector<int>(n + 1, 0));
    vector<vector<char>> direction(m + 1, vector<char>(n + 1, ' '));

    int maxScore = 0;
    int maxI = 0, maxJ = 0;




    //filling in the matrix
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int match = matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1]);
            int deleteGap = matrix[i - 1][j] + gapPenalty;
            int insertGap = matrix[i][j - 1] + gapPenalty;

            matrix[i][j] = max(0, max(match, max(deleteGap, insertGap)));

           
            //tracking direction for traceback
            if (matrix[i][j] == 0) {
                direction[i][j] = ' '; //stop
            }
            else if (matrix[i][j] == match) {
                direction[i][j] = 'D'; //diagonal (match/mismatch)
            }
            else if (matrix[i][j] == deleteGap) {
                direction[i][j] = 'U'; //up (gap in second sequence)
            }
            else {
                direction[i][j] = 'L'; //left (gap in first sequence)
            }

            //tracking the maximum score
            if (matrix[i][j] > maxScore) {
                maxScore = matrix[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    //traceback from the maximum score
    Alignment result;
    result.score = maxScore;
    traceback(matrix, seq1, seq2, maxI, maxJ, result);

    return result;
}

int Aligner::score(char a, char b) const {
    return (toupper(a) == toupper(b)) ? matchScore : mismatchScore;
}

void Aligner::traceback(const std::vector<std::vector<int>>& matrix,
    const std::string& seq1, const std::string& seq2,
    int maxI, int maxJ, Alignment& result) {
    string aligned1, aligned2, matchLine;
    int i = maxI, j = maxJ;

    //tracing back until there is a score 0
    while (i > 0 && j > 0 && matrix[i][j] > 0) {
        char a = seq1[i - 1];
        char b = seq2[j - 1];

        //determining the move direction by checking scores
        int current = matrix[i][j];
        int diagonal = matrix[i - 1][j - 1];
        int up = matrix[i - 1][j];
        int left = matrix[i][j - 1];

        if (current == diagonal + score(a, b)) {
            //diagonal move (match/mismatch)
            aligned1 = a + aligned1;
            aligned2 = b + aligned2;
            matchLine = (toupper(a) == toupper(b) ? '|' : ' ') + matchLine;
            i--;
            j--;
        }
        else if (current == up + gapPenalty) {
            //up move (gap in second sequence)
            aligned1 = a + aligned1;
            aligned2 = '-' + aligned2;
            matchLine = ' ' + matchLine;
            i--;
        }
        else if (current == left + gapPenalty) {
            //left move (gap in first sequence)
            aligned1 = '-' + aligned1;
            aligned2 = b + aligned2;
            matchLine = ' ' + matchLine;
            j--;
        }
        else {
            break; //should not happen
        }
    }

    result.sequence1 = aligned1;
    result.sequence2 = aligned2;
    result.matchLine = matchLine;
    result.start1 = i;      
    result.end1 = maxI - 1; 
    result.start2 = j;      
    result.end2 = maxJ - 1; 
}

void Aligner::printAlignment(const Alignment& alignment) {
    cout << "Alignment results:" << endl;
    cout << "Alignment score: " << alignment.score << endl;
    cout << "Similarity: " << fixed << setprecision(1) << calculateSimilarity(alignment) << "%" << endl;
    cout << "Sequence 1: " << alignment.start1 << " - " << alignment.end1 << endl;
    cout << "Sequence 2: " << alignment.start2 << " - " << alignment.end2 << endl;
    cout << endl;

    //displaying in blocks of 60 characters for readability
    int blockSize = 60;
    int totalLength = alignment.sequence1.length();

    for (int start = 0; start < totalLength; start += blockSize) {
        int end = min(start + blockSize, totalLength);

        cout << "Seq1: " << alignment.sequence1.substr(start, end - start) << endl;
        cout << "      " << alignment.matchLine.substr(start, end - start) << endl;
        cout << "Seq2: " << alignment.sequence2.substr(start, end - start) << endl;
        cout << endl;
    }

    //statistics
    int matches = count(alignment.matchLine.begin(), alignment.matchLine.end(), '|');
    int gaps = count(alignment.sequence1.begin(), alignment.sequence1.end(), '-') + count(alignment.sequence2.begin(), alignment.sequence2.end(), '-');

    cout << "Alignment Statistics:" << endl;
    cout << "Length: " << alignment.sequence1.length() << " bases" << endl;
    cout << "Matches: " << matches << endl;
    cout << "Gaps: " << gaps << endl;
    cout << "Identity: " << fixed << setprecision(1) << (matches * 100.0 / alignment.sequence1.length()) << "%" << endl;
}

double Aligner::calculateSimilarity(const Alignment& alignment) {
    if (alignment.sequence1.empty()) return 0.0;

    int matches = count(alignment.matchLine.begin(), alignment.matchLine.end(), '|');
    return (matches * 100.0) / alignment.sequence1.length();
}