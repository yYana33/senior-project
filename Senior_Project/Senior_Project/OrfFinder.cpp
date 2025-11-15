#include "OrfFinder.h"
#include "DNASequence.h"
#include <iostream>
#include <algorithm>

using namespace std;

//defining stop codons
const vector<string> OrfFinder::STOP_CODONS = { "TAA", "TAG", "TGA" };

vector<ORF> OrfFinder::findORFs(const std::string& sequence, TrieIndex& trie) {
    vector<ORF> allORFs;
    //3 forward frames
    for (int frame = 1; frame <= 3; frame++) {
        string frameSeq = sequence.substr(frame - 1); //starts from positions 0, 1, or 2
        scanFrame(frameSeq, sequence, frame - 1, frame, trie, allORFs);
    }

    //3 reverse frames
    string revComp = DNASequence::reverseComplement(sequence);
    for (int frame = 1; frame <= 3; frame++) {
        string frameSeq = revComp.substr(frame - 1);
        scanFrame(frameSeq, sequence, frame - 1, -frame, trie, allORFs);
    }

    return allORFs;
}

void OrfFinder::scanFrame(const std::string& frameSequence, const std::string& originalSequence,
    int frameOffset, int frameNumber, TrieIndex& trie, std::vector<ORF>& results) {
    
    for (size_t pos = 0; pos + 2 < frameSequence.length(); pos += 3) { string codon = frameSequence.substr(pos, 3); //moving through the frame in steps of 3 (codons)
        

        if (codon == "ATG") {
            for (size_t stopPos = pos + 3; stopPos + 2 < frameSequence.length(); stopPos += 3) {  //if this is a start codon look for the next stop codon
                string stopCodon = frameSequence.substr(stopPos, 3);

                bool isStop = false;
                for (const auto& stop : STOP_CODONS) {
                    if (stopCodon == stop) {
                        isStop = true;
                        break;
                    }
                }

                if (isStop) {
                    //found an ORF from pos to stopPos+2
                    int orfLength = (stopPos - pos) + 3;

                    if (orfLength >= 100) {
                        bool hasInternalStop = false;
                        for (size_t checkPos = pos + 3; checkPos < stopPos; checkPos += 3) {
                            string checkCodon = frameSequence.substr(checkPos, 3);
                            for (const auto& stop : STOP_CODONS) {
                                if (checkCodon == stop) {
                                    hasInternalStop = true;
                                    break;
                                }
                            }
                            if (hasInternalStop) break;
                        }

                        if (!hasInternalStop) {
                            ORF orf;
                            orf.sequence = frameSequence.substr(pos, orfLength);
                            orf.length = orfLength;
                            orf.frame = frameNumber;

                            //coordinate mapping
                            if (frameNumber > 0) {
                                orf.start = frameOffset + pos;
                                orf.end = frameOffset + stopPos + 2;
                            }
                            else {
                                orf.start = originalSequence.length() - (frameOffset + stopPos + 3);
                                orf.end = originalSequence.length() - (frameOffset + pos) - 1;
                                if (orf.start > orf.end) {
                                    std::swap(orf.start, orf.end); //swap if needed for now - will fix 
                                }
                            }

                            results.push_back(orf);
                        }
                    }
                    break;
                }
            }
        }
    }
}


void OrfFinder::printORFs(const std::vector<ORF>& orfs) {
    cout << "Found " << orfs.size() << " ORFs:" << endl;

    for (size_t i = 0; i < orfs.size(); i++) {
        const ORF& orf = orfs[i];
        cout << "ORF " << (i + 1) << ":" << endl;
        cout << "  Frame: " << orf.frame << endl;
        cout << "  Position: " << orf.start << " - " << orf.end << endl;
        cout << "  Length: " << orf.length << " nucleotides (" << (orf.length / 3) << " codons)" << endl;
        cout << "  Sequence: " << orf.sequence.substr(0, 30); //showing first 30 bases - test
        if (orf.sequence.length() > 30) {
            cout << "...";
        }
        cout << endl << endl;
    }
}