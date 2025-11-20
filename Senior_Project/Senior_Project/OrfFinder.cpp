#include "OrfFinder.h"
#include "DNASequence.h"
#include "Gene.h"
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

    //separating the forward and reverse ORFs for clearer display
    vector<ORF> forwardORFs, reverseORFs;
    for (const auto& orf : orfs) {
        if (orf.frame > 0) { forwardORFs.push_back(orf); }
        else { reverseORFs.push_back(orf);}
    }

    if (!forwardORFs.empty()) {
        cout << "\nFORWARD strand ORFs:" << endl;
        for (size_t i = 0; i < forwardORFs.size(); i++) {
            const ORF& orf = forwardORFs[i];
            cout << "ORF " << (i + 1) << " (frame +" << orf.frame << "):" << endl;
            cout << "Position: " << orf.start << " - " << orf.end << endl;
            cout << "Length: " << orf.length << " nt, " << (orf.length / 3) << " aa" << endl;
            cout << "Start: " << orf.sequence.substr(0, 9) << "..." << endl;
            cout << "End: ..." << orf.sequence.substr(orf.sequence.length() - 9) << endl << endl;
        }
    }

    if (!reverseORFs.empty()) {
        cout << "REVERSE strand ORFs:" << endl;
        for (size_t i = 0; i < reverseORFs.size(); i++) {
            const ORF& orf = reverseORFs[i];
            cout << "ORF " << (i + 1) << " (frame " << orf.frame << "):" << endl;
            cout << "Position: " << orf.start << " - " << orf.end << endl;
            cout << "Length: " << orf.length << " nt, " << (orf.length / 3) << " aa" << endl;
            cout << "Start: " << orf.sequence.substr(0, 9) << "..." << endl;
            cout << "End: ..." << orf.sequence.substr(orf.sequence.length() - 9) << endl << endl;
        }
    }

    cout << "Summary: " << forwardORFs.size() << " forward, " << reverseORFs.size() << " reverse ORFs found." << endl;
}

void OrfFinder::findAndAddORFs(DNASequence& sequence, TrieIndex& trie) {
    vector<ORF> orfs = findORFs(sequence.getSequence(), trie);

    // converting each ORF to a Gene feature 
    for (const auto& orf : orfs) { 
        int gcCount = 0; //calculating the cg content for now, could be changed
        for (char c : orf.sequence) {
            if (toupper(c) == 'G' || toupper(c) == 'C') {
                gcCount++;
            }
        }

        double gcContent = orf.sequence.empty() ? 0.0 : (double)gcCount / orf.sequence.length();

        // Gene feature 
        auto gene = make_unique<Gene>(orf.start, orf.end, orf.frame, gcContent);
        sequence.addFeature(std::move(gene));
    }

    cout << "Added " << orfs.size() << " ORFs as features to sequence" << endl;
}

std::vector<ORF> OrfFinder::filterByLength(const std::vector<ORF>& orfs, int minLength) {
    std::vector<ORF> filtered;
    for (const auto& orf : orfs) {
        if (orf.length >= minLength) {
            filtered.push_back(orf);
        }
    }
    return filtered;
}

std::vector<ORF> OrfFinder::filterByFrame(const std::vector<ORF>& orfs, int frame) {
    std::vector<ORF> filtered;
    for (const auto& orf : orfs) {
        if (orf.frame == frame) {
            filtered.push_back(orf);
        }
    }
    return filtered;
}

std::vector<ORF> OrfFinder::removeOverlaps(const std::vector<ORF>& orfs) {
    if (orfs.empty()) return orfs;


    std::vector<ORF> sorted = orfs;
    //Sortiing by start position
    std::sort(sorted.begin(), sorted.end(), [](const ORF& a, const ORF& b) { return a.start < b.start;});

    std::vector<ORF> nonOverlapping;
    nonOverlapping.push_back(sorted[0]);

    for (size_t i = 1; i < sorted.size(); i++) {
        const ORF& current = sorted[i];
        const ORF& last = nonOverlapping.back();

        //checking for overlaps
        if (current.start > last.end) { nonOverlapping.push_back(current);}
        else {
            //If overlap is found - keep the longer ORF
            if (current.length > last.length) { nonOverlapping.back() = current; }
        }
    }

    return nonOverlapping;
}