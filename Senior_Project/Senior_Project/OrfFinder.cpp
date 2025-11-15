#include "OrfFinder.h"
#include "DNASequence.h"
#include <iostream>

using namespace std;

//defining stop codons
const vector<string> OrfFinder::STOP_CODONS = { "TAA", "TAG", "TGA" };

vector<ORF> OrfFinder::findORFs(const std::string& sequence, TrieIndex& trie) {
    vector<ORF> allORFs;

    for (int frame = 1; frame <= 3; frame++) {
        string frameSeq = sequence.substr(frame - 1); //scanning 3 forward frames and starting from position 0, 1, or 2
        scanFrame(frameSeq, frame - 1, frame, trie, allORFs);
    }

    //the reverse frames, negative nums
    string revComp = DNASequence::reverseComplement(sequence);
    for (int frame = 1; frame <= 3; frame++) {
        string frameSeq = revComp.substr(frame - 1);
        scanFrame(frameSeq, frame - 1, -frame, trie, allORFs); 
    }

    return allORFs;
}

void OrfFinder::scanFrame(const std::string& frameSequence, int frameOffset,
    int frameNumber, TrieIndex& trie, std::vector<ORF>& results) {

    for (size_t pos = 0; pos + 2 < frameSequence.length(); pos += 3) {  //Moving through the frame in steps of 3 
        string codon = frameSequence.substr(pos, 3);

        if (codon == "ATG") { //checking if this is a start codon, then look for the next stop codon
            for (size_t stopPos = pos + 3; stopPos + 2 < frameSequence.length(); stopPos += 3) {
                string stopCodon = frameSequence.substr(stopPos, 3);

                bool isStop = false;
                for (const auto& stop : STOP_CODONS) {
                    if (stopCodon == stop) {
                        isStop = true;
                        break;
                    }
                }

                if (isStop) {
                    int orfLength = (stopPos - pos) + 3; // +3 to include the stop codon

                    //only keep ORFs longer than 100 nucleotides (about 33 codons) which is a common practice in biology since the smaller ORFs are considered insignificant 
                    if (orfLength >= 100) {
                        ORF orf;
                        orf.sequence = frameSequence.substr(pos, orfLength);
                        orf.length = orfLength;
                        orf.frame = frameNumber;
                        
                        //converting back to original sequence coordinates
                        if (frameNumber > 0) {
                            //forward frame
                            orf.start = frameOffset + pos;
                            orf.end = frameOffset + stopPos + 2;
                        }
                        else {
                            //reverse frame - wrong for now
                            orf.start = pos;  
                            orf.end = stopPos + 2;
                        }

                        results.push_back(orf);
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