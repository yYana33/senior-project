#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "FastaParser.h"
#include "SearchAlgorithm.h"

using namespace std;

int main() {
    cout << "DNA Analyzer Visualizer" << endl;
    cout << "\n" << endl;
    cout << "real FASTA file test" << endl;

    //file 
    string fastaFilePath = "C:/Users/temp/Desktop/senior_project/Senior_Project/Senior_Project/data/sample.fasta";

    cout << "\n1. Testing with the FASTA file: " << fastaFilePath << endl;

    auto sequence = FastaParser::parseFromFile(fastaFilePath);
    if (sequence) {
        sequence->printSummary();

        //additional useful info for now 
        cout << "\nAdditional Stats:" << endl;
        cout << "\n" << endl;

        //counting the  bases + their frequencies
        int a_count = 0, t_count = 0, c_count = 0, g_count = 0;
        const string& seq = sequence->getSequence();

        for (char base : seq) {
            switch (toupper(base)) {
            case 'A': a_count++; break;
            case 'T': t_count++; break;
            case 'C': c_count++; break;
            case 'G': g_count++; break;
            }
        }

        cout << "A: " << a_count << " (" << (a_count * 100.0 / seq.length()) << "%)" << endl;
        cout << "T: " << t_count << " (" << (t_count * 100.0 / seq.length()) << "%)" << endl;
        cout << "C: " << c_count << " (" << (c_count * 100.0 / seq.length()) << "%)" << endl;
        cout << "G: " << g_count << " (" << (g_count * 100.0 / seq.length()) << "%)" << endl;

        double gc_content = (g_count + c_count) * 100.0 / seq.length();
        cout << "GC Content: " << gc_content << "%" << endl;

        cout << "\n2. Testing Boyer-Moore Search Algorithm:" << endl;
        cout << "\n" << endl;

        //testing with common codons
        vector<string> testPatterns = {
            "ATG",    //start codon
            "TAA",    //stop codon
            "TAG",    //stop codon  
            "TGA",    //stop codon
            "GGCC",   //common restriction site for the enzyme HaeIII
            "ATAT",   //simple repeat
            "CG"      //CpG islands
        };

        for (const auto& pattern : testPatterns) {
            if (pattern.length() <= seq.length()) {
                auto positions = SearchAlgorithm::boyerMooreSearch(seq, pattern);
                SearchAlgorithm::printSearchResults(positions, pattern);
            }
        }

        //letting user test their own pattern
        cout << "\n3. Custom Pattern Search:" << endl;
        cout << "\n" << endl;
        cout << "Enter a DNA pattern to search for (or 'quit' to exit): ";

        string userPattern;
        cin >> userPattern;

        if (userPattern != "quit") {
            //validating it's a DNA pattern
            bool validPattern = true;
            for (char c : userPattern) {
                if (c != 'A' && c != 'T' && c != 'C' && c != 'G' &&
                    c != 'a' && c != 't' && c != 'c' && c != 'g') {
                    validPattern = false;
                    break;
                }
            }

            if (validPattern && userPattern.length() <= seq.length()) {
                //uppercase for consistency
                string upperPattern = userPattern;
                for (char& c : upperPattern) {
                    c = toupper(c);
                }

                auto customPositions = SearchAlgorithm::boyerMooreSearch(seq, upperPattern);
                SearchAlgorithm::printSearchResults(customPositions, upperPattern);
            }
            else {
                cout << "Invalid pattern or pattern is too long/too short." << endl;
            }
        }

    }
    else {
        cout << "ERROR: Failed to parse the FASTA file!" << endl;
        return 1;
    }

    cout << "\nPress Enter to exit...";
    cin.ignore(); 
    cin.get();   

    return 0;
}

