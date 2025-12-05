#include "Controller.h"
#include "./core/Feature.h"
#include "./algorithms/SearchAlgorithm.h"
#include "SequenceVisualizer.h"
#include "./features/SearchMatch.h"
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <QDebug>
#include <QFileInfo>
#include <QFile>

using namespace std;

Controller::Controller(QObject* parent) : QObject(parent) {
}

bool Controller::loadSequence(const QString& filePath) {

    auto newSequence = FastaParser::parseFromFile(filePath.toStdString());

    if (newSequence) {
        currentSequence = std::move(newSequence);
        qDebug() << "Successfully loaded sequence! Length:" << currentSequence->length();

        calculateStatistics();
        updateTrieIndex();
        runAutoSearch();

        return true;
    } else {
        return false;
    }
}

double Controller::getGCPercentage() const {
    if (!currentSequence) return 0.0;

    const std::string& seq = currentSequence->getSequence();
    int gc = 0;
    for (char c : seq) {
        char b = std::toupper(c);
        if (b == 'G' || b == 'C') gc++;
    }
    return (seq.empty() ? 0.0 : (100.0 * gc / seq.size()));
}

void Controller::calculateStatistics() {
    if (!currentSequence) return;

    const std::string& sequence = currentSequence->getSequence();
    baseCounts.clear();

    std::vector<std::pair<std::string, int>> counts = {
        {"A", 0}, {"T", 0}, {"C", 0}, {"G", 0}
    };

    for (char c : sequence) {
        char upperC = toupper(c);
        if (upperC == 'A') counts[0].second++;
        else if (upperC == 'T') counts[1].second++;
        else if (upperC == 'C') counts[2].second++;
        else if (upperC == 'G') counts[3].second++;
    }

    baseCounts = counts;
    qDebug() << "Statistics calculated for sequence of length" << sequence.length();
}

QString Controller::getStatisticsInfo() const {
    if (!currentSequence) return "No sequence loaded";

    std::stringstream ss;
    ss << "Sequence Statistics:\n";

    int totalBases = currentSequence->length();

    for (const auto& baseCount : baseCounts) {
        double percentage = (baseCount.second * 100.0) / totalBases;
        ss << baseCount.first << ": " << baseCount.second << " (" << std::fixed << std::setprecision(2) << percentage << "%)\n";
    }

    int gcCount = 0;
    for (const auto& baseCount : baseCounts) {
        if (baseCount.first == "G" || baseCount.first == "C") {
            gcCount += baseCount.second;
        }
    }
    double gcPercentage = (gcCount * 100.0) / totalBases;
    ss << "GC content: " << std::fixed << std::setprecision(2) << gcPercentage << "%\n";
    ss << "Total length: " << totalBases << " bases\n";

    return QString::fromStdString(ss.str());
}

void Controller::runAutoSearch() {
    if (!currentSequence) return;

    std::vector<std::string> autoPatterns = {
        "ATG", "TAA", "TAG", "TGA", "GGCC", "ATAT", "CG"
    };

    std::stringstream ss;
    ss << "Automatic search for common biological motifs:\n";

    for (const auto& pattern : autoPatterns) {
        std::vector<int> positions = SearchAlgorithm::boyerMooreSearch(
            currentSequence->getSequence(), pattern);

        ss << "Pattern '" << pattern << "' found at positions: ";
        if (positions.empty()) {
            ss << "None";
        } else {
            //showing only the first 10 positionsto avoid overwhelming display
            size_t maxToShow = std::min(positions.size(), size_t(10));
            for (size_t i = 0; i < maxToShow; i++) {
                ss << positions[i];
                if (i < maxToShow - 1) {
                    ss << ", ";
                }
            }
            if (positions.size() > maxToShow) {
                ss << "... (and " << (positions.size() - maxToShow) << " more)";
            }
        }
        ss << " (Total: " << positions.size() << " occurrences)\n";
    }

    autoSearchResults = QString::fromStdString(ss.str());
}

bool Controller::loadSecondSequence(const QString& filePath) {

    auto newSequence = FastaParser::parseFromFile(filePath.toStdString());

    if (newSequence) {
        secondSequence = std::move(newSequence);
        return true;
    } else {
        return false;
    }
}

int Controller::getSearchMatchCount() const {
    if (!currentSequence) return 0;

    int count = 0;
    const auto& features = currentSequence->getFeatures();
    for (const auto& feature : features) {
        if (feature && feature->getType() == "search_match") {
            count++;
        }
    }
    return count;
}

int Controller::getGeneCount() const {
    if (!currentSequence) return 0;

    int count = 0;
    const auto& features = currentSequence->getFeatures();
    for (const auto& feature : features) {
        if (feature->getType() == "gene") {
            count++;
        }
    }
    return count;
}

void Controller::searchPattern(const QString& pattern) {
    if (!currentSequence) {
        return;
    }

    //uppercase for consistency
    std::string patternStr = pattern.toUpper().toStdString();
    currentPattern = pattern.toUpper();
    qDebug() << "Uppercase pattern:" << QString::fromStdString(patternStr);

    //Validating pattern char
    for (char c : patternStr) {
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
            emit errorOccurred(QString("Invalid DNA character '%1' in search pattern.\n"
            "Only A, T, C, and G are allowed.").arg(c));
            return;
        }
    }

    //Validating pattern length
    if (patternStr.length() < 2) {
        emit errorOccurred("Search pattern must be at least 2 characters long!");
        return;
    }

    if (patternStr.length() > 20) {
        emit errorOccurred("Search pattern cannot exceed 20 characters!");
        return;
    }

    SearchAlgorithm::searchAndAddMatches(*currentSequence, patternStr);

    emit searchCompleted(pattern, currentSequence->getFeatures().size());
}

void Controller::findORFs() {
    if (!currentSequence) {
        return;
    }

    OrfFinder::findAndAddORFs(*currentSequence, trieIndex);

    const auto& features = currentSequence->getFeatures();
    for (const auto& feature : features) {
        qDebug() << "Feature:" << QString::fromStdString(feature->getDisplayName()) << "at" << feature->getStart() << "-" << feature->getEnd();
    }
}

void Controller::alignSequences() {
    if (!currentSequence || !secondSequence) {
        return;
    }

    lastAlignment = aligner.align(
        currentSequence->getSequence(),
        secondSequence->getSequence());
    hasAlignmentResult = true;

    //aligned strings to QString for visualization
    QString aligned1 = QString::fromStdString(lastAlignment.sequence1);
    QString aligned2 = QString::fromStdString(lastAlignment.sequence2);

    alignmentBlocks = computeAlignmentBlocks(aligned1, aligned2);
}

std::vector<SequenceVisualizer::AlignBlock>
Controller::computeAlignmentBlocks(const QString& seq1Aligned, const QString& seq2Aligned)
{
    std::vector<SequenceVisualizer::AlignBlock> blocks;

    int n = seq1Aligned.size();
    if (n == 0) return blocks;

    int blockStart = 0;
    int matches = 0;
    int length = 0;

    for (int i = 0; i <= n; ++i)
    {
        bool end = (i == n);

        bool valid = !end && seq1Aligned[i] != '-' && seq2Aligned[i] != '-';
        bool match = valid && (seq1Aligned[i] == seq2Aligned[i]);

        if (!end)
        {
            if (valid)
            {
                length++;
                if (match) matches++;
            }
        }

        if (end || !valid)
        {
            if (length > 0)
            {
                double sim = (double)matches / length;

                blocks.push_back({
                    blockStart,
                    i - 1,
                    sim
                });
            }

            blockStart = i + 1;
            matches = 0;
            length = 0;
        }
    }

    return blocks;
}

std::vector<Controller::ORFRegion> Controller::getORFs() const {
    std::vector<ORFRegion> out;

    if (!currentSequence) return out;

    const auto& features = currentSequence->getFeatures();

    for (const auto& f : features) {
        if (f->getType() == "gene") {

            const Gene* g = dynamic_cast<const Gene*>(f.get());
            if (!g) continue;

            int frame = g->getReadingFrame();
            int strand = (frame > 0 ? +1 : -1);

            out.push_back({
                f->getStart(),
                f->getEnd(),
                strand
            });
        }
    }

    return out;
}


void Controller::clearResults() {
    if (currentSequence) {
        currentSequence->clearFeatures();
    }
}

void Controller::clearSearchResults() {
    if (currentSequence) {
        currentSequence->removeFeaturesOfType("search_match");
    }
}

void Controller::clearORFResults() {
    if (currentSequence) {
        currentSequence->removeFeaturesOfType("gene");
    }
}

void Controller::clearAlignmentResults() {
    hasAlignmentResult = false;
    lastAlignment = Alignment();
}

void Controller::clearAllResults() {
    clearSearchResults();
    clearORFResults();
    clearAlignmentResults();
}

bool Controller::hasSearchResults() const {
    if (!currentSequence) return false;
    const auto& features = currentSequence->getFeatures();
    for (const auto& feature : features) {
        if (feature->getType() == "search_match") {
            return true;
        }
    }
    return false;
}

bool Controller::hasORFResults() const {
    if (!currentSequence) return false;
    const auto& features = currentSequence->getFeatures();
    for (const auto& feature : features) {
        if (feature->getType() == "gene") {
            return true;
        }
    }
    return false;
}

bool Controller::hasAlignmentResults() const {
    return hasAlignmentResult;
}

QString Controller::getSequenceInfo() const {
    if (!currentSequence) return "No sequence loaded";

    stringstream ss;
    ss << "Sequence: " << currentSequence->getHeader() << "\n";
    ss << "Length: " << currentSequence->length() << " bases\n";
    ss << "Valid DNA: " << (currentSequence->isValidDNA() ? "Yes" : "No");

    return QString::fromStdString(ss.str());
}

QString Controller::getFeaturesInfo() const {
    if (!currentSequence) return "No sequence loaded";

    const auto& features = currentSequence->getFeatures();
    if (features.empty()) return "No results found";

    stringstream ss;
    ss << "Found " << features.size() << " features:\n\n";

    for (const auto& feature : features) {
        ss << "• " << feature->getDisplayName() << "\n";
        ss << "  " << feature->getDetails() << "\n\n";
    }

    return QString::fromStdString(ss.str());
}

QString Controller::getAlignmentInfo() const {
    if (!currentSequence || !secondSequence) {
        return "Please load two sequences for an alignment.";
    }

    if (!hasAlignmentResult) {
        return "Two sequences loaded. Use the Align button to compare them.";
    }

    std::stringstream ss;

    ss << "SEQUENCE ALIGNMENT RESULT\n\n";
    ss << "Alignment Score: " << lastAlignment.score << "\n";
    ss << "Similarity: " << Aligner::calculateSimilarity(lastAlignment) << "%\n";
    ss << "Aligned Length: " << lastAlignment.sequence1.length() << " bases\n";
    ss << "Sequence 1 Region: " << lastAlignment.start1 << " - " << lastAlignment.end1 << "\n";
    ss << "Sequence 2 Region: " << lastAlignment.start2 << " - " << lastAlignment.end2 << "\n\n";

    ss << "Complete alignment:\n";

    int blockSize = 80;
    int totalLength = lastAlignment.sequence1.length();
    int blockNumber = 1;

    for (int start = 0; start < totalLength; start += blockSize) {
        int end = std::min(start + blockSize, totalLength);

        ss << "Block " << blockNumber << ":\n";
        ss << "Seq1: " << lastAlignment.sequence1.substr(start, end - start) << "\n";
        ss << "      " << lastAlignment.matchLine.substr(start, end - start) << "\n";
        ss << "Seq2: " << lastAlignment.sequence2.substr(start, end - start) << "\n";

        int blockMatches = std::count(lastAlignment.matchLine.begin() + start, lastAlignment.matchLine.begin() + end, '|');

        double blockIdentity =
            (end - start > 0) ? (blockMatches * 100.0 / (end - start)) : 0.0;

        ss << "      [Block identity: " << std::fixed << std::setprecision(1) << blockIdentity << "%]\n\n";

        blockNumber++;
    }

    int totalMatches = std::count(lastAlignment.matchLine.begin(), lastAlignment.matchLine.end(), '|');

    ss << "OVERALL STATISTICS:\n";
    ss << "Total Alignment Length: " << totalLength << " bases\n";
    ss << "Matches: " << totalMatches << "\n";
    ss << "Gaps: "
       << std::count(lastAlignment.sequence1.begin(), lastAlignment.sequence1.end(), '-')
        + std::count(lastAlignment.sequence2.begin(), lastAlignment.sequence2.end(), '-') << "\n";

    double overallIdentity = totalMatches * 100.0 / totalLength;

    ss << "Overall Identity: " << std::fixed << std::setprecision(1) << overallIdentity << "%\n";

    return QString::fromStdString(ss.str());
}

void Controller::updateTrieIndex() {
    if (currentSequence) {
        trieIndex.buildIndex(currentSequence->getSequence());
        cout << "Trie index updated for sequence of length " << currentSequence->length() << endl;
    } else {
        cout << "Cannot update trie, there is no sequence" << endl;
    }
}
