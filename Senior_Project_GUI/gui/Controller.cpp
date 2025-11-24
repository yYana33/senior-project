#include "Controller.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <QDebug>
#include <QFileInfo>
#include <QFile>

using namespace std;

Controller::Controller(QObject* parent) : QObject(parent) {
    qDebug() << "Controller created";
}

bool Controller::loadSequence(const QString& filePath) {
    qDebug() << "Attempting to load file:" << filePath;
    qDebug() << "File exists:" << QFile::exists(filePath);

    auto newSequence = FastaParser::parseFromFile(filePath.toStdString());

    if (newSequence) {
        currentSequence = std::move(newSequence);
        qDebug() << "Successfully loaded sequence! Length:" << currentSequence->length();
        updateTrieIndex();
        return true;
    } else {
        qDebug() << "Failed to load sequence!";
        return false;
    }
}

bool Controller::loadSecondSequence(const QString& filePath) {
    qDebug() << "Attempting to load second file for alignment:" << filePath;
    qDebug() << "File exists:" << QFile::exists(filePath);

    auto newSequence = FastaParser::parseFromFile(filePath.toStdString());

    if (newSequence) {
        secondSequence = std::move(newSequence);
        qDebug() << "Successfully loaded second sequence! Length:" << secondSequence->length();
        return true;
    } else {
        qDebug() << "Failed to load second sequence!";
        return false;
    }
}

void Controller::searchPattern(const QString& pattern) {
    if (!currentSequence) {
        qDebug() << "Cannot search - no sequence loaded!";
        return;
    }

    qDebug() << "STARTING SEARCH";
    qDebug() << "Searching for pattern:" << pattern;
    qDebug() << "Current features before search:" << currentSequence->getFeatures().size();

    string patternStr = pattern.toStdString();
    SearchAlgorithm::searchAndAddMatches(*currentSequence, patternStr);

    qDebug() << "Current features after search:" << currentSequence->getFeatures().size();

    const auto& features = currentSequence->getFeatures();
    for (const auto& feature : features) {
        qDebug() << "Feature:" << QString::fromStdString(feature->getDisplayName()) << "at" << feature->getStart() << "-" << feature->getEnd();
    }

    qDebug() << "SEARCH COMPLETED";
}

void Controller::findORFs() {
    if (!currentSequence) {
        qDebug() << "Cannot find ORFs - no sequence loaded!";
        return;
    }

    qDebug() << "STARTING ORF FINDING";
    qDebug() << "Current features before ORF:" << currentSequence->getFeatures().size();

    OrfFinder::findAndAddORFs(*currentSequence, trieIndex);

    qDebug() << "Current features after ORF:" << currentSequence->getFeatures().size();

    const auto& features = currentSequence->getFeatures();
    for (const auto& feature : features) {
        qDebug() << "Feature:" << QString::fromStdString(feature->getDisplayName()) << "at" << feature->getStart() << "-" << feature->getEnd();
    }
    qDebug() << "ORF finding COMPLETED";
}

void Controller::alignSequences() {
    if (!currentSequence || !secondSequence) {
        qDebug() << "Cannot align - need two sequences loaded!";
        return;
    }

    qDebug() << "STARTING ALIGNMENT";
    qDebug() << "Sequence 1 length:" << currentSequence->length();
    qDebug() << "Sequence 2 length:" << secondSequence->length();

    //a progress bar for long sequences
    if (currentSequence->length() > 1000 || secondSequence->length() > 1000) {
        qDebug() << "Performing alignment on large sequences...";
    }

    lastAlignment = aligner.align(currentSequence->getSequence(), secondSequence->getSequence());
    hasAlignmentResult = true;

    qDebug() << "Alignment completed with score:" << lastAlignment.score;
    qDebug() << "Alignment length:" << lastAlignment.sequence1.length();
    qDebug() << "Similarity:" << Aligner::calculateSimilarity(lastAlignment) << "%";

    if (lastAlignment.sequence1.length() > 1000) {
        qDebug() << "Long alignment generated";
    }
}
QString Controller::getAlignmentInfo() const {
    if (!currentSequence || !secondSequence) {
        return "Please load two sequences for an alignment.";
    }

    if (!hasAlignmentResult) {
        return "Two sequences loaded. Use the Align button to compare them.";
    }

    stringstream ss;
    ss << "SEQUENCE ALIGNMENT RESULT\n\n";
    ss << "Alignment Score: " << lastAlignment.score << "\n";
    ss << "Similarity: " << Aligner::calculateSimilarity(lastAlignment) << "%\n";
    ss << "Aligned Length: " << lastAlignment.sequence1.length() << " bases\n";
    ss << "Sequence 1 Region: " << lastAlignment.start1 << " - " << lastAlignment.end1 << "\n";
    ss << "Sequence 2 Region: " << lastAlignment.start2 << " - " << lastAlignment.end2 << "\n\n";

    ss << "Complete alignment:\n";

    // Displaying alignment in blocks of 80 characters
    int blockSize = 80;
    int totalLength = lastAlignment.sequence1.length();
    int blockNumber = 1;

    for (int start = 0; start < totalLength; start += blockSize) {
        int end = min(start + blockSize, totalLength);

        ss << "Block " << blockNumber << ":\n";
        ss << "Seq1: " << lastAlignment.sequence1.substr(start, end - start) << "\n";
        ss << "      " << lastAlignment.matchLine.substr(start, end - start) << "\n";
        ss << "Seq2: " << lastAlignment.sequence2.substr(start, end - start) << "\n";

        int blockMatches = count(lastAlignment.matchLine.begin() + start, lastAlignment.matchLine.begin() + end, '|');
        double blockIdentity = (end - start > 0) ? (blockMatches * 100.0 / (end - start)) : 0.0;

        ss << "      [Block identity: " << fixed << setprecision(1) << blockIdentity << "%]\n\n";
        blockNumber++;
    }

    //overall stats
    int totalMatches = count(lastAlignment.matchLine.begin(), lastAlignment.matchLine.end(), '|');
    int totalGaps = count(lastAlignment.sequence1.begin(), lastAlignment.sequence1.end(), '-') + count(lastAlignment.sequence2.begin(), lastAlignment.sequence2.end(), '-');

    ss << "OVERALL STATISTICS:\n";
    ss << "Total Alignment Length: " << totalLength << " bases\n";
    ss << "Matches: " << totalMatches << "\n";
    ss << "Mismatches: " << (totalLength - totalMatches - totalGaps) << "\n";
    ss << "Gaps: " << totalGaps << "\n";
    ss << "Overall Identity: " << fixed << setprecision(1)
       << (totalMatches * 100.0 / totalLength) << "%\n";
    ss << "Alignment Score: " << lastAlignment.score << "\n";

    return QString::fromStdString(ss.str());
}

void Controller::clearResults() {
    if (currentSequence) {
        currentSequence->clearFeatures();
    }
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
    if (features.empty()) return "No features found";

    stringstream ss;
    ss << "Found " << features.size() << " features:\n\n";

    for (const auto& feature : features) {
        ss << "• " << feature->getDisplayName() << "\n";
        ss << "  " << feature->getDetails() << "\n\n";
    }

    return QString::fromStdString(ss.str());
}

void Controller::updateTrieIndex() {
    if (currentSequence) {
        trieIndex.buildIndex(currentSequence->getSequence());
        cout << "Trie index updated for sequence of length " << currentSequence->length() << endl;
    } else {
        cout << "Cannot update trie - no sequence loaded" << endl;
    }
}
