#pragma once
#include <QObject>
#include <QString>
#include <memory>
#include "./core/DNASequence.h"
#include "./algorithms/TrieIndex.h"
#include "./algorithms/SearchAlgorithm.h"
#include "./algorithms/OrfFinder.h"
#include "./features/Gene.h"
#include "./features/SearchMatch.h"
#include "../io/FastaParser.h"
#include "./algorithms/Aligner.h"
#include "../gui/SequenceVisualizer.h"

class Controller : public QObject {
    Q_OBJECT

signals:
    void errorOccurred(const QString& errorMessage);
    void searchCompleted(const QString& pattern, size_t matchCount);

public:
    const DNASequence* getSequencePointer() const {
        return currentSequence.get();
    }

    std::vector<int> getSearchMatchPositions() const {
        std::vector<int> positions;
        if (!currentSequence) return positions;

        const auto& features = currentSequence->getFeatures();
        for (const auto& feature : features) {
            if (feature->getType() == "search_match") {
                positions.push_back(feature->getStart());
            }
        }
        return positions;
    }

    struct ORFRegion {
        int start;
        int end;
        int strand;   //+1 or -1
    };

    int getPatternLength() const {
        return currentPattern.length();
    }

    std::vector<ORFRegion> getORFs() const;

    struct AlignBlock {
        int start;
        int end;
        double similarity;
    };

    std::vector<SequenceVisualizer::AlignBlock> getAlignmentBlocks() const {
        return alignmentBlocks;
    }

    //block computation
    std::vector<SequenceVisualizer::AlignBlock> computeAlignmentBlocks(
        const QString& seq1Aligned,
        const QString& seq2Aligned);

    explicit Controller(QObject* parent = nullptr);

    double getGCPercentage() const;

    bool loadSequence(const QString& filePath);
    bool loadSecondSequence(const QString& filePath);
    void runAutoSearch();
    void searchPattern(const QString& pattern);
    void findORFs();
    void alignSequences();
    void clearResults();

    int getSearchMatchCount() const;
    int getGeneCount() const;

    void clearSearchResults();
    void clearORFResults();
    void clearAlignmentResults();
    void clearAllResults();

    bool hasSearchResults() const;
    bool hasORFResults() const;
    bool hasAlignmentResults() const;

    QString getStatisticsInfo() const;
    QString getSequenceInfo() const;
    QString getFeaturesInfo() const;
    QString getAlignmentInfo() const;

    bool hasSequence() const { return currentSequence != nullptr; }
    bool hasSecondSequence() const { return secondSequence != nullptr; }

    QString getAutoSearchResults() const { return autoSearchResults; }
    const DNASequence* getCurrentSequence() const { return currentSequence.get(); }

private:
    QString currentPattern;
    std::unique_ptr<DNASequence> currentSequence;
    std::unique_ptr<DNASequence> secondSequence;
    std::vector<std::pair<std::string, int>> baseCounts;
    void calculateStatistics();

    Alignment lastAlignment;
    bool hasAlignmentResult = false;

    //alignment bar
    std::vector<SequenceVisualizer::AlignBlock> alignmentBlocks;

    TrieIndex trieIndex{3};
    Aligner aligner;
    QString autoSearchResults;

    void updateTrieIndex();
};
