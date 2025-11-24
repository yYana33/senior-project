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
#include "./io/FastaParser.h"
#include "./algorithms/Aligner.h"

class Controller : public QObject {
    Q_OBJECT

public:
    explicit Controller(QObject* parent = nullptr);

    bool loadSequence(const QString& filePath);
    bool loadSecondSequence(const QString& filePath);
    void searchPattern(const QString& pattern);
    void findORFs();
    void alignSequences();
    void clearResults();

    QString getSequenceInfo() const;
    QString getFeaturesInfo() const;
    QString getAlignmentInfo() const;
    bool hasSequence() const { return currentSequence != nullptr; }
    bool hasSecondSequence() const { return secondSequence != nullptr;}

private:
    std::unique_ptr<DNASequence> currentSequence;
    std::unique_ptr<DNASequence> secondSequence;
    Alignment lastAlignment;
    bool hasAlignmentResult = false;
    TrieIndex trieIndex{3};
    Aligner aligner;

    void updateTrieIndex();
};
