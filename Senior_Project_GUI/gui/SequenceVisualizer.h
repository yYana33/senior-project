#pragma once
#include <QWidget>
#include <QPainter>
#include <QMouseEvent>
#include <vector>
#include "./core/DNASequence.h"

class SequenceVisualizer : public QWidget {
    Q_OBJECT

public:

    //visualization modes
    enum VisualizationMode {
        DETAILED_VIEW,
        FEATURE_VIEW,
        ALIGNMENT_VIEW
    };

    //alignment block structure
    struct AlignBlock {
        int start;          // position in first seq
        int end;            // position in first seq
        double similarity;  // from 0.0 to 1.0
    };

    struct ORFRegion {
        int start;
        int end;
        int strand;   // +1 = forward, -1 = reverse
    };

    explicit SequenceVisualizer(QWidget* parent = nullptr);

    void setSequence(const DNASequence* sequence);
    void setSearchMatches(const std::vector<int>& positions, int patternLength);
    void setORFs(const std::vector<ORFRegion>& orfs);

    void setAlignmentBlocks(const std::vector<AlignBlock>& blocks);

    void clearVisualization();

    //navigation/zoom
    void setZoomLevel(double zoom);
    void setStartPosition(int position);
    double getZoomLevel() const { return zoomLevel; }

    //mode switch
    void updateMode(bool hasMatches, bool hasORFs, bool hasAlignment);

signals:
    void zoomLevelChanged(double zoom);
    void viewModeChanged(VisualizationMode mode);

protected:
    void paintEvent(QPaintEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;

private:
    const DNASequence* currentSequence = nullptr;
    std::vector<int> searchMatchPositions;
    int patternLength = 0;
    std::vector<ORFRegion> orfRegions;
    std::vector<AlignBlock> alignmentBlocks;

    //display state
    VisualizationMode mode = DETAILED_VIEW;
    double zoomLevel = 1.0;
    int startPosition = 0;
    int visibleBases = 0;
    int xToBaseIndex(int x) const;
    void mouseMoveEvent(QMouseEvent* event) override;

    void drawDetailedView(QPainter& painter, int x, int y, int w, int h);
    void drawFeatureView(QPainter& painter, int x, int y, int w, int h);
    void drawAlignmentView(QPainter& painter, int x, int y, int w, int h);
    void drawLegend(QPainter& painter, int x, int y);

    void updateVisibleRange();
    int xToPos(int x, int totalWidth) const;

    //colours - UCSC based
    QColor baseA = QColor("#00A000");
    QColor baseT = QColor("#DD0000");
    QColor baseC = QColor("#0000CC");
    QColor baseG = QColor("#FF9900");

    QColor searchColor = QColor(255, 0, 0, 180);
    QColor orfColor = QColor(0, 180, 0, 150);
};
