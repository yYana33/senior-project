#include "SequenceVisualizer.h"
#include <QToolTip>
#include <algorithm>
#include <cmath>

SequenceVisualizer::SequenceVisualizer(QWidget* parent) : QWidget(parent)
{
    setMinimumHeight(180);
    setMouseTracking(true);
}

void SequenceVisualizer::setSequence(const DNASequence* seq) {
    currentSequence = seq;
    startPosition = 0;
    updateVisibleRange();
    update();
}

void SequenceVisualizer::setSearchMatches(const std::vector<int>& pos, int len) {
    searchMatchPositions = pos;
    patternLength = len;
    update();
}

void SequenceVisualizer::setORFs(const std::vector<ORFRegion>& orfs) {
    orfRegions = orfs;
    update();
}

void SequenceVisualizer::setAlignmentBlocks(const std::vector<AlignBlock>& blocks) {
    alignmentBlocks = blocks;
    update();
}

void SequenceVisualizer::clearVisualization() {
    searchMatchPositions.clear();
    orfRegions.clear();
    alignmentBlocks.clear();
    patternLength = 0;
    update();
}

void SequenceVisualizer::updateMode(bool hasMatches, bool hasORFs, bool hasAlignment) {
    VisualizationMode newMode;

    //Search/ORF results override alignment
    if (hasMatches || hasORFs) {
        newMode = FEATURE_VIEW;
    }
    else if (hasAlignment) {
        newMode = ALIGNMENT_VIEW;
    }
    else {
        newMode = DETAILED_VIEW;
    }

    if (mode != newMode) {
        mode = newMode;
        emit viewModeChanged(mode);
    }

    update();
}

void SequenceVisualizer::setZoomLevel(double zoom) {
    zoomLevel = std::clamp(zoom, 0.2, 20.0);
    updateVisibleRange();
    emit zoomLevelChanged(zoomLevel);
    update();
}

void SequenceVisualizer::setStartPosition(int position) {
    if (!currentSequence) return;

    int seqLen = currentSequence->length();
    int maxStart = std::max(0, seqLen - visibleBases);

    startPosition = std::clamp(position, 0, maxStart);
    update();
}

void SequenceVisualizer::paintEvent(QPaintEvent* event) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    painter.fillRect(rect(), Qt::white);

    if (!currentSequence) {
        painter.drawText(rect(), Qt::AlignCenter, "No sequence loaded");
        return;
    }

    int margin = 10;
    int legendHeight = 30;     //to reserve space for legend
    int gap = 12;              //space between visualization and legend

    int visualizationHeight = height() - legendHeight - gap - margin * 2;

    //main visualization based on mode
    switch (mode) {
    case DETAILED_VIEW:
        drawDetailedView(painter, margin, margin, width() - margin * 2, visualizationHeight);
        break;

    case FEATURE_VIEW:
        drawFeatureView(painter, margin, margin, width() - margin * 2, visualizationHeight);
        break;

    case ALIGNMENT_VIEW:
        drawAlignmentView(painter, margin, margin, width() - margin * 2, visualizationHeight);
        break;
    }

    drawLegend(painter, margin, margin + visualizationHeight + gap);
}

void SequenceVisualizer::drawDetailedView(QPainter& painter, int x, int y, int w, int h) {
    int seqLen = currentSequence->length();
    if (seqLen == 0) return;

    int baseW = std::max(1, w / visibleBases);

    painter.setPen(Qt::NoPen);

    for (int i = 0; i < visibleBases; ++i) {
        int pos = startPosition + i;
        if (pos >= seqLen) break;

        char c = currentSequence->getSequence()[pos];
        QColor col = Qt::gray;

        switch (toupper(c)) {
        case 'A': col = baseA; break;
        case 'T': col = baseT; break;
        case 'C': col = baseC; break;
        case 'G': col = baseG; break;
        }

        painter.setBrush(col);
        int xx = x + i * baseW;
        painter.drawRect(xx, y, baseW, h);
    }
}

void SequenceVisualizer::drawFeatureView(QPainter& painter, int x, int y, int w, int h) {
    painter.setBrush(QColor(220,220,220));
    painter.setPen(Qt::NoPen);
    painter.drawRect(x, y + h/2 - 5, w, 10);

    double scale = (double)w / visibleBases;

    //search matches
    painter.setBrush(searchColor);
    for (int pos : searchMatchPositions) {
        if (pos < startPosition || pos > startPosition + visibleBases) continue;
        int xx = x + (int)((pos - startPosition) * scale);
        painter.drawRect(xx, y + h/2 - 10, 6, 20);
    }

    //ORFs
    int forwardY = y + h/2 + 20;    //lower lane
    int reverseY = y + h/2 - 35;    //upper lane

    QColor orfFill  = QColor(0, 180, 0, 120);
    QColor orfBorder = QColor(0, 120, 0);

    for (const auto& orf : orfRegions) {
        int start = std::max(orf.start, startPosition);
        int end   = std::min(orf.end, startPosition + visibleBases);

        if (end <= start) continue;

        int xx = x + int((start - startPosition) * scale);
        int ww = int((end - start) * scale);

        int laneY = (orf.strand == +1 ? forwardY : reverseY);

        painter.setBrush(orfFill);
        painter.setPen(Qt::NoPen);
        painter.drawRect(xx, laneY, ww, 12);

        painter.setBrush(Qt::NoBrush);
        painter.setPen(QPen(orfBorder, 2));
        painter.drawRect(xx, laneY, ww, 12);

        painter.setPen(QPen(orfBorder, 2));
        if (orf.strand == +1) {
            // forward
            painter.drawLine(xx, laneY + 6, xx + ww, laneY + 6);
            painter.drawLine(xx + ww - 6, laneY + 3, xx + ww, laneY + 6);
        } else {
            // reverse
            painter.drawLine(xx, laneY + 6, xx + ww, laneY + 6);
            painter.drawLine(xx + 6, laneY + 3, xx, laneY + 6);
        }
    }
}

void SequenceVisualizer::drawAlignmentView(QPainter& painter, int x, int y, int w, int h)
{
    if (alignmentBlocks.empty()) return;

    int seqLen = currentSequence->length();
    if (seqLen <= 0) return;

    int midY = y + h/2;
    int barH = 12;

    painter.setPen(QPen(Qt::black, 2));
    painter.drawLine(x, midY, x + w, midY);

    for (const auto& block : alignmentBlocks) {
        double s = block.similarity;
        int sx = x + (block.start / (double)seqLen) * w;
        int ex = x + (block.end   / (double)seqLen) * w;
        int ww = std::max(2, ex - sx);

        QColor c;
        if (s >= 0.90) c = QColor(0,180,0);
        else if (s >= 0.75) c = QColor(220,180,0);
        else if (s >= 0.50) c = QColor(255,140,0);
        else c = QColor(200,0,0);

        painter.setBrush(c);
        painter.setPen(Qt::NoPen);
        painter.drawRect(sx, midY - barH/2, ww, barH);

        painter.setPen(Qt::black);
    }
}

void SequenceVisualizer::drawLegend(QPainter& painter, int startX, int startY) {
    painter.setPen(Qt::black);
    painter.setFont(QFont("Arial", 8));

    int box = 14;
    int spacing = 10;

    int row1Y = startY;
    int row2Y = startY + 22;   //second line below

    // row 1
    QStringList baseLabels = {"A", "T", "C", "G"};
    QColor baseColorsList[4] = { baseA, baseT, baseC, baseG };

    int x = startX;
    for (int i = 0; i < 4; i++) {
        painter.setBrush(baseColorsList[i]);
        painter.setPen(Qt::NoPen);
        painter.drawRect(x, row1Y, box, box);

        painter.setPen(Qt::black);
        painter.drawText(x + box + 3, row1Y + box - 2, baseLabels[i]);

        x += 60; // spacing between A, T, C, G groups
    }

    // row 2
    x = startX;

    // Search matches
    painter.setBrush(searchColor);
    painter.setPen(Qt::NoPen);
    painter.drawRect(x, row2Y, box, box);
    painter.setPen(Qt::black);
    painter.drawText(x + box + 3, row2Y + box - 2, "Search Match");
    x += 120;

    // ORFs
    painter.setBrush(QColor(0, 180, 0, 150));
    painter.setPen(Qt::NoPen);
    painter.drawRect(x, row2Y, box, box);
    painter.setPen(Qt::black);
    painter.drawText(x + box + 3, row2Y + box - 2, "ORF (Forward/Reverse)");
    x += 145;

    // Alignment blocks
    painter.setBrush(QColor(70, 130, 180, 150));
    painter.setPen(Qt::NoPen);
    painter.drawRect(x, row2Y, box, box);
    painter.setPen(Qt::black);
    painter.drawText(x + box + 3, row2Y + box - 2, "Alignment Similarity: 90–100% - Green, 70–90% - Yellow, 50-70% - Orange, 0–50% - Red");
}

void SequenceVisualizer::updateVisibleRange() {
    if (!currentSequence) {
        visibleBases = 0;
        return;
    }

    int seqLen = currentSequence->length();
    int maxVisible = width() / 3;
    int scaled = (int)(maxVisible * zoomLevel);

    visibleBases = std::min(seqLen, scaled);
}

void SequenceVisualizer::resizeEvent(QResizeEvent*) {
    updateVisibleRange();
}

int SequenceVisualizer::xToPos(int x, int totalWidth) const {
    if (!currentSequence || visibleBases <= 0 || totalWidth <= 0)
        return startPosition;

    double rel = double(x) / totalWidth;
    if (rel < 0.0) rel = 0.0;
    if (rel > 1.0) rel = 1.0;

    int pos = startPosition + int(rel * visibleBases);

    if (pos < 0) pos = 0;
    int maxPos = currentSequence->length() - 1;
    if (pos > maxPos) pos = maxPos;
    return pos;
}

int SequenceVisualizer::xToBaseIndex(int x) const {
    if (!currentSequence || visibleBases == 0)
    return -1;

    int margin = 10;
    int drawableWidth = width() - margin * 2;

    double rel = double(x - margin) / drawableWidth;
    if (rel < 0 || rel > 1) return -1;

    return startPosition + int(rel * visibleBases);
}


void SequenceVisualizer::mousePressEvent(QMouseEvent* event) {
    if (!currentSequence) return;
    int pos = xToPos(event->pos().x(), width());
    setStartPosition(pos - visibleBases/2);
}

void SequenceVisualizer::mouseMoveEvent(QMouseEvent* event) {
    if (!currentSequence) return;

    int baseIndex = xToBaseIndex(event->pos().x());
    if (baseIndex < 0 || baseIndex >= currentSequence->length())
    return;

    char base = currentSequence->getSequence()[baseIndex];
    QString tooltip;

    switch (mode) {
    case DETAILED_VIEW:
        tooltip = QString("Position: %1\nBase: %2")
        .arg(baseIndex)
        .arg(base);
        break;

    case FEATURE_VIEW:
        tooltip = QString("Position: %1").arg(baseIndex);
        break;

    case ALIGNMENT_VIEW: {
        tooltip = QString("Position: %1").arg(baseIndex);

        //block similarity info
        for (const auto& b : alignmentBlocks) {
            if (baseIndex >= b.start && baseIndex <= b.end) {
                tooltip += QString("\nSimilarity: %1%")
                .arg(int(b.similarity * 100));
                break;
            }
        }
        break;
    }
    }

    QToolTip::showText(event->globalPos(), tooltip, this);
}

void SequenceVisualizer::wheelEvent(QWheelEvent* event) {
    if (!currentSequence) return;
    double factor = (event->angleDelta().y() > 0) ? 1.1 : 1/1.1;
    setZoomLevel(zoomLevel * factor);
}
