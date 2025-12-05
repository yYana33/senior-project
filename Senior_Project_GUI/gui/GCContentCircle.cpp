#include "GCContentCircle.h"
#include <QPainter>
#include <QPen>
#include <QFont>

GCContentCircle::GCContentCircle(QWidget* parent)
    : QWidget(parent)
{
    setMinimumSize(120, 120);
}

void GCContentCircle::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    int w = width();
    int h = height();
    int size = qMin(w, h) - 10;

    QRectF rect((w - size) / 2, (h - size) / 2, size, size);

    p.setPen(QPen(QColor(200, 200, 200), 12));
    p.drawArc(rect, 0 * 16, 360 * 16);

    p.setPen(QPen(QColor(50, 150, 255), 12));
    int angle = int((gcPercent / 100.0) * 360);
    p.drawArc(rect, 90 * 16, -angle * 16);

    p.setPen(Qt::black);
    p.setFont(QFont("Arial", 14, QFont::Bold));
    p.drawText(rect, Qt::AlignCenter, QString("%1% GC").arg(gcPercent, 0, 'f', 1));
}
