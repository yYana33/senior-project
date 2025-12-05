#pragma once
#include <QWidget>

class GCContentCircle : public QWidget {
    Q_OBJECT

public:
    explicit GCContentCircle(QWidget* parent = nullptr);

    void setGC(double gc) {
        gcPercent = gc;
        update();
    }

protected:
    void paintEvent(QPaintEvent* event) override;

private:
    double gcPercent = 0.0;
};
