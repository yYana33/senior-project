#pragma once
#include "Feature.h"

class Gene : public Feature {
private:
    int readingFrame;
    double gcContent;

public:
    Gene(int start, int end, int readingFrame, double gcContent = 0.0, const std::string& description = "");

    int getReadingFrame() const { return readingFrame; }
    double getGCContent() const { return gcContent; }

    std::string getDisplayName() const override;
    std::string getDetails() const override;
};