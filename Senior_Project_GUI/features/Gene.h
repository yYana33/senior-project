#pragma once
#include "./core/Feature.h"

class Gene : public Feature {
private:
    int readingFrame;
    double gcContent;

public:
    Gene(int start, int end, int readingFrame, int strand, double gcContent = 0.0, const std::string& description = "");

    int getReadingFrame() const { return readingFrame; }
    double getGCContent() const { return gcContent; }
    int getStrand() const { return strand; }

    std::string getDisplayName() const override;
    std::string getDetails() const override;
};
