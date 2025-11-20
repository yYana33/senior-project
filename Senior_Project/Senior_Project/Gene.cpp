#include "Gene.h"
#include <sstream>

using namespace std;

Gene::Gene(int start, int end, int readingFrame, double gcContent, const std::string& description) : Feature(start, end, "gene", description), readingFrame(readingFrame), gcContent(gcContent) {}

string Gene::getDisplayName() const {
    stringstream ss;
    ss << "Gene (Frame " << readingFrame << ")";
    return ss.str();
}

string Gene::getDetails() const {
    stringstream ss;

    ss << "Gene: " << start << "-" << end << " (Length: " << getLength() << " bp, Frame: " << readingFrame << ", GC: " << (gcContent * 100) << "%)";
    if (!description.empty()) {
        ss << " - " << description;
    }

    return ss.str();
}
