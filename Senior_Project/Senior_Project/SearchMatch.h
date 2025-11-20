#pragma once
#include "Feature.h"

class SearchMatch : public Feature {
private:
    std::string pattern;

public:
    SearchMatch(int start, int end, const std::string& pattern, const std::string& description = "");

    const std::string& getPattern() const { return pattern; }

    std::string getDisplayName() const override;
    std::string getDetails() const override;
};