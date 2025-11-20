#include "SearchMatch.h"
#include <sstream>

using namespace std;

SearchMatch::SearchMatch(int start, int end, const std::string& pattern, const std::string& description) : Feature(start, end, "search_match", description), pattern(pattern) {}

string SearchMatch::getDisplayName() const {
    return "Search: '" + pattern + "'";
}

string SearchMatch::getDetails() const {
    stringstream ss;

    ss << "Pattern '" << pattern << "' found at " << start << "-" << end << " (Length: " << getLength() << " bp)";

    return ss.str();
}