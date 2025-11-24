#pragma once
#include <string>

class Feature {
protected:
    int start;
    int end;
    std::string type;
    std::string description;

public:

    Feature(int start, int end, const std::string& type, const std::string& description = "");
    virtual ~Feature() = default;

    int getStart() const { return start; }
    int getEnd() const { return end; }
    const std::string& getType() const { return type; }
    const std::string& getDescription() const { return description; }
    int getLength() const { return end - start + 1; }

    //methods for polymorphism
    virtual std::string getDisplayName() const = 0;
    virtual std::string getDetails() const = 0;
};


