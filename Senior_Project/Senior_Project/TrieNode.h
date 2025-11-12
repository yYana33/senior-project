#pragma once
#include <unordered_map>
#include <vector>

class TrieNode {
private:
    std::unordered_map<char, TrieNode*> children;
    std::vector<int> positions; //positions where k-mer ends

public:
    TrieNode();
    ~TrieNode();

    void addChild(char nucleotide);//adding a child node for a character
    TrieNode* getChild(char nucleotide);

    bool hasChild(char nucleotide);

    void addPosition(int position);

    const std::vector<int>& getPositions() const;

    const std::unordered_map<char, TrieNode*>& getChildren() const;//for traversal
}; 
