#include "TrieNode.h"

TrieNode::TrieNode() { }

TrieNode::~TrieNode() {
    for (auto& pair : children) {
        delete pair.second;
    }
}

void TrieNode::addChild(char nucleotide) {
    if (children.find(nucleotide) == children.end()) {
        children[nucleotide] = new TrieNode(); //creating new node if it doesn't exist
    }
}

TrieNode* TrieNode::getChild(char nucleotide) {
    auto it = children.find(nucleotide);
    if (it != children.end()) {
        return it->second;
    }
    return nullptr;
}

bool TrieNode::hasChild(char nucleotide) {
    return children.find(nucleotide) != children.end();
}

void TrieNode::addPosition(int position) {
    positions.push_back(position);
}

const std::vector<int>& TrieNode::getPositions() const {
    return positions;
}

const std::unordered_map<char, TrieNode*>& TrieNode::getChildren() const {
    return children;
}