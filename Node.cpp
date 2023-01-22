//
// Created by Jort Rodenburg on 12/4/22.
//

#include "Node.h"
#include <vector>

std::string Node::ToString() {
    return "";
}

/// AREA NODE
std::vector<Edge> edges = std::vector<Edge>();
void AreaNode::AddEdge(Node origin, Node destination) {
    Edge edge = Edge(origin, destination);
    edges.push_back(edge);
}
