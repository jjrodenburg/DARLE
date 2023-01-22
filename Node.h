//
// Created by Jort Rodenburg on 12/4/22.
//

#ifndef ALGORITHM_NODE_H
#define ALGORITHM_NODE_H

#include "LabPixel.h"

class Node {
public:
    virtual std::string ToString();
};

class AreaNode : public Node {
public:
    int Location;
    double AvgEdelta;
    CIEDE2000::LAB LabEncoding;

    AreaNode(int location, int avgEdelta, CIEDE2000::LAB labEncoding): Location(location), AvgEdelta(avgEdelta),
                                                                       LabEncoding(labEncoding){};
    AreaNode() = default;

    void AddEdge(Node origin, Node destination);
};

class Edge {
public:
    Node source;
    Node destination;

    Edge(Node origin, Node target) : source(origin), destination(target){}
};


#endif //ALGORITHM_NODE_H
