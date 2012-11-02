#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include "vector3.h"
#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>

class PowerVertex;
class PowerEdge;

class GraphNode {
public:
    PowerVertex * pVert;
    double x, y, z;
    double radius;
    bool boundary;
    std::vector<uint> edges;
};

class GraphEdge{
public:
    uint v1, v2;
    PowerEdge * pEdge;
    double weight;
};

using namespace lemon;

class Graph
{
public:
    std::vector<GraphNode> nodes;
    std::vector<GraphEdge> edges;
    std::vector<ListGraph::Node> nodeMap;
    std::vector<ListGraph::Edge> edgeMap;
    ListGraph::NodeMap<GraphNode*> * revNodeMap;
    ListGraph::EdgeMap<GraphEdge*> * revEdgeMap;

    ListGraph::EdgeMap<double> * weights;
    ListGraph* lemonGraph;

    bool initialized;

    Graph();

    void clear();
    void initLemonGraph();
    void runDijkstra(int start, std::vector<double> *cost);
    double runDijkstra(int start, int target, std::vector<GraphNode*> *path);
    void writeGraph(const char* file);
    void writeGraphCRD(const char* file);
};

#endif // GRAPH_H
