#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include "vector3.h"
#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>

#define l2norm(x, y, z) sqrt((x)*(x)+(y)*(y)+(z)*(z))

class PowerVertex;
class PowerEdge;

class GraphNode {
public:
    PowerVertex * pVert;
    uint index;
    double x, y, z;
    double radius;
    bool boundary;
    std::vector<uint> edges;
};

class GraphEdge{
public:
    uint v1, v2;
    uint index;
    PowerEdge * pEdge;
    double weight;

    double length(std::vector<GraphNode> * nodes){
        GraphNode n1 = nodes->at(v1);
        GraphNode n2 = nodes->at(v2);
        return l2norm(n1.x-n2.x, n1.y-n2.y, n1.z-n2.z);
    }

    Vector3 direction(std::vector<GraphNode> * nodes){
        GraphNode n1 = nodes->at(v1);
        GraphNode n2 = nodes->at(v2);
        Vector3 vec (n2.x-n1.x, n2.y-n1.y, n2.z-n1.z);
        vec.Normalize();
        return vec;
    }
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
    double runDijkstra(int start, int target, std::vector<GraphNode*> *pathNodes,
                       std::vector<GraphEdge*> *pathEdges);
    void writeGraph(const char* file);
    void writeGraphCRD(const char* file);
};

#endif // GRAPH_H
