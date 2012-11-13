/***************************************************************************
 *   Copyright (C) 2010 by talha bin masood                                *                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

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
    bool selected;

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

    uint getOtherVertex(int v){
        return v==v1? v2: v1;
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
    double maxEdgeWeight;

    Graph();

    void clear();
    void initLemonGraph();
    void fillWeightMap(ListGraph::EdgeMap<double>* weightMap);
    int runDijkstra(int start, std::vector<std::vector<GraphNode*> > *pathsNodes,
                     std::vector<std::vector<GraphEdge*> > *pathsEdges);
    double runDijkstra(int start, int target, std::vector<GraphNode*> *pathNodes,
                       std::vector<GraphEdge*> *pathEdges);
    bool runDijkstraEscape(int start, std::vector<GraphNode*> *pathNodes,
                                 std::vector<GraphEdge*> *pathEdges);
    bool runDijkstraEscapeOneIter(int start, std::vector<GraphNode*> *pathNodes,
                                 std::vector<GraphEdge*> *pathEdges, ListGraph::EdgeMap<double>* weightMap);
    bool runDijkstraEscapeRepeated(int start, int maxIter, std::vector<std::vector<GraphNode*> > *pathsNodes,
                            std::vector<std::vector<GraphEdge*> > *pathsEdges);
    void writeGraph(const char* file);
    void writeGraphCRD(const char* file);
    void writeGraphOFF(const char* file);
    void writePathOFF(const char* file, std::vector<GraphNode*> *pathNodes);
    void writePathCRD(const char* file, std::vector<GraphNode*> *pathNodes);
    void writeAllPathsCRD(const char* file, std::vector<std::vector<GraphNode*> > *pathsNodes);
};

#endif // GRAPH_H
