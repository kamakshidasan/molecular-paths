#include "graph.h"
#include <stdio.h>
#include <iostream>
#include <powerdiagram.h>

using namespace lemon;

Graph::Graph()
{
    try{
        lemonGraph = new ListGraph();
    }catch(...){
        std::cout << "Not initialized" << std::endl;
    }
    revNodeMap = NULL;
    revEdgeMap = NULL;
    weights = NULL;
    initialized = false;
}

void Graph::clear(){
    if(initialized){
        if (revNodeMap!=NULL) {
            delete revNodeMap;
            revNodeMap = NULL;
        }
        if (revEdgeMap!=NULL) {
            delete revEdgeMap;
            revEdgeMap = NULL;
        }
        if (weights!=NULL) {
            delete weights;
            weights = NULL;
        }
        nodeMap.clear();
        edgeMap.clear();
        lemonGraph->clear();
        initialized = false;
    }
    nodes.clear();
    edges.clear();
}

void Graph::initLemonGraph(){

    revNodeMap = new ListGraph::NodeMap<GraphNode*>(*lemonGraph);

    for(uint i=0;i<nodes.size();i++){
        ListGraph::Node node = lemonGraph->addNode();
        nodeMap.push_back(node);
        (*revNodeMap)[node] = &nodes[i];
    }

    revEdgeMap = new ListGraph::EdgeMap<GraphEdge*>(*lemonGraph);
    weights = new ListGraph::EdgeMap<double>(*lemonGraph);

    for(uint i=0;i<edges.size();i++){
        GraphEdge gEdge = edges[i];
        ListGraph::Edge edge = lemonGraph->addEdge(nodeMap[gEdge.v1], nodeMap[gEdge.v2]);
        edgeMap.push_back(edge);
        (*revEdgeMap)[edge] = &edges[i];
        (*weights)[edge] = gEdge.weight;
    }

    initialized = true;
}

typedef ListGraph::EdgeMap<double> CostMap;

int Graph::runDijkstra(int start, std::vector<std::vector<GraphNode*> > *pathsNodes,
                        std::vector<std::vector<GraphEdge*> > *pathsEdges){
    ListGraph::NodeMap<double> dist(*lemonGraph);
//    dijkstra(*lemonGraph, *weights).distMap(dist).run(nodeMap[start]);

    std::cout << nodes[start].edges.size() << std::endl;

    Dijkstra<ListGraph, CostMap> dijkstra(*lemonGraph, *weights);
    dijkstra.distMap(dist);
    dijkstra.init();
    dijkstra.run(nodeMap[start]);

    pathsNodes->clear();
    pathsEdges->clear();

    int shortest = -1;
    double leastCost = 0;
    int index = 0;

    for(uint i=0;i<nodeMap.size();i++){
        if(nodes[i].boundary && dijkstra.reached(nodeMap[i])){
            Path<ListGraph> p = dijkstra.path(nodeMap[i]);
            std::vector<GraphNode*> nodeList;
            std::vector<GraphEdge*> edgeList;
            nodeList.push_back(&nodes[start]);
            for (Path<ListGraph>::ArcIt it(p); it != INVALID; ++it) {
                ListGraph::Arc e = it;
                ListGraph::Node succ = lemonGraph->target(e);
                nodeList.push_back((*revNodeMap)[succ]);
                edgeList.push_back((*revEdgeMap)[e]);
            }
            pathsNodes->push_back(nodeList);
            pathsEdges->push_back(edgeList);
            double cost = dist[nodeMap[i]];
            if(index == 0){
                leastCost = cost;
                shortest = index;
            }else if(cost<leastCost){
                leastCost = cost;
                shortest = index;
            }
        }
    }

    return shortest;
}

double Graph::runDijkstra(int start, int target, std::vector<GraphNode*> *pathNodes,
                          std::vector<GraphEdge*> *pathEdges){
    if(start<0 || start>=nodes.size() || target<0 || target>=nodes.size()){
        return -1;
    }

    ListGraph::NodeMap<double> dist(*lemonGraph);

    Dijkstra<ListGraph, CostMap> dijkstra(*lemonGraph, *weights);
    dijkstra.distMap(dist);
    dijkstra.init();
    dijkstra.run(nodeMap[start], nodeMap[target]);

    pathNodes->clear();
    pathEdges->clear();

    if(dijkstra.reached(nodeMap[target])){
        pathNodes->push_back(&nodes[start]);
        Path<ListGraph> p = dijkstra.path(nodeMap[target]);
        for (Path<ListGraph>::ArcIt it(p); it != INVALID; ++it) {
            ListGraph::Arc e = it;
            ListGraph::Node succ = lemonGraph->target(e);
            pathNodes->push_back((*revNodeMap)[succ]);
            pathEdges->push_back((*revEdgeMap)[e]);
        }
        return dist[nodeMap[target]];
    } else {
        return -1;
    }
}

void Graph::writeGraph(const char* file){
    FILE *fp = fopen(file,"w");
    fprintf(fp, "%d %d\n", nodes.size(), edges.size());
    for(uint i=0; i<nodes.size();i++){
        GraphNode node = nodes[i];
        fprintf(fp, "%f %f %f %f %d\n", node.pVert->center.X, node.pVert->center.Y,
                node.pVert->center.Z, node.pVert->powerDistance, node.boundary);
    }
    for(uint i=0; i<edges.size();i++){
        GraphEdge edge = edges[i];
        fprintf(fp, "%d %d %f\n", edge.v1, edge.v2, edge.weight);
    }
    fclose(fp);
}


void Graph::writeGraphCRD(const char* file){
    FILE *fp = fopen(file,"w");
    fprintf(fp, "  %d\n", nodes.size());
    fprintf(fp, "#   i        X          Y          Z        R         Epsilon     Sigma     Charge      ASP       Atm name    Res name    Chain      Res #\n");
    for(uint i=0; i<nodes.size();i++){
        GraphNode node = nodes[i];
        double radius = node.radius;
//        if (node.boundary || radius>2){
//            radius = 0.1;
//        }
        fprintf(fp, " \t%d \t%f \t%f \t%f \t%f \t0.000 \t0.000 \t0.000 \t0.000 \t%s \tGLY \tA \t1\n", (i+1),
                node.x, node.y, node.z, radius, node.boundary?"C":"N");
    }
    fprintf(fp, "#TER\n");
    fprintf(fp, "#  i       Nexclude    Exclude list ...\n");
    for(uint i=0; i<nodes.size();i++){
        GraphNode node = nodes[i];
        int numEdges = 0;
        for (uint j=0;j<node.edges.size(); j++){
            GraphEdge edge = edges[node.edges[j]];
            uint other = (i == edge.v1) ? edge.v2 : edge.v1;
            if(other>i) numEdges++;
        }
        fprintf(fp, " \t%d \t%d", (i+1), numEdges);
        for (uint j=0;j<node.edges.size(); j++){
            GraphEdge edge = edges[node.edges[j]];
            uint other = (i == edge.v1) ? edge.v2 : edge.v1;
            if(other>i){
                other++;
                fprintf(fp, " \t%d", other);
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
