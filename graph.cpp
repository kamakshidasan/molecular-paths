#include "graph.h"
#include <stdio.h>
#include <iostream>
#include <powerdiagram.h>
#include <set>

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
        if(i==0 || maxEdgeWeight < gEdge.weight){
            maxEdgeWeight = gEdge.weight;
        }
    }

    initialized = true;
}

void Graph::fillWeightMap(ListGraph::EdgeMap<double>* weightMap){
    for(uint i=0;i<edges.size();i++){
        GraphEdge gEdge = edges[i];
        ListGraph::Edge edge = edgeMap[i];
        (*weightMap)[edge] = gEdge.weight;
    }
}

typedef ListGraph::EdgeMap<double> CostMap;

int Graph::runDijkstra(int start, std::vector<std::vector<GraphNode*> > *pathsNodes,
                        std::vector<std::vector<GraphEdge*> > *pathsEdges){
    ListGraph::NodeMap<double> dist(*lemonGraph);
//    dijkstra(*lemonGraph, *weights).distMap(dist).run(nodeMap[start]);

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
            if(index == 0 || cost<leastCost){
                leastCost = cost;
                shortest = index;
            }
            index++;
        }
    }
    return shortest;
}

bool Graph::runDijkstraEscape(int start, std::vector<GraphNode*> *pathNodes,
                             std::vector<GraphEdge*> *pathEdges){
    return runDijkstraEscapeOneIter(start, pathNodes, pathEdges, weights);
}

bool Graph::runDijkstraEscapeOneIter(int start, std::vector<GraphNode*> *pathNodes,
                             std::vector<GraphEdge*> *pathEdges, ListGraph::EdgeMap<double>* weightMap){
    ListGraph::NodeMap<double> dist(*lemonGraph);

    Dijkstra<ListGraph, CostMap> dijkstra(*lemonGraph, *weightMap);
    dijkstra.distMap(dist);
    dijkstra.init();
    dijkstra.run(nodeMap[start]);

    pathNodes->clear();
    pathEdges->clear();

    double leastCost = 0;
    bool first = true;
    int shortest = -1;

    for(uint i=0;i<nodeMap.size();i++){
        if(nodes[i].boundary && dijkstra.reached(nodeMap[i])){
            double cost = dist[nodeMap[i]];
            if(first){
                leastCost = cost;
                shortest = i;
                first = false;
            }else if(cost<leastCost){
                leastCost = cost;
                shortest = i;
            }
        }
    }

    if(shortest!=-1){
        pathNodes->push_back(&nodes[start]);
        Path<ListGraph> p = dijkstra.path(nodeMap[shortest]);
        for (Path<ListGraph>::ArcIt it(p); it != INVALID; ++it) {
            ListGraph::Arc e = it;
            ListGraph::Node succ = lemonGraph->target(e);
            pathNodes->push_back((*revNodeMap)[succ]);
            pathEdges->push_back((*revEdgeMap)[e]);
        }
        return true;
    }
    return false;
}

bool Graph::runDijkstraEscapeRepeated(int start, int maxIter, std::vector<std::vector<GraphNode*> > *pathsNodes,
                        std::vector<std::vector<GraphEdge*> > *pathsEdges){
    ListGraph::EdgeMap<double> weightMap(*lemonGraph);
    fillWeightMap(&weightMap);

    pathsNodes->clear();
    pathsEdges->clear();

    int iter = 0;
    for(iter=0;iter<maxIter;iter++){
        std::vector<GraphNode*> nodeList;
        std::vector<GraphEdge*> edgeList;
        bool foundPath = runDijkstraEscapeOneIter(start, &nodeList, &edgeList, &weightMap);
        if(!foundPath){
            std::cout << "Escape Path not found -- " << iter << std::endl;
            break;
        }
        pathsNodes->push_back(nodeList);
        pathsEdges->push_back(edgeList);
        // assign very high weights to already found path edges
        for(int i=0;i<edgeList.size();i++){
            GraphEdge* gEdge = edgeList[i];
            weightMap[edgeMap[gEdge->index]] = maxEdgeWeight;
        }
    }
    return iter>0;
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

void Graph::writePathCRD(const char* file, std::vector<GraphNode*> *pathNodes){
    FILE *fp = fopen(file,"w");
    fprintf(fp, "  %d\n", pathNodes->size());
    fprintf(fp, "#   i        X          Y          Z        R         Epsilon     Sigma     Charge      ASP       Atm name    Res name    Chain      Res #\n");
    for(uint i=0; i<pathNodes->size();i++){
        GraphNode* node = pathNodes->at(i);
        double radius = node->radius;
        fprintf(fp, " \t%d \t%f \t%f \t%f \t%f \t0.000 \t0.000 \t0.000 \t0.000 \t%s \tGLY \tA \t1\n", (i+1),
                node->x, node->y, node->z, radius, node->boundary?"C":"N");
    }
    fprintf(fp, "#TER\n");
    fprintf(fp, "#  i       Nexclude    Exclude list ...\n");
    for(uint i=0; i<pathNodes->size()-1;i++){
        fprintf(fp, " \t%d \t1 \t%d\n", (i+1), (i+2));
    }
    if(pathNodes->size()>0){
        fprintf(fp, " \t%d \t0\n", pathNodes->size());
    }
    fclose(fp);
}

void Graph::writeAllPathsCRD(const char* file, std::vector<std::vector<GraphNode*> > *pathsNodes){
    FILE *fp = fopen(file,"w");
    fprintf(fp, "              \n");
    fprintf(fp, "#   i        X          Y          Z        R         Epsilon     Sigma     Charge      ASP       Atm name    Res name    Chain      Res #\n");
    std::vector<std::set<int> > edges;
    std::vector<int> map;
    for(int i=0;i<nodes.size();i++){
        map.push_back(-1);
    }
    int index = 1;
    for(int i=0;i<pathsNodes->size();i++){
        uint size = pathsNodes->at(i).size();
        for(uint j=0; j<size; j++){
            GraphNode* node = pathsNodes->at(i)[j];
            if(map[node->index]==-1){
                double radius = node->radius;
                fprintf(fp, " \t%d \t%f \t%f \t%f \t%f \t0.000 \t0.000 \t0.000 \t0.000 \t%s \tGLY \tA \t1\n", index,
                        node->x, node->y, node->z, radius, j==0? "O": (node->boundary?"C":"N"));
                map[node->index] = index-1;
                std::set<int> edgeIndices;
                edges.push_back(edgeIndices);
                index++;
            }
            if(j!=0){
                GraphNode* prev = pathsNodes->at(i)[j-1];
                edges[map[prev->index]].insert(map[node->index]);
            }
        }
    }
    fprintf(fp, "#TER\n");
    fprintf(fp, "#  i       Nexclude    Exclude list ...\n");
    for(int i=0;i<edges.size();i++){
        fprintf(fp, " \t%d \t%d", (i+1), edges[i].size());
        for (std::set<int>::iterator it=edges[i].begin() ; it != edges[i].end(); it++){
            fprintf(fp, " \t%d", (*it)+1);
        }
        fprintf(fp, "\n");
    }
    fseek(fp, 1, SEEK_SET);
    fprintf(fp, "%d", edges.size());
    fclose(fp);
}

