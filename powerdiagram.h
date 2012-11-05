#ifndef POWERDIAGRAM_H
#define POWERDIAGRAM_H

#include <vector>
#include <QtOpenGL>
#include "deluanaycomplex.h"
#include "tetrahedron.h"
#include "triangle.h"
#include "graph.h"
#include <iostream>

class PowerVertex{
public:
    int tetIndex;
    Vector3 center;
    double powerDistance;
    std::vector<int> neigbours;
    bool inside;
};

const uint INSIDE = 1;
const uint INTERSECTING = 2;
const uint OUTSIDE = 0;
const uint INFINITE = 4;

class PowerEdge{
public:
    int v1, v2;
    int triIndex;
    uint edgeType;
    Vector3 intersect, triIntersect;
    double leastPowerDistance;
    double weight;

    double edgeLength(std::vector <PowerVertex> *vertices, bool useIntersect){
        PowerVertex pv1 = vertices->at(v1);
        Vector3 p1 = pv1.center;
        Vector3 p2 = intersect;
        if (edgeType == INSIDE) useIntersect = false;
        if (!useIntersect && !((edgeType & INFINITE)==INFINITE)){
            p2 = vertices->at(v2).center;
        }
        if(useIntersect && (edgeType == INTERSECTING)){
            if(!pv1.inside){
                p1 = vertices->at(v2).center;
            }
        }
        double len;
        Vector3 diff;
        Vector3::DiffVector(&diff, &p1, &p2);
        Vector3::DotProduct(&diff, &diff, &len);
        len = sqrt(len);
        return len;
    }

    void setLeastPowerDistance(double leastPD, std::vector <PowerVertex> *vertices, bool useIntersect){
        leastPowerDistance = leastPD;
        double length = edgeLength(vertices, useIntersect);
        weight = leastPowerDistance == 0? 9999999 : (length / leastPowerDistance);
    }

    int getOtherVertex(int v){
        return v == v1? v2: v1;
    }

    bool isIncident(int v){
        return v==v1 || v==v2;
    }
};

class PowerDiagram
{
private:
    double min[3];
    double max[3];

public:
    std::vector<PowerVertex> vertices;
    std::vector<PowerEdge> edges;
    DeluanayComplex* delCplx;
    std::vector <Vertex> & vertList;
    Graph currentGraph;
    std::vector<GraphNode*> pathNodes;
    std::vector<std::vector<GraphNode*> > pathsNodes;

    bool singlePath;

    PowerDiagram(DeluanayComplex* delCplx, std::vector <Vertex> &vertlist,
                 double min[], double max[]);
    void printGraph();
    void makeDisplayList(bool complementSpacePD, bool onlyInsideVerts, bool pruneIsolatedVerts, bool intersectEdges);
    void render(bool showPath);
    int getEdgeTo(int v1, int v2);
    void constructGraph(bool considerAlpha);
    void writeGraph(bool considerAlpha, const char* filename);
    void savePathCRD(const char* filename);
    bool findShortestPath(int start, int end, QVector<double>* X, QVector<double>* Y,
                          double * length, double *minY, double *maxY);
    int findShortestEscapePaths(int start, int steps, bool repeated, int maxIter, std::vector<QVector<double> >* Xs,
                                 std::vector<QVector<double> >* Ys,
                                 std::vector<double> * lengths, std::vector<double> *minYs,
                                 std::vector<double> *maxYs);
    bool findShortestEscapePath(int start,QVector<double>* X, QVector<double>* Y,
                                              double * length, double *minY, double *maxY);
    void getPathWeights(std::vector<GraphNode*> *pathNodes,  std::vector<GraphEdge*> *pathEdges,
                        QVector<double>* X, QVector<double>* Y, double * length, double *minY, double *maxY);
    ~PowerDiagram();
};

#endif // POWERDIAGRAM_H
