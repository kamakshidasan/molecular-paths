#include "powerdiagram.h"
#include "lightmaterial.h"
#include "metric.h"
#include <list>
#include <GL/glu.h>
#include <skinsurface.h>

bool inside(Vector3 v, double min[], double max[], double pad){
    return (v.X>(min[0]-pad) && v.X<(max[0]+pad) &&
            v.Y>(min[1]-pad) && v.Y<(max[1]+pad) &&
            v.Z>(min[2]-pad) && v.Z<(max[2]+pad));
}

std::vector<Triangle> getConvexHull(DeluanayComplex* delcx, std::vector<Vertex> &vertlist){
    std::vector<Triangle> cHull;
    for(int i=1; i< delcx->DeluanayTrigs.size(); i++){
        Triangle tri = delcx->DeluanayTrigs[i];
        if(tri.Hull){
           cHull.push_back(tri);
        }
    }

    // orient convex hull triangles

    for(int k=0;k<cHull.size();k++){
        Triangle tri = cHull[k];
        if(tri.ReverseLink1==0){
            std::cout << "Reverse Link of triangle not initialized : " << std::endl;
            continue;
        }
        Tetrahedron tet = delcx->DeluanayTet[tri.ReverseLink1];
        int otherVertIndex = -1;
        for(int i=1;i<5;i++){
            bool contains = false;
            for(int j=1;j<4;j++){
                if(tet.Corners[i]==tri.Corners[j]){
                    contains = true;
                    break;
                }
            }
            if(!contains){
                otherVertIndex = i;
                break;
            }
        }
        if(otherVertIndex ==-1){
            std::cout << "Other Vertex not found : " << std::endl;
            continue;
        }
        Vertex other = vertlist[tet.Corners[otherVertIndex]];
        if(tri.Corners[1]==-1){
            std::cout << "Triangle Corner Problem : " << std::endl;
            continue;
        }
        Vertex triVert = vertlist[tri.Corners[1]];
        Vector3 otherV = other.getCoordVector();
        Vector3 triV = triVert.getCoordVector();
        Vector3 diff;
        Vector3::DiffVector(&diff, &triV, &otherV);
        double dot;
        Vector3::DotProduct(&diff, tri.Normal, &dot);
        if(dot<=0){
            tri.Normal->Negate();
        }
    }

    /*
    std::vector<bool> mark;
    for(int i=0;i<cHull->size();i++){
        mark.push_back(false);
    }
    mark[0]=true;
    std::list<int> bfsList;
    bfsList.push_back(0);
    while(!bfsList.empty()){
        int triIndex = bfsList.pop_front();
        Triangle tri = cHull->at(triIndex);
        Tetrahedron tet = delcx->DeluanayTrigs[tri.ReverseLink1];
    }
    */

    return cHull;
}

bool inside(std::vector<Triangle> * cHull, std::vector<Vertex> &vertlist, Vector3 v){
    int positive =0;
    for(int i = 0; i< cHull->size();i++){
        Triangle t = cHull->at(i);
        if(t.Corners[1]==-1){
            continue;
        }
        Vertex a = vertlist[t.Corners[1]];
        Vector3 aVec = a.getCoordVector();
        Vector3 diff;
        Vector3::DiffVector(&diff, &aVec, &v);
        double dot;
        Vector3::DotProduct(&diff, t.Normal, &dot);
        if(dot>=0)positive++;
    }
//    std::cout << positive << " (" << cHull->size() << ")" << std::endl;
    if(positive==(cHull->size())||positive==0)return true;
    return false;
}

bool intersect(Triangle tri, std::vector<Vertex> &vertlist, Vector3* pt1, Vector3* pt2, Vector3* I){
    Vector3 dir;
    Vector3::DiffVector(&dir, pt2, pt1);
    Vector3 v0 = vertlist[tri.Corners[1]].getCoordVector();
    Vector3 v1 = vertlist[tri.Corners[2]].getCoordVector();
    Vector3 v2 = vertlist[tri.Corners[3]].getCoordVector();
    Vector3 u, v, n;
    Vector3::DiffVector(&u, &v1, &v0);
    Vector3::DiffVector(&v, &v2, &v0);
    Vector3::CrossProduct(&n, &u, &v);

    Vector3 w0;
    Vector3::DiffVector(&w0, pt1, &v0);

    double a, b;
    Vector3::DotProduct(&n, &w0 ,&a);
    a = -a;
    Vector3::DotProduct(&n, &dir, &b);

    if (fabs(b) < 0.000001) {       // segment is parallel to triangle plane
        if (a == 0)                 // segment lies in triangle plane
           return false;
        else return false;              // segment disjoint from plane
    }

    bool rayIntersect = true;
    // get intersect point of ray with triangle plane
    double r = a / b;
    if (r < 0.0)                  // ray goes away from triangle
       rayIntersect = false;                  // => no intersect
    if (r > 1.0)
       rayIntersect = false;

    Vector3 inc;
    Vector3::Scale(&inc, &dir, r);
    Vector3::Sum(I, pt1, &inc);

    if(!rayIntersect) return false;

    // is I inside T?
    double uu, uv, vv, wu, wv, D;
    Vector3::DotProduct(&u, &u, &uu);
    Vector3::DotProduct(&u, &v, &uv);
    Vector3::DotProduct(&v, &v, &vv);
    Vector3 w;
    Vector3::DiffVector(&w, I, &v0);
    Vector3::DotProduct(&w, &u, &wu);
    Vector3::DotProduct(&w, &v, &wv);
    D = uv * uv - uu * vv;

    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)        // I is outside T
        return false;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
        return false;

    return true;                   // I is in T
}

bool intersect(std::vector<Triangle> * cHull, std::vector<Vertex> &vertlist, Vector3* pt1, Vector3* pt2, Vector3* I){
    for(int i = 0; i< cHull->size();i++){
        Triangle tri = cHull->at(i);
        if(intersect(tri, vertlist, pt1, pt2, I)){
            return true;
        }
    }
    return false;
}

int getCommonTri(Tetrahedron* t1, Tetrahedron* t2){
    for(int i=1;i<5;i++){
        for(int j=1;j<5;j++){
            if((t1->TetLink[i]!=-1) && (t1->TetLink[i] == t2->TetLink[j])){
                return t1->TetLink[i];
            }
        }
    }
    return -1;
}

int getCHullTri(Tetrahedron* t, DeluanayComplex* delCplx){
    for(int i=1;i<5;i++){
        if((t->TetLink[i]!=-1) ){
            Triangle tri = delCplx->DeluanayTrigs[t->TetLink[i]];
            if(tri.Hull){
                return t->TetLink[i];
            }
        }
    }
    return -1;
}

PowerDiagram::PowerDiagram(DeluanayComplex* delCplx, std::vector<Vertex> &vertlist,
                           double min[], double max[]) : vertList(vertlist)
{
    this->delCplx = delCplx;
    for(int i=0;i<3;i++){
        this->min[i] = min[i];
        this->max[i] = max[i];
    }

    std::vector<Triangle> cHull = getConvexHull(delCplx, vertlist);

    int index = 0;
    std::vector <int> map;
    map.push_back(-1);
    for(int i=1; i<delCplx-> DeluanayTet.size(); i++){
        Tetrahedron tet = delCplx->DeluanayTet[i];
        PowerVertex vert;
        vert.tetIndex = i;
        if(tet.Status){
            vert.center = tet.powerVert(vertlist).getCoordVector();
            Vertex ball1 = vertlist[tet.Corners[1]];
            Vector3 ballCenter = ball1.getCoordVector();
            Vector3 diff;
            double dsq;
            Vector3::DiffVector(&diff, &ballCenter, &vert.center);
            Vector3::DotProduct(&diff, &diff, &dsq);
            vert.powerDistance = dsq - (ball1.Radius*ball1.Radius);

            if(inside(vert.center, min, max, 0.01) && inside (&cHull, vertList, vert.center)){
                vert.inside = true;
            } else {
                vert.inside = false;
            }
            vertices.push_back(vert);
            //vertices[index].setCenter(vert.center);
            map.push_back(index);
            index++;
        } else {
            map.push_back(-1);
        }
    }

    Metric metric(vertlist, vertlist.size());

    int edgeIndex = 0;
    for(int i=0;i<vertices.size();i++){
        PowerVertex* vert = &vertices[i];
        Tetrahedron tet = delCplx->DeluanayTet[vert->tetIndex];
        for(int j = 1 ; j < 5;j++){
            int x = map[tet.Neighbours[j]];
            if(x > i){
                Tetrahedron nb = delCplx->DeluanayTet[tet.Neighbours[j]];
                int commonTri = getCommonTri(&tet, &nb);
                PowerEdge edge;
                edge.v1 = i;
                edge.v2 = x;
                edge.triIndex = commonTri;
                if(vertices[i].inside && vertices[x].inside){
                    edge.edgeType = INSIDE;
                } else if (vertices[i].inside ^ vertices[x].inside){
                    edge.edgeType = INTERSECTING;
                    if(!vertices[i].inside){
                        edge.v1 = x;
                        edge.v2 = i;
                    }
                    intersect(&cHull, vertlist, &(vertices[i].center), &(vertices[x].center), &edge.intersect);
                } else{
                    edge.edgeType = OUTSIDE;
                }

                Triangle tri = delCplx->DeluanayTrigs[commonTri];
                bool intersects = intersect(tri, vertlist, &(vertices[i].center),
                          &(vertices[x].center), &edge.triIntersect);
                Vertex t1 = vertlist[tri.Corners[1]];
                Vector3 tv1 = t1.getCoordVector();
                Vector3 diff1;
                double d1sq;
                Vector3::DiffVector(&diff1, &tv1, &edge.triIntersect);
                Vector3::DotProduct(&diff1, &diff1, &d1sq);
                double pd1;
                pd1 = d1sq - (t1.Radius*t1.Radius);
                if(intersects){
                    edge.setLeastPowerDistance(pd1, &vertices, true);
                } else {
                    double leastPD = (vertices[i].powerDistance < vertices[x].powerDistance) ?
                                vertices[i].powerDistance : vertices[x].powerDistance;
                    edge.setLeastPowerDistance(leastPD, &vertices, true);
                }
                edges.push_back(edge);
                vert -> neigbours.push_back(edgeIndex);
                edgeIndex++;
            } else if (x == -1){
                PowerEdge edge;
                edge.v1 = i;
                edge.v2 = x;
                edge.triIndex = getCHullTri(&tet, delCplx);
                if(edge.triIndex==-1)
                    continue;
                if(vert->inside){
                    edge.edgeType = INTERSECTING | INFINITE;
                } else {
                    edge.edgeType = OUTSIDE | INFINITE;
                }

                Triangle tri = delCplx->DeluanayTrigs[edge.triIndex];
                edge.triIntersect = metric.Center3(vertlist, tri.Corners[1], tri.Corners[2], tri.Corners[3]);
                edge.intersect = edge.triIntersect;

                Vertex ball1 = vertlist[tri.Corners[1]];
                Vector3 ballCenter = ball1.getCoordVector();
                Vector3 diff;
                double dsq;
                Vector3::DiffVector(&diff, &ballCenter, &edge.intersect);
                Vector3::DotProduct(&diff, &diff, &dsq);
                edge.setLeastPowerDistance(dsq - (ball1.Radius*ball1.Radius), &vertices, true);
//                edge.intersect = delCplx->DeluanayTrigs[edge.triIndex].getCentroid(vertlist);

                if(vert->inside){
                    edges.push_back(edge);
                    vert -> neigbours.push_back(edgeIndex);
                    edgeIndex++;
                }
            } else {
                int eIndx = getEdgeTo(x, i);
                if(eIndx!=-1){
                    vert-> neigbours.push_back(eIndx);
                } else {
                    std::cerr << "Edge somehow not inserted earlier!" << std::endl;
                }
            }
        }
    }

    constructGraph(true);
    singlePath = true;
    startVert = targetVert = -1;

/*
    for(int i=0;i<edges.size();i++){
        PowerEdge edge = edges[i];
        std::cout << i << " \tType: " << edge.edgeType << " \tPD: " << edge.leastPowerDistance
                  << " \t\tWeight: " << edge.weight << " -- " <<
                     (edge.leastPowerDistance<0 && edge.edgeType==INSIDE) << std::endl;
    }
*/
}

GLuint listID = -1, pathListID = -1, pathSpheresListID = -1;

void PowerDiagram::render(bool showPowerDiag, bool showPath){
    if(showPowerDiag && glIsList(listID) == GL_TRUE){
        glCallList(listID);
        glEnable(GL_COLOR_MATERIAL);
        if(startVert>=0 && startVert<vertices.size()){
            glPointSize(10);
            glColor3d(0,1,1);
            glBegin(GL_POINTS);
            PowerVertex v = vertices[startVert];
            glVertex3d(v.center.X, v.center.Y, v.center.Z);
            glEnd();
        }
        if(targetVert>=0 && targetVert<vertices.size()){
            glPointSize(10);
            glColor3d(1,0,1);
            glBegin(GL_POINTS);
            PowerVertex v = vertices[targetVert];
            glVertex3d(v.center.X, v.center.Y, v.center.Z);
            glEnd();
        }
    }
    if(showPath && glIsList(pathListID) == GL_TRUE){
        glCallList(pathListID);
    }
    if(showPath && glIsList(pathSpheresListID) == GL_TRUE){
//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glCallList(pathSpheresListID);
//        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
}

void PowerDiagram::makeDisplayList(bool complementSpacePD, bool onlyInsideVerts,
                                   bool pruneIsolatedVerts, bool intersectEdges){
    if(glIsList(listID) == GL_TRUE){
        glDeleteLists(listID, 1);
    }
    listID = glGenLists(1);
    glNewList(listID, GL_COMPILE);

    glEnable(GL_COLOR_MATERIAL);

    glLineWidth(1);
    glColor3d(1,1,0);
    glBegin(GL_LINES);

    std::vector<bool> renderedIncidentEdge;
    for(int i=0;i<vertices.size();i++){
        renderedIncidentEdge.push_back(false);
    }

    for(int i=0;i<edges.size();i++){
        Vector3 v1 = vertices[edges[i].v1].center;
        switch(edges[i].edgeType){
            case INTERSECTING:
                glColor3d(0.1, 0.5, 1);
                break;
            case OUTSIDE:
                glColor3d(1, 0, 1);
                break;
            case INTERSECTING | INFINITE:
            case OUTSIDE | INFINITE:
                glColor3d(0, 0, 1);
                break;
            default:
                glColor3d(1, 1, 0);
        }
        Triangle tri = delCplx -> DeluanayTrigs[edges[i].triIndex];

        if((edges[i].edgeType & INFINITE)==0){
            // Finite edges

            Vector3 v2 = vertices[edges[i].v2].center;
    //        Tetrahedron t1 = delCplx->DeluanayTet[vertices[edges[i].v1].tetIndex];
    //        Tetrahedron t2 = delCplx->DeluanayTet[vertices[edges[i].v2].tetIndex];
    //        if(!complementSpacePD || (t1.AlphaStatus<1 && t2.AlphaStatus<1)){

            if(!onlyInsideVerts || edges[i].edgeType == INSIDE){
                if(!complementSpacePD || tri.AlphaStatus < 1){
                    glVertex3d(v1.X, v1.Y, v1.Z);
                    glVertex3d(v2.X, v2.Y, v2.Z);
                    renderedIncidentEdge[edges[i].v1] = true;
                    renderedIncidentEdge[edges[i].v2] = true;
//                    std::cout << "Edge: " << i << " \t" << edges[i].leastPowerDistance << std::endl;
                }
            }

            if(intersectEdges){
                if(onlyInsideVerts && edges[i].edgeType == INTERSECTING){
                    if(!complementSpacePD || tri.AlphaStatus < 1){
                        PowerVertex pv1 = vertices[edges[i].v1];
                        Vector3 I = edges[i].intersect;
                        Vector3 v = pv1.center;
                        glVertex3d(v.X, v.Y, v.Z);
                        glVertex3d(I.X, I.Y, I.Z);
                        renderedIncidentEdge[edges[i].v1] = true;
//                        std::cout << "Edge: " << i << " \t" << edges[i].leastPowerDistance << std::endl;
                    }
                }
            }
        }else{
            // Infinite Edges
            if(intersectEdges){
                if(!complementSpacePD || tri.AlphaStatus < 1){
                    Vector3 I = edges[i].triIntersect;
                    glVertex3d(v1.X, v1.Y, v1.Z);
                    glVertex3d(I.X, I.Y, I.Z);
                    renderedIncidentEdge[edges[i].v1] = true;
//                    std::cout << "Edge: " << i << " \t" << edges[i].leastPowerDistance << std::endl;
                }
            }
        }
    }
    glEnd();

    glColor3d(1,0,0);
    glPointSize(6);
    glBegin(GL_POINTS);
    for(int i=0;i<vertices.size();i++){
        Vector3 v = vertices[i].center;
        Tetrahedron t = delCplx->DeluanayTet[vertices[i].tetIndex];
        if (vertices[i].inside) {
            glColor3d(1,0,0);
        } else if (t.Hull) {
            glColor3d(0,0,1);
        } else {
            glColor3d(0,1,1);
        }
        if(!onlyInsideVerts || vertices[i].inside){
            if(!complementSpacePD || t.AlphaStatus<1){
//                std::cout << "Vertex: " << i << " \t" << vertices[i].powerDistance << std::endl;
                if(!pruneIsolatedVerts || renderedIncidentEdge[i]){
                    glVertex3d(v.X, v.Y, v.Z);
                }
            }
        }
    }
    glEnd();
    glDisable(GL_COLOR_MATERIAL);
    glEndList();
}

void convertToByte(uint i, GLubyte* bytes){
    i = i+1;
    bytes[0] = i & 0xff;
    bytes[1] = (i>>8) & 0xff;
    bytes[2] = (i>>16) & 0xff;
//    std::cout << i << "  " << (uint)bytes[0] << "  " << (uint)bytes[1] << "  " << (uint)bytes[2] << "  " << std::endl;
}

uint convertToInt(GLubyte* bytes){
    uint r = bytes[0];
    uint g = bytes[1];
    uint b = bytes[2];
    uint i = r + (g<<8) + (b<<16);
    return i;
}

void PowerDiagram::drawNodesForPicking() {
    glDisable(GL_DITHER);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);

    glPointSize(10);
    glBegin(GL_POINTS);
    GLubyte color[3];
    for(uint i=0;i<vertices.size();i++){
        Vector3 v = vertices[i].center;
        Tetrahedron t = delCplx->DeluanayTet[vertices[i].tetIndex];
        if(vertices[i].inside && t.AlphaStatus<1){
            convertToByte(i, color);
            glColor3ub(color[0], color[1], color[2]);
            glVertex3d(v.X, v.Y, v.Z);
        }
    }
    glEnd();
}

int PowerDiagram::processPick(int cursorX, int cursorY)
{
    GLint viewport[4];
    GLubyte pixel[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    glReadPixels(cursorX,viewport[3]-cursorY,1,1,GL_RGB ,GL_UNSIGNED_BYTE,(void *)pixel);
    int index = convertToInt(pixel);
    return index-1;
}

int PowerDiagram::getEdgeTo(int v1, int v2){
    PowerVertex pv1 = vertices[v1];
    for(int i=0;i<pv1.neigbours.size();i++){
        if(edges.at(pv1.neigbours[i]).isIncident(v2)){
            return pv1.neigbours[i];
        }
    }
    return -1;
}

void PowerDiagram::writeGraph(bool considerAlpha, const char* filename){
//    currentGraph.writeGraphCRD(filename);
    currentGraph.writeGraphOFF(filename);
    std::string cmd = "python /home/tbmasood/Downloads/common_modules/tritools.py off_to_vtp ";
    system(cmd.append(filename).data());
}

void PowerDiagram::savePathCRD(const char* filename){
    if(singlePath && !pathNodes.empty()){
        currentGraph.writePathCRD(filename, &pathNodes);
//        currentGraph.writePathOFF(filename, &pathNodes);
//        std::string cmd = "python /home/tbmasood/Downloads/common_modules/tritools.py off_to_vtp ";
//        system(cmd.append(filename).data());
    } else if(!singlePath && !pathsNodes.empty()){
        currentGraph.writeAllPathsCRD(filename, &pathsNodes);
    }
}

void PowerDiagram::constructGraph(bool considerAlpha){
    currentGraph.clear();
    vertexMap.clear();
    int index = 0;
    for(uint i=0; i<vertices.size(); i++){
        PowerVertex* vert = &vertices[i];
        Tetrahedron t = delCplx->DeluanayTet[vert->tetIndex];
        if(vert->inside && (!considerAlpha || (t.AlphaStatus < 1))){
            GraphNode node;
            node.index = index;
            node.boundary = false;
            node.x = vert->center.X;
            node.y = vert->center.Y;
            node.z = vert->center.Z;
            node.radius = vert->powerDistance <= 0 ? 0 : sqrt(vert->powerDistance);
            node.pVert = vert;
            currentGraph.nodes.push_back(node);
            vertexMap.push_back(index);
            index++;
        }else{
            vertexMap.push_back(-1);
        }
    }
    int edgeIndex = 0;
    for(int i=0;i<edges.size();i++){
        PowerEdge * pEdge = &edges[i];
        Triangle tri = delCplx -> DeluanayTrigs[pEdge->triIndex];
        if(!considerAlpha || (tri.AlphaStatus<1)){
            if(pEdge->edgeType == INSIDE){
                GraphEdge edge;
                edge.index = edgeIndex;
                edge.pEdge = pEdge;
                edge.v1 = vertexMap[pEdge->v1];
                edge.v2 = vertexMap[pEdge->v2];

                currentGraph.nodes[edge.v1].edges.push_back(edgeIndex);
                currentGraph.nodes[edge.v2].edges.push_back(edgeIndex);

                edge.weight = pEdge->weight;
                currentGraph.edges.push_back(edge);
                edgeIndex++;
            } else if(pEdge->edgeType == INTERSECTING){
                GraphEdge edge;
                edge.index = edgeIndex;
                edge.pEdge = pEdge;
                edge.v1 = vertexMap[pEdge->v1];

                GraphNode node;
                node.index = index;
                node.boundary = true;
                node.x = pEdge->intersect.X;
                node.y = pEdge->intersect.Y;
                node.z = pEdge->intersect.Z;
                // check -- this may not be correct
                node.radius = pEdge->leastPowerDistance <= 0 ? 0 : sqrt(pEdge->leastPowerDistance);
                node.pVert = &vertices[pEdge->v2];
                node.edges.push_back(edgeIndex);
                currentGraph.nodes.push_back(node);
                edge.v2 = index;
                index++;

                currentGraph.nodes[edge.v1].edges.push_back(edgeIndex);

                edge.weight = pEdge->weight;
                currentGraph.edges.push_back(edge);
                edgeIndex++;
            } else if(pEdge->edgeType == (INTERSECTING | INFINITE)){
                GraphEdge edge;
                edge.index = edgeIndex;
                edge.pEdge = pEdge;
                edge.v1 = vertexMap[pEdge->v1];

                GraphNode node;
                node.index = index;
                node.boundary = true;
                node.x = pEdge->intersect.X;
                node.y = pEdge->intersect.Y;
                node.z = pEdge->intersect.Z;
                node.radius = pEdge->leastPowerDistance <= 0 ? 0 : sqrt(pEdge->leastPowerDistance);
                node.edges.push_back(edgeIndex);
                currentGraph.nodes.push_back(node);
                edge.v2 = index;
                index++;

                currentGraph.nodes[edge.v1].edges.push_back(edgeIndex);

                edge.weight = pEdge->weight;
                currentGraph.edges.push_back(edge);
                edgeIndex++;
            }
        }
    }
    currentGraph.initLemonGraph();
}

void drawPath(std::vector<GraphNode*> *pathNodes, std::vector<GraphEdge*> *pathEdges,
              GLUquadric* quad, bool useSelected){
//  glBegin(GL_LINES);
    for(int i=0;i<pathNodes->size()-1;i++){
        if(useSelected){
            if(pathEdges->at(i)->selected){
                glColor3d(1, 0, 1);
            }else{
                glColor3d(0, 1, 1);
            }
        }
        GraphNode* curr = pathNodes->at(i);
        GraphNode* next = pathNodes->at(i+1);
        Vector3 n(next->x - curr->x, next->y - curr->y, next->z - curr->z);
        double len;
        Vector3::DotProduct(&n, &n, &len);
        len = sqrt(len);
        n.Normalize();
        glPushMatrix();
        double div = sqrt(n.X * n.X + n.Z * n.Z);
        if(div==0){
            div = 0.0000001;
        }
        double m[] = {n.Z/div, 0, -n.X/div, 0,
                     -(n.X*n.Y)/div , div, -(n.Y*n.Z)/div, 0,
                      n.X, n.Y, n.Z, 0,
                     curr->x, curr->y, curr->z, 1};
        glMultMatrixd(m);
        gluCylinder(quad, 0.05, 0.05, len ,5, 5);
        glPopMatrix();
//      glVertex3d(curr->x, curr->y, curr->z);
//      glVertex3d(next->x, next->y, next->z);
    }
//  glEnd();

    glPointSize(8);
//        glBegin(GL_POINTS);
    for(int i=0;i<pathNodes->size();i++){
        if(i==0){
            glColor3d(0.8,0.2,0);
        } else if (i == pathNodes->size()-1){
            glColor3d(0,0.5,0);
        } else {
            glColor3d(0,0,1);
        }
        GraphNode* curr = pathNodes->at(i);
        glPushMatrix();
        glTranslated(curr->x, curr->y, curr->z);
        gluSphere(quad, 0.1, 6, 6);
        glPopMatrix();
//      glVertex3d(curr->x, curr->y, curr->z);
    }
//  glEnd();
}

void initPathList(std::vector<GraphNode*> *pathNodes){
    if(glIsList(pathSpheresListID) == GL_TRUE){
        glDeleteLists(pathSpheresListID, 1);
    }
    pathSpheresListID = glGenLists(1);
    glNewList(pathSpheresListID, GL_COMPILE);

    glEnable(GL_COLOR_MATERIAL);

//        glColor3d(1,0,1);

//        for(int i=0;i<pathNodes.size();i++){
//            GraphNode* curr = pathNodes[i];
//            glPushMatrix();
//            glTranslated(curr->x, curr->y, curr->z);
//            gluSphere(quad, sqrt(curr->pVert->powerDistance), 8, 8);
//            glPopMatrix();
//        }

    FILE *fp = fopen("pathskin","w");
    fprintf(fp,"%d\n",pathNodes->size());
    fprintf(fp,"#junk\n");
    for(int i=0;i<pathNodes->size();i++){
        GraphNode* curr = pathNodes->at(i);
        fprintf(fp,"%d %f %f %f %f\n",i+1,curr->x, curr->y, curr->z,
                curr->radius);
    }
    fclose(fp);
    system("./smesh pathskin -s skin.off -t skin.tet");

    double center[] = {0,0,0};
    SkinSurface skin(center, 1, 1);
    skin.Read("skin_lev0.off",0);
    skin.Process();
    skin.Draw(true, true);

    glDisable(GL_COLOR_MATERIAL);

    glEndList();
}

#include<set>
void initMultiplePathList(std::vector<std::vector<GraphNode*> > *pathsNodes){
    if(glIsList(pathSpheresListID) == GL_TRUE){
        glDeleteLists(pathSpheresListID, 1);
    }
    pathSpheresListID = glGenLists(1);
    glNewList(pathSpheresListID, GL_COMPILE);

    glEnable(GL_COLOR_MATERIAL);
    std::vector<GraphNode*> pathNodes;
    std::set<int> unique;

    for(int i=0;i<pathsNodes->size();i++){
        uint size = pathsNodes->at(i).size();
        for(uint j=0; j<size; j++){
            GraphNode* node = pathsNodes->at(i)[j];
            if(unique.find(node->index) == unique.end()){
                unique.insert(node->index);
                pathNodes.push_back(node);
            }
        }
    }

    FILE *fp = fopen("pathskin","w");
    fprintf(fp,"%d\n",pathNodes.size());
    fprintf(fp,"#junk\n");
    for(int i=0;i<pathNodes.size();i++){
        GraphNode* curr = pathNodes.at(i);
        fprintf(fp,"%d %f %f %f %f\n",i+1,curr->x, curr->y, curr->z,
                curr->radius);
    }
    fclose(fp);
    system("./smesh pathskin -s skin.off -t skin.tet");

    double center[] = {0,0,0};
    SkinSurface skin(center, 1, 1);
    skin.Read("skin_lev0.off",0);
    skin.Process();
    skin.Draw(true, true);

    glDisable(GL_COLOR_MATERIAL);

    glEndList();
}

bool PowerDiagram::findShortestPath(QVector<double>* X, QVector<double>* Y,
                                    double * length, double *minY, double *maxY){
    if(startVert<0 || targetVert<0)
        return false;
    int start = vertexMap[startVert];
    int target = vertexMap[targetVert];
    std::vector<GraphEdge*> pathEdges;
    double cost = currentGraph.runDijkstra(start, target, &pathNodes, &pathEdges);
    std::cout << cost << std::endl;
    if(!pathNodes.empty()){
        getPathWeights(&pathNodes, &pathEdges, X, Y, length, minY, maxY);

        if(glIsList(pathListID) == GL_TRUE){
            glDeleteLists(pathListID, 1);
        }
        pathListID = glGenLists(1);
        glNewList(pathListID, GL_COMPILE);

        glEnable(GL_COLOR_MATERIAL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glLineWidth(3);
        glColor3d(1,0,1);

        GLUquadric* quad = gluNewQuadric();
        gluQuadricNormals(quad, GLU_SMOOTH);
        gluQuadricOrientation(quad, GLU_OUTSIDE);

        drawPath(&pathNodes, &pathEdges, quad, false);

        glLineWidth(1);
        glPointSize(6);

        glDisable(GL_COLOR_MATERIAL);

        glEndList();

        initPathList(&pathNodes);
        singlePath = true;

        return true;
    }else{
        return false;
    }
}

int PowerDiagram::findShortestEscapePaths(int steps, bool repeated, int maxIter,
                                    std::vector<QVector<double> >* Xs, std::vector<QVector<double> >* Ys,
                                    std::vector<double> * lengths, std::vector<double> *minYs,
                                    std::vector<double> *maxYs){
    if(startVert<0)
        return false;
    int start = vertexMap[startVert];

    std::vector<std::vector<GraphEdge*> > pathsEdges;

    int shortest = -1;
    if(repeated){
        bool found = currentGraph.runDijkstraEscapeRepeated(start, maxIter, &pathsNodes, &pathsEdges);
        if(found){
           shortest = 0;
        }
    } else {
        shortest = currentGraph.runDijkstra(start, &pathsNodes, &pathsEdges);
    }

    for(int k=0;k < pathsEdges.size();k++){
        std::vector<GraphEdge*> pathEdges = pathsEdges[k];
        for(int j=0;j<pathEdges.size();j++){
            pathEdges[j]->selected = false;
        }
    }
    if(shortest!=-1){
        for(int j=0;j<pathsEdges[shortest].size();j++){
            pathsEdges[shortest][j]->selected = true;
        }
    }

    for(int k=0;k < pathsNodes.size();k++){
        QVector<double> X(steps), Y(steps);
        Xs->push_back(X);
        Ys->push_back(Y);
        lengths->push_back(0);
        minYs->push_back(0);
        maxYs->push_back(0);
    }
    if(!pathsNodes.empty()){

        if(glIsList(pathListID) == GL_TRUE){
            glDeleteLists(pathListID, 1);
        }
        pathListID = glGenLists(1);
        glNewList(pathListID, GL_COMPILE);

        glEnable(GL_COLOR_MATERIAL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        GLUquadric* quad = gluNewQuadric();
        gluQuadricNormals(quad, GLU_SMOOTH);
        gluQuadricOrientation(quad, GLU_OUTSIDE);

        for(int k=0;k < pathsNodes.size();k++){
            std::vector<GraphNode*> pathNodes = pathsNodes[k];
            std::vector<GraphEdge*> pathEdges = pathsEdges[k];
            if(!pathNodes.empty()){
                if(k==shortest){
                    glColor3d(1,0,1);
                }else{
                    glColor3d(0,1,1);
                }
                drawPath(&pathNodes, &pathEdges, quad, true);
            }
            getPathWeights(&pathNodes, &pathEdges, &(Xs->at(k)), &(Ys->at(k)), &(lengths->at(k)),
                           &(minYs->at(k)), &(maxYs->at(k)));
        }

        glLineWidth(1);
        glPointSize(6);

        glDisable(GL_COLOR_MATERIAL);

        glEndList();

        if(repeated){
            initMultiplePathList(&pathsNodes);
        }else if(glIsList(pathSpheresListID) == GL_TRUE){
            glDeleteLists(pathSpheresListID, 1);
        }
        singlePath = false;
    }
    return shortest;
}

bool PowerDiagram::findShortestEscapePath(QVector<double>* X, QVector<double>* Y,
                                          double * length, double *minY, double *maxY){
    if(startVert<0)
        return false;
    int start = vertexMap[startVert];

    std::vector<GraphEdge*> pathEdges;
    bool pathFound = currentGraph.runDijkstraEscape(start, &pathNodes, &pathEdges);
    if(pathFound && !pathNodes.empty()){
        getPathWeights(&pathNodes, &pathEdges, X, Y, length, minY, maxY);

        GLUquadric* quad = gluNewQuadric();
        gluQuadricNormals(quad, GLU_SMOOTH);
        gluQuadricOrientation(quad, GLU_OUTSIDE);

        if(glIsList(pathListID) == GL_TRUE){
            glDeleteLists(pathListID, 1);
        }
        pathListID = glGenLists(1);
        glNewList(pathListID, GL_COMPILE);

        glEnable(GL_COLOR_MATERIAL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glLineWidth(3);
        glColor3d(1,0,1);

        drawPath(&pathNodes, &pathEdges, quad, false);

        glLineWidth(1);
        glPointSize(6);

        glDisable(GL_COLOR_MATERIAL);

        glEndList();

        initPathList(&pathNodes);
        singlePath = true;
    }
    return pathFound;
}

void PowerDiagram::getPathWeights(std::vector<GraphNode*> *pathNodes,  std::vector<GraphEdge*> *pathEdges,
                            QVector<double>* X, QVector<double>* Y, double * length, double *minY, double *maxY)
{
    int steps = X->size();
    double totalLength = 0;
    for(uint i = 0; i< pathEdges->size();i++){
        GraphEdge* gEdge = pathEdges->at(i);
        totalLength += gEdge->length(&currentGraph.nodes);
    }
    *length = totalLength;
    assert (pathEdges->size()+1 == pathNodes->size());
    assert (totalLength > 0 && steps > 0);
    double alreadyCovered = 0;
    double stepSize = totalLength/steps;
    int index = 0;
    for(uint i = 0; i < pathEdges->size();i++){
        GraphEdge* gEdge = pathEdges->at(i);
        GraphNode* n1 = pathNodes->at(i);
        Vector3 n1Vec(n1->x, n1->y, n1->z);
        Vector3 dir = gEdge->direction(&currentGraph.nodes);

        Triangle tri = delCplx->DeluanayTrigs[gEdge->pEdge->triIndex];
        Vertex t1 = vertList[tri.Corners[1]];
        Vector3 tv1 = t1.getCoordVector();

        double  edgeLen = gEdge->length(&currentGraph.nodes);
        double curr = stepSize - alreadyCovered;
        if(i == 0){
            curr = 0;
        }
        bool containedPoint = false;
        while(curr <= edgeLen){
            containedPoint = true;
            if(index == steps){
                break;  // for safety
            }
            Vector3 pt;
            Vector3 toAdd;
            if(n1->index == gEdge->v1){
                Vector3::Scale(&toAdd, &dir, curr);
                Vector3::Sum(&pt, &n1Vec, &toAdd);
            } else {
                Vector3::Scale(&toAdd, &dir, -curr);
                Vector3::Sum(&pt, &n1Vec, &toAdd);
            }
            Vector3 diff;
            double dsq;
            Vector3::DiffVector(&diff, &tv1, &pt);
            Vector3::DotProduct(&diff, &diff, &dsq);
            double pd1;
            pd1 = dsq - (t1.Radius*t1.Radius);
//            pd1 = pd1<0? 0 : sqrt(pd1);

            if(index == 0){
                *minY = *maxY = pd1;
            } else {
                *minY = (*minY < pd1)? *minY : pd1;
                *maxY = (*maxY > pd1)? *maxY : pd1;
            }
            (*X)[index] = index * stepSize;
            (*Y)[index] = pd1;
            index++;
            curr+=stepSize;
        }
        if(containedPoint){
            alreadyCovered = edgeLen - curr + stepSize;
        } else {
            alreadyCovered += edgeLen;
        }
    }
}

void PowerDiagram::printGraph(){
    std::cout << std::endl << "Printing Graph" << std::endl;
    std::cout << "Vertices:" << std::endl;
    for(int i=0;i<vertices.size();i++){
        PowerVertex vert = vertices[i];
        std::cout << vert.tetIndex << " --- ";
        for(int j=0; j<vert.neigbours.size(); j++){
            std::cout << vert.neigbours[j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Edges:" << std::endl;
    for(int i=0;i<edges.size();i++){
        PowerEdge edge = edges[i];
        std::cout << edge.v1 << " "<< edge.v2 << std::endl;
    }
}

PowerDiagram::~PowerDiagram()
{
    edges.clear();
    vertices.clear();
}

