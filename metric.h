/***************************************************************************
 *   Copyright (C) 2010 by raghavendra,,,                                  *
 *   raghavendra@incognito                                                 *
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

#ifndef METRIC_H
#define METRIC_H


#include <vector3.h>

#include <sos.h>
#include <vector>
#include <vertex.h>
#include <gmp.h>
#include <size.h>

typedef double Matrix[4][4];

class Metric
{
private:
    double PI;
    std::vector<double> Radius2;
    Size *size;
    Sos *sos;
public:
    Metric(std::vector<Vertex> &vetexList, int vcount);
    ~Metric();

    double Det2(Matrix A);
    double Det3(Matrix A);
    double Det4(Matrix A);

    void Angle(double ang, Vector3* T, Vector3* u, Vector3* v);

    int Ccw(std::vector<Vertex> vertexList,int i,int j,int k,int l);

    Vector3 Center2(std::vector<Vertex> &vertexList, int i, int j);
    Vector3 Center3(std::vector<Vertex> &vertexList, int i, int j, int k);
    Vector3 Center4(std::vector<Vertex> &vertexList, int i, int j, int k, int l);

    Vector3 TriangleDual(std::vector<Vertex> &vertexList, int i, int j, int k);

    double AngleDihedral(Vector3 s, Vector3 t, Vector3 u, Vector3 v);
    double AngleSolid(std::vector<Vertex> &vertexList, int i, int j, int k, int l);

    double DiskRadius(std::vector<Vertex> &vertexList, int i, int j);
    double DiskLength(std::vector<Vertex> &vertexList, int i, int j);
    double DiskArea(std::vector<Vertex> &vertexList, int i, int j);

    double SegmentHeight(std::vector<Vertex> &vertexList, int i, int j, int k);
    double SegmentAngle(std::vector<Vertex> &vertexList, int i, int j, int k);
    double SegmentLength(std::vector<Vertex> &vertexList, int i, int j, int k);
    double SegmentArea(std::vector<Vertex> &vertexList, int i, int j, int k);

    double Segment2Angle(std::vector<Vertex> &vertexList, int i, int j, int k, int l);
    double Segment2Length(std::vector<Vertex> &vertexList, int i, int j, int k, int l);
    double Segment2Area(std::vector<Vertex> &vertexList, int i, int j, int k, int l);

    double CapHeight(std::vector<Vertex> &vertexList, int i, int j);
    double CapArea(std::vector<Vertex> &vertexList, int i, int j);
    double CapVolume(std::vector<Vertex> &vertexList, int i, int j);

    double Cap2Volume(std::vector<Vertex> &vertexList, int i, int j, int k);
    double Cap2Area(std::vector<Vertex> &vertexList, int i, int j, int k);

    double Cap3Volume(std::vector<Vertex> &vertexList, int i, int j, int k, int l);
    double Cap3Area(std::vector<Vertex> &vertexList, int i, int j, int k, int l);

    double SectorArea(std::vector<Vertex> &vertexList, int i, int j, int k, int l);
    double SectorVolume(std::vector<Vertex> &vertexList, int i, int j, int k, int l);

    double WedgeLength(std::vector<Vertex> &vertexList, int i, int j, int k, int l);
    double WedgeArea(std::vector<Vertex> &vertexList, int i, int j, int k, int l);
    double WedgeVolume(std::vector<Vertex> &vertexList, int i, int j, int k, int l);

    double PawnLength(std::vector<Vertex> &vertexList, int i, int j, int k);
    double PawnArea(std::vector<Vertex> &vertexList, int i, int j, int k);
    double PawnVolume(std::vector<Vertex> &vertexList, int i, int j, int k);

    double BallVolume(std::vector<Vertex> &vertexList, int idx);
    double BallArea(int idx);

    double Ball2Length(std::vector<Vertex> &vertexList, int i, int j);
    double Ball2Volume(std::vector<Vertex> &vertexList, int i, int j);
    double Ball2Area(std::vector<Vertex> &vertexList, int i, int j);

    double Ball3Length(std::vector<Vertex> &vertexList, int i, int j, int k);
    double Ball3Volume(std::vector<Vertex> &vertexList, int i, int j, int k);
    double Ball3Area(std::vector<Vertex> &vertexList, int i, int j, int k);

    double Ball4Length(std::vector<Vertex> &vertexList, int i, int j, int k, int l);
    double Ball4Volume(std::vector<Vertex> &vertexList, int i, int j, int k, int l);
    double Ball4Area(std::vector<Vertex> &vertexList, int i, int j, int k, int l);

    double TetraVolume(double a[], double b[], double c[], double d[]);
};

#endif // METRIC_H
