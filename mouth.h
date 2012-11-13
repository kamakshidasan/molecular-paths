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

#ifndef MOUTH_H
#define MOUTH_H

#include<QTableWidgetItem>
#include <QTableWidget>
#include <QLineEdit>
#include <GL/glew.h>
#include <cmath>
#include <vector>
#include <disjointset.h>
#include <vector3.h>
#include <deluanaycomplex.h>
#include <lightmaterial.h>

class Mouth
{
    private:
        double PI;
        void TwoCircle(double ra, double rb, double dist, double *surf);
        void TriangleAng(double a[], double b[], double c[], double *ang1, double *ang2, double *ang3, double *dist1, double *dist2, double *dist3);
        void InitMaterial(int i);
        void DrawTrig(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int l);

    public:
        int nmouths;
        int MouthNVertex;
        double TotSurfArea;

        std::vector<int> MouthPocket;
        DisJointSet *unionFind;
        std::vector<std::vector<int> > AllMouths;
        std::vector<double> SurfMouth;
        std::vector<int> MouthVertex;
        std::vector<int> vertices;

        Mouth();
        ~Mouth();

        int mouthPocket(int i, DeluanayComplex *delcx);
        void CleanUp(DeluanayComplex *delcx);
        void FindMouths(std::vector<std::vector<int> > &pockets, DeluanayComplex *delcx);
        void MeasureMouths(DeluanayComplex *delcx, std::vector<Vertex> & vertexList);
        void FillTable(QTableWidget *table,QLineEdit *mouSurf);
        void Render(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int persistence);
        void RenderMouthOfPocket(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int persistence,int pocIndex);
};

#endif // MOUTH_H
