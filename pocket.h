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

#ifndef POCKET_H
#define POCKET_H

#include<QTableWidgetItem>
#include <QTableWidget>
#include <QLineEdit>
#include <QtOpenGL>
#include <vector>
#include <spacefillmeasure.h>
#include <vertex.h>
#include <deluanaycomplex.h>
#include <disjointset.h>
#include <lightmaterial.h>
#include <skinsurface.h>
#include <volume.h>
#include <alphacomplex.h>

//#include <GL/glut.h>

class Pocket
{
    private:
        int Rank;
        double PI;
        void DrawTet(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int l,bool mouth);
        void InitMaterial(int i);
        SkinSurface *skin;
        double pcenter[6];
        double pscale;
        std::vector<int> pocPersistence;

    public:
        int npockets;
        double TotVolume;
        double TotSurfArea;
        std::vector<std::vector<int> > AllPockets;
        int PocketNVertex;
        std::vector<int> PocketVertex;
        std::vector<double> SurfPocket;
        std::vector<double> VolPocket;
        std::vector<int> vertices;
        std::vector<double> Radius2;
        SpaceFillMeasure *spm;
        DisJointSet *unionFind;
        Volume *vol;

        Pocket(std::vector<Vertex> & vertexList,double center[],double scale,int trcount,int ecount);
        ~Pocket();

        bool InComplex(DeluanayComplex *delcx,std::vector<Vertex> &vertexList,unsigned int ftype, int rank,int i);

        void FindPockets(DeluanayComplex *delcx, std::vector <int> &sortedTet);
        std::vector< std::vector<int> > GetPockets();
        void PocketProperties(DeluanayComplex *delcx, std::vector<Vertex> &vertexList,int rank);
        void MeasurePockets(DeluanayComplex *delcx, std::vector<Vertex> & vertexList);
        void FillTable(QTableWidget *table,std::vector<int> & PocketNMouths,QLineEdit *pocVol,QLineEdit *pocSurf);
        void Render(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,std::vector<int> & PocketNMouths,bool allPockets,bool onlyPockets,bool onlyVoids,bool skinSurface,int persistence,int rank,bool smoothShading,bool skinWireFrame,bool mouth,bool pocketWireFrame);
        void RenderSingle(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,std::vector<int> & PocketNMouths,bool allPockets,bool onlyPockets,bool onlyVoids,bool skinSurface,int persistence,int rank,int pocIndex,bool smoothShading,bool skinWireFrame,bool mouth,bool pocketWireFrame);
        void RenderTwo(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,std::vector<int> & PocketNMouths,bool allPockets,bool onlyPockets,bool onlyVoids,bool skinSurface,int persistence,int rank,int pocIndex,bool smoothShading,bool skinWireFrame,bool mouth,bool pocketWireFrame);
};

#endif // POCKET_H
