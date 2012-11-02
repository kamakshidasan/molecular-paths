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

#ifndef VOLUME_H
#define VOLUME_H

#include <vector>
#include <QLineEdit>
#include <cmath>
#include <spacefillmeasure.h>
#include <deluanaycomplex.h>
#include <trigvol.h>
#include <linktrig.h>
#include <metric.h>

class Volume
{
    private:
        double PI;
//        AlfMasterNode* mlNode;

    public:
        double TotVolume;
        double TotSurfaceArea;
        double Ssolv;
        double Vsolv;
        std::vector<double> IndVolume;
        std::vector<double> IndSurfArea;
        std::vector<TrigVol> DerivVol;
        std::vector<TrigVol> DerivSurf;

        std::vector<double> trig_eps;//(ntrig_max)
        std::vector<TrigVol> trig_sh;//(3,ntrig_max)
        std::vector<TrigVol> trig_dual1;//(3,ntrig_max)
        std::vector<TrigVol> trig_dual2;//(3,ntrig_max)
        std::vector<TrigVol> trig_ang;

        std::vector<int> nlink_trig;
        std::vector<LinkTrig> link_trig;

        std::vector<double> distpair;//(nedge_max)
        std::vector<double> distpair2;//(nedge_max)

        std::vector<double> edgeCoef;

        std::vector<double> Radius2;

        std::vector<double> CoefAsp;

        SpaceFillMeasure *spm;

        Metric *met;

        int test;

        Volume();
        Volume(std::vector<Vertex> &vertexList,int vcount,int trcount,int ecount);
        ~Volume();

        void VertexProperties(std::vector<Vertex> &vertexList,int idx,FILE *fp1);
        void EdgeProperties(std::vector<Vertex> &vertexList,int i,int j,FILE *fp1);
        void TriangleProperties(std::vector<Vertex> &vertexList,int i,int j,int k,FILE *fp1);
        void TetProperties(std::vector<Vertex> &vertexList,int i,int j,int k,int l,FILE *fp1);
        void FindTotal(std::vector<Vertex> &vertexList,double *TotVol, double *TotSurf);

        void DoTetraVertex(std::vector<Vertex> &verteList,int i,int j,int k,int l,int m,double *surf,double *vol);
        void DoTetraEdge(std::vector<Vertex> &vertexList,int i,int j,int k, int l,int m,double *surf,double *vol);
        void DoTetraTriangle(std::vector<Vertex> &vertexList ,int i,int j,int k,int l,int m,double *surf,double *vol);
        double DoTetraVolume(std::vector<Vertex> &vertexList,int i,int j,int k,int l);

        void MeasureVolume(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,QLineEdit *totVol,QLineEdit *totSurf,int option);
        //void MeasureVolume(AlphaComplex *alcx,std::vector<Vertex> & vertexList,QLineEdit *totVol,QLineEdit *totSurf,int option,int Rank);
        //void FindVolume(AlphaComplex *alcx,std::vector<Vertex> & vertexList,QLineEdit *totVol,QLineEdit *totSurf,int option,int Rank);
};

#endif // VOLUME_H
