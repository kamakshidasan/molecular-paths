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

#include "mouth.h"

Mouth::Mouth()
{
        PI = 4 * atan(1);
        unionFind = new DisJointSet();
}

Mouth::~Mouth()
{
}

/*!
    \fn Mouth::TwoCircle(double ra, double rb, double dist, double *surf)
 */
void Mouth::TwoCircle(double ra, double rb, double dist, double *surf)
{
        double ra2, rb2, dist2;
        double val, eps, ang1, ang2, s;

        eps = 1.0e-8;

        ra2 = ra * ra;
        rb2 = rb * rb;
        dist2 = dist * dist;

        val = (dist2 + ra2 - rb2) / (2 * dist * ra);
        if (val >= 1) val = 1.0 - eps;
        if (val <= -1) val = -1.0 + eps;
        ang1 = acos(val);

        val = (dist2 + rb2 - ra2) / (2 * dist * rb);
        if (val >= 1) val = 1.0 - eps;
        if (val <= -1) val = -1.0 + eps;
        ang2 = acos(val);

        s = (-dist + ra + rb) * (dist + ra - rb) * (dist - ra + rb) * (dist + ra + rb);
        s = 0.5 * sqrt(s);

        *surf = ra2 * ang1 + rb2 * ang2 - s;
}

/*!
    \fn Mouth::TriangleAng(double a[], double b[], double c[], double *ang1, double *ang2, double *ang3, double *dist1, double *dist2, double *dist3)
 */
void Mouth::TriangleAng(double a[], double b[], double c[], double *ang1, double *ang2, double *ang3, double *dist1, double *dist2, double *dist3)
{
        double x, val, eps;

        Vector3 *u1 = new Vector3();
        Vector3 *u2 = new Vector3();
        Vector3 *u3 = new Vector3();

        Vector3 *va = new Vector3(a[1], a[2], a[3]);
        Vector3 *vb = new Vector3(b[1], b[2], b[3]);
        Vector3 *vc = new Vector3(c[1], c[2], c[3]);

        eps = 1.0e-8;

        Vector3::DiffVector(u1,vb,va);
        Vector3::DotProduct(u1, u1,dist1);
        *dist1 = sqrt(*dist1);

        Vector3::DiffVector(u2,vc,va);
        Vector3::DotProduct(u2, u2,dist2);
        *dist2 = sqrt(*dist2);

        Vector3::DiffVector(u3,vc,vb);
        Vector3::DotProduct(u3, u3,dist3);
        *dist3 = sqrt(*dist3);

        Vector3::DotProduct(u1, u2,&val);
        x = val / (*dist1 * *dist2);
        if (x >= 1) x = 1.0 - eps;
        if (x <= -1) x = -1.0 + eps;
        *ang1 = acos(x) / (2 * PI);

        Vector3::DotProduct(u1, u3,&val);
        x = val / (*dist1 * *dist3);
        if (x >= 1) x = 1.0 - eps;
        if (x <= -1) x = -1.0 + eps;
        *ang2 = acos(x) / (2 * PI);

        *ang3 = 0.5 - *ang1 - *ang2;
}

/*!
    \fn Mouth::mouthPocket(int i, DeluanayComplex *delcx)
 */
int Mouth::mouthPocket(int i, DeluanayComplex *delcx)
{
        MouthPocket.push_back(delcx->DeluanayTrigs[AllMouths[i][0]].PocIndex);
        return (delcx->DeluanayTrigs[AllMouths[i][0]].PocIndex);
}

/*!
    \fn Mouth::CleanUp(DeluanayComplex delcx)
 */
void Mouth::CleanUp(DeluanayComplex *delcx)
{
        for (uint i = 0; i < AllMouths.size(); i++)
        {
                for (uint j = 0; j < AllMouths[i].size(); j++)
                {
                        int index = AllMouths[i][j];
                        delcx->DeluanayTrigs[index].ufKey = -1;
                        delcx->DeluanayTrigs[index].PocIndex = -1;
                }
        }
}

/*!
    \fn Mouth::FindMouths(std::vector<std::vector<int> > &pockets, DeluanayComplex *delcx)
 */
void Mouth::FindMouths(std::vector<std::vector<int> > &pockets, DeluanayComplex *delcx)
{
        uint i, j;
        int k, pindex, itrig, jtrig, neighbour;
        int iset, jset;
        int pairi_1, pairi_2, pairi_3;
        int pairj_1, pairj_2, pairj_3;
        unionFind->Clear();
        for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
        {
                delcx->DeluanayTrigs[i].IsMouth = false;
        }
        for (i = 0; i < pockets.size(); i++)
        {
                for (j = 0; j < pockets[i].size(); j++)
                {
                        pindex = pockets[i][j];
                        for (k = 1; k <= 4; k++)
                        {
                                neighbour = delcx->DeluanayTet[pindex].Neighbours[k];
                                itrig = delcx->DeluanayTet[pindex].TetLink[k];

                                if (neighbour == 0)                              //on hull
                                {
                                        if (delcx->DeluanayTrigs[itrig].AlphaStatus == 0)
                                        {
                                                delcx->DeluanayTrigs[itrig].ufKey = unionFind->ElementCount();
                                                unionFind->Add(itrig);
                                                delcx->DeluanayTrigs[itrig].PocIndex = i;
                                        }
                                }
                                else
                                {
                                        if ((delcx->DeluanayTet[pindex].PocIndex != delcx->DeluanayTet[neighbour].PocIndex) && delcx->DeluanayTrigs[itrig].AlphaStatus == 0)
                                        {
                                                delcx->DeluanayTrigs[itrig].ufKey = unionFind->ElementCount();
                                                unionFind->Add(itrig);
                                                delcx->DeluanayTrigs[itrig].PocIndex = i;
                                        }
                                }
                        }
                }
        }
        for (i = 0; i < unionFind->ElementCount(); i++)
        {
                itrig = unionFind->GetElementAt(i);
                pairi_1 = delcx->DeluanayTrigs[itrig].TrigLink[1];
                pairi_2 = delcx->DeluanayTrigs[itrig].TrigLink[2];
                pairi_3 = delcx->DeluanayTrigs[itrig].TrigLink[3];

                for (j = i + 1; j < unionFind->ElementCount(); j++)
                {
                        iset = unionFind->FindSet(i);        //important as iset can change
                        jset = unionFind->FindSet(j);
                        if (iset == jset) continue;

                        jtrig = unionFind->GetElementAt(j);

                        pairj_1 = delcx->DeluanayTrigs[jtrig].TrigLink[1];
                        pairj_2 = delcx->DeluanayTrigs[jtrig].TrigLink[2];
                        pairj_3 = delcx->DeluanayTrigs[jtrig].TrigLink[3];

                        if ((pairi_1 == pairj_1) || (pairi_1 == pairj_2) || (pairi_1 == pairj_3) ||
                                                  (pairi_2 == pairj_1) || (pairi_2 == pairj_2) || (pairi_2 == pairj_3) ||
                                                  (pairi_3 == pairj_1) || (pairi_3 == pairj_2) || (pairi_3 == pairj_3))
                        {
                                unionFind->Union(iset, jset);
                        }
                }
        }
        AllMouths = unionFind->Consolidate();
        for (i = 0; i < pockets.size(); i++)
        {
                for (j = 0; j < pockets[i].size(); j++)
                {
                        pindex = pockets[i][j];
                        delcx->DeluanayTet[pindex].PocIndex = -1;
                }
        }
        for (i = 0; i < AllMouths.size(); i++)
        {
                for (j = 0; j < AllMouths[i].size(); j++)
                {
                        pindex = AllMouths[i][j];
                        delcx->DeluanayTrigs[pindex].IsMouth = true;
                }
        }
        nmouths = AllMouths.size();
}

/*!
    \fn Mouth::MeasureMouths(DeluanayComplex *delcx, std::vector<Vertex> & vertexList)
 */
void Mouth::MeasureMouths(DeluanayComplex *delcx, std::vector<Vertex> & vertexList)
{
        int i, j, k, m;
        int imouth, idx, id;

        double ra, rb, rc;
        double surfa = 0.0, surf, s;
        double ang1 = 0.0, ang2 = 0.0, ang3 = 0.0;
        double dist1 = 0.0, dist2 = 0.0, dist3 = 0.0;

        double a[4];
        double b[4];
        double c[4];

        SurfMouth.reserve(nmouths + 1);

        for (imouth = 0; imouth < nmouths; imouth++)
        {
                SurfMouth[imouth] = 0;
                for (id = 0; id < AllMouths[imouth].size(); id++)
                {
                        idx = AllMouths[imouth][id];
                        i = delcx->DeluanayTrigs[idx].Corners[1];
                        j = delcx->DeluanayTrigs[idx].Corners[2];
                        k = delcx->DeluanayTrigs[idx].Corners[3];

                        ra = vertexList[i].Radius;
                        rb = vertexList[j].Radius;
                        rc = vertexList[k].Radius;

                        for (m = 1; m <= 3; m++)
                        {
                                a[m] = vertexList[i].Coordinates[m];
                                b[m] = vertexList[j].Coordinates[m];
                                c[m] = vertexList[k].Coordinates[m];
                        }

                        TriangleAng(a, b, c, &ang1, &ang2, &ang3, &dist1, &dist2, &dist3);

                        s = (dist1 + dist2 + dist3) / 2;
                        surf = s * (s - dist1) * (s - dist2) * (s - dist3);
                        surf = sqrt(surf);

                        SurfMouth[imouth] += surf;

                        //First check all vertices of the triangle
                        if (vertexList[i].AlphaStatus == 0)
                        {
                                surfa = PI * ra * ra * ang1;
                                SurfMouth[imouth] -= surfa;
                        }

                        if (vertexList[j].AlphaStatus == 0)
                        {
                                surfa = PI * rb * rb * ang2;
                                SurfMouth[imouth] -= surfa;
                        }

                        if (vertexList[k].AlphaStatus == 0)
                        {
                                surfa = PI * rc * rc * ang3;
                                SurfMouth[imouth] -= surfa;
                        }

                        //Now check all edges
                        if (dist1 <= (ra + rb))
                        {
                                TwoCircle(ra, rb, dist1, &surfa);
                                SurfMouth[imouth] += 0.5 * surfa;
                        }

                        if (dist2 <= (ra + rc))
                        {
                                TwoCircle(ra, rc, dist2, &surfa);
                                SurfMouth[imouth] += 0.5 * surfa;
                        }

                        if (dist3 <= (rb + rc))
                        {
                                TwoCircle(rb, rc, dist3, &surfa);
                                SurfMouth[imouth] += 0.5 * surfa;
                        }

                        if (SurfMouth[imouth] < 0) SurfMouth[imouth] = 0;
                }
        }
}

/*!
    \fn Mouth::FillTable(QTableWidget *table,QLineEdit *mouSurf)
 */
void Mouth::FillTable(QTableWidget *table,QLineEdit *mouSurf)
{
        QString insertionString;
        QTableWidgetItem *insertionItem;
        TotSurfArea = 0;

        table->setRowCount(nmouths);
        for(int i = 0;i<nmouths;i++)
        {
                QString &tempstr = insertionString.setNum(i+1);
                insertionItem = new QTableWidgetItem(tempstr);
                table->setItem (i,0,insertionItem);

                insertionItem = new QTableWidgetItem();
                insertionItem->setCheckState(Qt::Checked);
                table->setItem(i,1,insertionItem);

                tempstr = insertionString.setNum(MouthPocket[i]+1);
                insertionItem = new QTableWidgetItem(tempstr);
                table->setItem(i,2,insertionItem);

                tempstr = insertionString.setNum(SurfMouth[i],'f',4);
                insertionItem = new QTableWidgetItem(tempstr);
                table->setItem(i,3,insertionItem);

                TotSurfArea+=SurfMouth[i];
        }
        QString &ran = insertionString.setNum(TotSurfArea,'f',5);
        mouSurf->setText(ran);
}

/*!
    \fn Mouth::DrawTrig(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int i)
 */
void Mouth::DrawTrig(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int i)
{
        glBegin(GL_TRIANGLES);
        glNormal3d(delcx->DeluanayTrigs[i].Normal->X,delcx->DeluanayTrigs[i].Normal->Y,delcx->DeluanayTrigs[i].Normal->Z);
        glVertex3d(vertexList[delcx->DeluanayTrigs[i].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[i].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[i].Corners[1]].NormCoordinates[3]);
        glVertex3d(vertexList[delcx->DeluanayTrigs[i].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[i].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[i].Corners[2]].NormCoordinates[3]);
        glVertex3d(vertexList[delcx->DeluanayTrigs[i].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[i].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[i].Corners[3]].NormCoordinates[3]);
        glEnd();
}

/*!
    \fn Mouth::InitMaterial(int i)
 */
void Mouth::InitMaterial(int i)
{
        if (i == 0)
        {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[16]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[16]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[16]);
                glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[16]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);
        }
        else
        {
                //int m = i % 9+10;
                int m = i % 15+16;
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[m]);
                glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);
        }
}

/*!
    \fn Mouth::RenderMouthOfPocket(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int persistence,int pocIndex)
 */
void Mouth::RenderMouthOfPocket(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int persistence,int pocIndex)
{
        int j;
        for (int i = 0; i < nmouths; i++)
        {
                if(pocIndex - 1 == MouthPocket[i])
                {
                        InitMaterial(i);
                        for (j = 0; j < AllMouths[i].size(); j++)
                        {
                                if((delcx->DeluanayTrigs[AllMouths[i][j]].Persistence!=-1) && (delcx->DeluanayTrigs[AllMouths[i][j]].Persistence < persistence))
                                {
                                        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                                }
                                else
                                {
                                        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                                }
                                DrawTrig(delcx, vertexList, AllMouths[i][j]);
                        }
                        InitMaterial(0);
                }
        }
}

/*!
    \fn Mouth::Render(DeluanayComplex *delcx,std::vector<Vertex> & vertexList)
 */
void Mouth::Render(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int persistence)
{
        int j;
        for (int i = 0; i < nmouths; i++)
        {
                InitMaterial(i);
                for (j = 0; j < AllMouths[i].size(); j++)
                {
                        if((delcx->DeluanayTrigs[AllMouths[i][j]].Persistence!=-1) && (delcx->DeluanayTrigs[AllMouths[i][j]].Persistence < persistence))
                        {
                                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                        }
                        else
                        {
                                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                        }
                        DrawTrig(delcx, vertexList, AllMouths[i][j]);
                }
                InitMaterial(0);
        }
}
