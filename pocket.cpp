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

#include "pocket.h"

Pocket::Pocket(std::vector<Vertex> & vertexList,double center[],double scale,int trcount,int ecount)
{
        Rank = -1;
        PI = 4 * atan(1);
        spm = new SpaceFillMeasure();
        vol = new Volume(vertexList,vertexList.size (),trcount,ecount);
        unionFind = new DisJointSet();
        Radius2.reserve(vertexList.size());
        for (uint i = 1; i < vertexList.size(); i++)
        {
                Radius2[i] = pow(vertexList[i].Radius, 2.0);
        }
        for(uint i = 0;i<3;i++)
        {
                pcenter[i] = center[i];
        }
        pscale = scale;
        skin = new SkinSurface(center,scale,1);
}

Pocket::~Pocket()
{
}

/*!
    \fn Pocket::FindPockets(DeluanayComplex *delcx, std::vector <int> &sortedTet)
 */
void Pocket::FindPockets(DeluanayComplex *delcx, std::vector <int> &sortedTet)
{
        int i, j, k, itrig;
        unionFind->Clear();
        //FILE *fp = fopen("pocket union find.txt","w");
        uint TetSize = delcx->DeluanayTet.size() - delcx->redundantCount;
        pocPersistence.clear ();

        //fprintf(fp,"Tetsize = %d\n",TetSize);

        for (i = TetSize - 1; i > 0; i--)
        {
                if (delcx->DeluanayTet[sortedTet[i]].Depth == delcx->DeluanayTet.size()) continue;
          //      fprintf(fp," %d %d %d %d ",i,sortedTet[i],delcx->DeluanayTet[sortedTet[i]].AlphaStatus,delcx->DeluanayTet[sortedTet[i]].isValid);

                if (delcx->DeluanayTet[sortedTet[i]].AlphaStatus != -1)
                {
            //        fprintf(fp,"i = %d sortedTet[i] = %d cond1 fail\n",i,sortedTet[i]);
                    continue;
                }
              //  fprintf(fp,"\n");
                delcx->DeluanayTet[sortedTet[i]].ufKey = unionFind->ElementCount();
                unionFind->Add(sortedTet[i]);
        }
        for (i = TetSize - 1; i > 0; i--)
        {
                if (delcx->DeluanayTet[sortedTet[i]].ufKey == -1) continue;  //tet not in the uf structure

                for (j = 1; j <= 4; j++)
                {
                        k = delcx->DeluanayTet[sortedTet[i]].Neighbours[j];

                        if (delcx->DeluanayTet[k].ufKey == -1) continue;         //neighbour not in uf structure

                        itrig = delcx->DeluanayTet[sortedTet[i]].TetLink[j];

                        if (delcx->DeluanayTrigs[itrig].isValid == true) continue;  //useful during the modified filtration*/
                        if (delcx->DeluanayTrigs[itrig].AlphaStatus == 1)
                        {
                            //if (delcx->DeluanayTrigs[itrig].isValid == true)
                            //{
                //                fprintf(fp,"sortedTet[i] = %d k = %d itrig = %d cond2 fail\n",sortedTet[i],k,itrig);
                                continue;
                            //}
                        }
                  //      fprintf(fp,"sortedTet[i] = %d k = %d itrig = %d\n",sortedTet[i],k,itrig);
                        unionFind->Union(unionFind->FindSet(delcx->DeluanayTet[sortedTet[i]].ufKey), unionFind->FindSet(delcx->DeluanayTet[k].ufKey));
                }
        }

        AllPockets = unionFind->Consolidate();
        npockets = AllPockets.size();
        for (i = 0; i < AllPockets.size(); i++)
                {
                    //    fprintf(fp,"Pocket No %d =",i);
                        for (j = 0; j < AllPockets[i].size(); j++)
                        {
                                int pindex = AllPockets[i][j];
                                delcx->DeluanayTet[pindex].ufKey = -1;
                                delcx->DeluanayTet[pindex].PocIndex = i;
                      //          fprintf(fp," tet = %d persistence = %d alpha persistence = %lf\n",pindex,delcx->DeluanayTet[pindex].Persistence,delcx->DeluanayTet[pindex].AlphaPersistence);
                        }
                        //fprintf(fp,"\n");
                }
        //fclose(fp);

        //now calculate persistence or extended persistence of each pocket here.
        for(i = 0; i < AllPockets.size (); i++)
        {
            int persistence;
            if(AllPockets[i].size() == 1)                   //tet persistence is enough
            {
                int tindex = AllPockets[i][0];
                persistence = delcx->DeluanayTet[tindex].Persistence;
                pocPersistence.push_back (persistence);
            }
            else
            {   /*
                A tetrahedron is part of the pocket means it is not present in the alpha complex also only we need to check Rho.
                Also InComplex doesn't always work as we have not taken care of filtration modifcation in that function.
                */
                int youngestTrig = 0;
                int oldestTet = 0;
                for(j = 0; j < AllPockets[i].size(); j++)
                {
                    int trindex, tindex;
                    tindex = AllPockets[i][j];
                    if(oldestTet < delcx->DeluanayTet[tindex].Rho)
                    {
                        oldestTet = delcx->DeluanayTet[tindex].Rho;
                    }
                    for(k = 1; k < 5; k++)
                    {
                        trindex = delcx->DeluanayTet[tindex].TetLink[k];
                        if(delcx->DeluanayTrigs[trindex].AlphaStatus != 0)         //belongs to alpha complex
                        {
                            if(delcx->DeluanayTrigs[trindex].Rho)
                            {
                                if(youngestTrig < delcx->DeluanayTrigs[trindex].Rho)
                                {
                                    youngestTrig = delcx->DeluanayTrigs[trindex].Rho;
                                }
                            }
                            else
                            {
                                if(youngestTrig < delcx->DeluanayTrigs[trindex].Mu1)
                                {
                                    youngestTrig = delcx->DeluanayTrigs[trindex].Mu1;
                                }
                            }
                        }
                    }
                }
                assert(youngestTrig <= oldestTet);
                persistence = oldestTet - youngestTrig;
                pocPersistence.push_back (persistence);
            }
        }

        /*FILE *fp1 = fopen("pocpersist.txt","w");
        for(i = 0; i < pocPersistence.size (); i++)
        {
            fprintf(fp1,"pocket id = %d persistence = %d\n",i,pocPersistence[i]);
        }
        fclose(fp1);*/
}

/*!
    \fn Pocket::GetPockets()
 */
std::vector< std::vector<int> > Pocket::GetPockets()
{
    return AllPockets;
}

bool Pocket::InComplex (DeluanayComplex *delcx,std::vector<Vertex> &vertexList, unsigned int ftype, int rank, int i)
{
    switch (ftype)
        {
         case ALF_VERTEX:
          {
            if (vertexList[i].Rho)
              return (vertexList[i].Rho <= rank);
            else
              return (vertexList[i].Mu1 <= rank);
          }
         case ALF_EDGE:
          {
            if (delcx->DeluanayEdges[i].Rho)
              return (delcx->DeluanayEdges[i].Rho <= rank);
            else
              return (delcx->DeluanayEdges[i].Mu1 <= rank);
          }
         case ALF_TRIANGLE:
          {
            if (delcx->DeluanayTrigs[i].Rho)
              return (delcx->DeluanayTrigs[i].Rho <= rank);
            else
              return (delcx->DeluanayTrigs[i].Mu1 <= rank);
          }
         case ALF_TETRA:
          {
            return (delcx->DeluanayTet[i].Rho <= rank);
          }
        }
      return (false);
}

void Pocket::PocketProperties (DeluanayComplex *delcx, std::vector<Vertex> &vertexList,int rank)
{
    int i,j,k,l;
    int ipocket,id,idx,pair1, pair2, pair3, pair4, pair5, pair6;
    int trig1, trig2, trig3, trig4;
    double tetravolume,surf,volume;

    SurfPocket.reserve(npockets+1);
    VolPocket.reserve(npockets+1);

    for (i = 1; i < vertexList.size(); i++)
    {
            Radius2[i] = pow(vertexList[i].Radius, 2.0);
    }

    for (ipocket = 0; ipocket < npockets; ipocket++)
    {
            for (id = 0; id < AllPockets[ipocket].size(); id++)
            {
                idx = AllPockets[ipocket][id];

                i = delcx->DeluanayTet[idx].Corners[1];
                j = delcx->DeluanayTet[idx].Corners[2];
                k = delcx->DeluanayTet[idx].Corners[3];
                l = delcx->DeluanayTet[idx].Corners[4];

                trig1 = delcx->DeluanayTet[idx].TetLink[4];
                trig2 = delcx->DeluanayTet[idx].TetLink[3];
                trig3 = delcx->DeluanayTet[idx].TetLink[2];
                trig4 = delcx->DeluanayTet[idx].TetLink[1];

                pair1 = delcx->DeluanayTrigs[trig1].TrigLink[3];
                pair2 = delcx->DeluanayTrigs[trig1].TrigLink[2];
                pair4 = delcx->DeluanayTrigs[trig1].TrigLink[1];
                pair3 = delcx->DeluanayTrigs[trig2].TrigLink[2];
                pair5 = delcx->DeluanayTrigs[trig2].TrigLink[1];
                pair6 = delcx->DeluanayTrigs[trig3].TrigLink[1];

                 //find initial volume

                tetravolume = vol->DoTetraVolume (vertexList,i,j,k,l);
                VolPocket[ipocket] += tetravolume;

                if(InComplex (delcx,vertexList,ALF_VERTEX,rank,i))
                {
                    vol->DoTetraVertex (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] -= volume;
                    SurfPocket[ipocket] += surf;
                }
                if(InComplex (delcx,vertexList,ALF_VERTEX,rank,j))
                {
                    vol->DoTetraVertex (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] -= volume;
                    SurfPocket[ipocket] += surf;
                }
                if(InComplex (delcx,vertexList,ALF_VERTEX,rank,k))
                {
                    vol->DoTetraVertex (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] -= volume;
                    SurfPocket[ipocket] += surf;
                }
                if(InComplex (delcx,vertexList,ALF_VERTEX,rank,l))
                {
                    vol->DoTetraVertex (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] -= volume;
                    SurfPocket[ipocket] += surf;
                }

                if(InComplex (delcx,vertexList,ALF_EDGE,rank,pair1))
                {
                    vol->DoTetraEdge (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] += volume;
                    SurfPocket[ipocket] -= surf;
                }
                if(InComplex (delcx,vertexList,ALF_EDGE,rank,pair2))
                {
                    vol->DoTetraEdge (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] += volume;
                    SurfPocket[ipocket] -= surf;
                }
                if(InComplex (delcx,vertexList,ALF_EDGE,rank,pair3))
                {
                    vol->DoTetraEdge (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] += volume;
                    SurfPocket[ipocket] -= surf;
                }
                if(InComplex (delcx,vertexList,ALF_EDGE,rank,pair4))
                {
                    vol->DoTetraEdge (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] += volume;
                    SurfPocket[ipocket] -= surf;
                }
                if(InComplex (delcx,vertexList,ALF_EDGE,rank,pair5))
                {
                    vol->DoTetraEdge (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] += volume;
                    SurfPocket[ipocket] -= surf;
                }
                if(InComplex (delcx,vertexList,ALF_EDGE,rank,pair6))
                {
                    vol->DoTetraEdge (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] += volume;
                    SurfPocket[ipocket] -= surf;
                }

                if(InComplex (delcx,vertexList,ALF_TRIANGLE,rank,trig1))
                {
                    vol->DoTetraTriangle (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] -= volume;
                    SurfPocket[ipocket] += surf;
                }
                if(InComplex (delcx,vertexList,ALF_TRIANGLE,rank,trig2))
                {
                    vol->DoTetraTriangle (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] -= volume;
                    SurfPocket[ipocket] += surf;
                }
                if(InComplex (delcx,vertexList,ALF_TRIANGLE,rank,trig3))
                {
                    vol->DoTetraTriangle (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] -= volume;
                    SurfPocket[ipocket] += surf;
                }
                if(InComplex (delcx,vertexList,ALF_TRIANGLE,rank,trig4))
                {
                    vol->DoTetraTriangle (vertexList,i,j,k,l,1,&surf,&volume);
                    VolPocket[ipocket] -= volume;
                    SurfPocket[ipocket] += surf;
                }
            }
    }
    for (i = 0; i < npockets; i++)
    {
            if (SurfPocket[i] < 0) SurfPocket[i] = 0;
            if (VolPocket[i] < 0) VolPocket[i] = 0;
    }
}

/*!
    \fn Pocket::MeasurePockets(DeluanayComplex *delcx, std::vector<Vertex> & vertexList)
 */
void Pocket::MeasurePockets(DeluanayComplex *delcx, std::vector<Vertex> & vertexList)
{
        int i, j, k, l, m, option;
        int ipocket, idx, id, pair1, pair2, pair3, pair4, pair5, pair6;
        int trig1, trig2, trig3, trig4;

        double tetra_volume;
        double ra, rb, rc, rd, ra2, rb2, rc2, rd2;
        double wa, wb, wc, wd;
        double surfa = 0.0, surfb = 0.0, surfc = 0.0;
        double vola = 0.0, volb = 0.0, volc = 0.0, coef;
        double dist1 = 0.0, dist2 = 0.0, dist3 = 0.0, dist4 = 0.0, dist5 = 0.0, dist6 = 0.0;
        double d2_1 = 0.0, d2_2 = 0.0, d2_3 = 0.0, d2_4 = 0.0, d2_5 = 0.0, d2_6 = 0.0;
        double ang1 = 0.0, ang2 = 0.0, ang3 = 0.0, ang4 = 0.0, ang5 = 0.0, ang6 = 0.0;

        double a[4];
        double b[4];
        double c[4];
        double d[4];

        double dsurfa[4][3];
        double dsurfb[4][3];
        double dvola[4][3];
        double dvolb[4][3];

        std::vector<int> index;
        std::vector<int> rank;
        std::vector<double> temp;

        index.reserve(npockets+1);
        rank.reserve(npockets+1);
        temp.reserve(npockets+1);
        SurfPocket.reserve(npockets+1);
        VolPocket.reserve(npockets+1);

        option = 0;

        for (i = 1; i < vertexList.size(); i++)
        {
                Radius2[i] = pow(vertexList[i].Radius, 2.0);
        }

        for (i = 0;i<npockets+1;i++)
        {
            SurfPocket[i] = 0.0;
            VolPocket[i] = 0.0;
        }

       // FILE *fp = fopen("pocvol.txt","w");

        for (ipocket = 0; ipocket < npockets; ipocket++)
        {
                for (id = 0; id < AllPockets[ipocket].size(); id++)
                {
                        idx = AllPockets[ipocket][id];
                        //fprintf(fp,"tetindex = %d\n",idx);
                        i = delcx->DeluanayTet[idx].Corners[1];
                        j = delcx->DeluanayTet[idx].Corners[2];
                        k = delcx->DeluanayTet[idx].Corners[3];
                        l = delcx->DeluanayTet[idx].Corners[4];

                        ra = vertexList[i].Radius;
                        rb = vertexList[j].Radius;
                        rc = vertexList[k].Radius;
                        rd = vertexList[l].Radius;

                        ra2 = Radius2[i];
                        rb2 = Radius2[j];
                        rc2 = Radius2[k];
                        rd2 = Radius2[l];

                        wa = (1.0/2.0) * vertexList[i].Weight;
                        wb = (1.0/2.0) * vertexList[j].Weight;
                        wc = (1.0/2.0) * vertexList[k].Weight;
                        wd = (1.0/2.0) * vertexList[l].Weight;

                        trig1 = delcx->DeluanayTet[idx].TetLink[4];
                        trig2 = delcx->DeluanayTet[idx].TetLink[3];
                        trig3 = delcx->DeluanayTet[idx].TetLink[2];
                        trig4 = delcx->DeluanayTet[idx].TetLink[1];

                        pair1 = delcx->DeluanayTrigs[trig1].TrigLink[3];
                        pair2 = delcx->DeluanayTrigs[trig1].TrigLink[2];
                        pair4 = delcx->DeluanayTrigs[trig1].TrigLink[1];
                        pair3 = delcx->DeluanayTrigs[trig2].TrigLink[2];
                        pair5 = delcx->DeluanayTrigs[trig2].TrigLink[1];
                        pair6 = delcx->DeluanayTrigs[trig3].TrigLink[1];

                        for (m = 1; m <= 3; m++)
                        {
                                a[m] = vertexList[i].Coordinates[m];
                                b[m] = vertexList[j].Coordinates[m];
                                c[m] = vertexList[k].Coordinates[m];
                                d[m] = vertexList[l].Coordinates[m];
                        }

                        spm->Tetra6Dihedral(a, b, c, d, &ang1, &ang2, &ang4, &ang3, &ang5, &ang6);

                        tetra_volume = spm->TetraVolume(a, b, c, d);

                        //fprintf(fp,"tetravolume = %f\n",tetra_volume);
                        VolPocket[ipocket] += tetra_volume;
                        //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                        //fprintf(fp,"\n");
                        //check all vertices of the tetrahedron
                        if (vertexList[i].AlphaStatus == 0)
                        {
                                coef = 0.5 * (ang1 + ang2 + ang3 - 0.5);
                                surfa = 4 * PI * Radius2[i] * coef;
                                vola = surfa * vertexList[i].Radius / 3.0;

                                SurfPocket[ipocket] += surfa;
                                VolPocket[ipocket] -= vola;

                                //fprintf(fp,"corner1 volume = %lf\n",vola);
                                //fprintf(fp,"corner1 surface =%lf\n",surfa);
                                //fprintf(fp,"pocketvolume = %lf\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %lf\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (vertexList[j].AlphaStatus == 0)
                        {
                                coef = 0.5 * (ang1 + ang4 + ang5 - 0.5);
                                surfa = 4 * PI * Radius2[j] * coef;
                                vola = surfa * vertexList[j].Radius / 3.0;

                                SurfPocket[ipocket] += surfa;
                                VolPocket[ipocket] -= vola;

                                //fprintf(fp,"corner2 volume = %lf\n",vola);
                                //fprintf(fp,"corner2 surface =%lf\n",surfa);
                                //fprintf(fp,"pocketvolume = %lf\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %lf\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (vertexList[k].AlphaStatus == 0)
                        {
                                coef = 0.5 * (ang2 + ang4 + ang6 - 0.5);
                                surfa = 4 * PI * Radius2[k] * coef;
                                vola = surfa * vertexList[k].Radius / 3.0;

                                SurfPocket[ipocket] += surfa;
                                VolPocket[ipocket] -= vola;

                                //fprintf(fp,"corner3 volume = %lf\n",vola);
                                //fprintf(fp,"corner3 surface =%lf\n",surfa);
                                //fprintf(fp,"pocketvolume = %lf\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %lf\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (vertexList[l].AlphaStatus == 0)
                        {
                                coef = 0.5 * (ang3 + ang5 + ang6 - 0.5);
                                surfa = 4 * PI * Radius2[l] * coef;
                                vola = surfa * vertexList[l].Radius / 3.0;

                                SurfPocket[ipocket] += surfa;
                                VolPocket[ipocket] -= vola;

                                //fprintf(fp,"corner4 volume = %lf\n",vola);
                                //fprintf(fp,"corner4 surface =%lf\n",surfa);
                                //fprintf(fp,"pocketvolume = %lf\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %lf\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        //check all edges
                        if (delcx->DeluanayEdges[pair1].AlphaStatus == 1)
                        {
                                spm->Distance2(vertexList, i, j, &d2_1);
                                dist1 = sqrt(d2_1);

                                spm->TwoSphereVol(a, b, ra, ra2, rb, rb2, dist1, d2_1, &surfa, &surfb, &vola, &volb, dsurfa, dsurfb, dvola, dvolb, option);

                                SurfPocket[ipocket] -= ang1 * (surfa + surfb);
                                VolPocket[ipocket] += ang1 * (vola + volb);

                                //fprintf(fp,"angle1 = %f\n",ang1);
                                //fprintf(fp,"edge1 volume = %f\n",vola + volb);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (delcx->DeluanayEdges[pair2].AlphaStatus == 1)
                        {
                                spm->Distance2(vertexList, i, k, &d2_2);
                                dist2 = sqrt(d2_2);

                                spm->TwoSphereVol(a, c, ra, ra2, rc, rc2, dist2, d2_2, &surfa, &surfb, &vola, &volb, dsurfa, dsurfb, dvola, dvolb, option);

                                SurfPocket[ipocket] -= ang2 * (surfa + surfb);
                                VolPocket[ipocket] += ang2 * (vola + volb);

                                //fprintf(fp,"angle2 = %f\n",ang2);
                                //fprintf(fp,"edge2 volume = %f\n",vola + volb);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (delcx->DeluanayEdges[pair3].AlphaStatus == 1)
                        {
                                spm->Distance2(vertexList, i, l, &d2_3);
                                dist3 = sqrt(d2_3);

                                spm->TwoSphereVol(a, d, ra, ra2, rd, rd2, dist3, d2_3, &surfa, &surfb, &vola, &volb, dsurfa, dsurfb, dvola, dvolb, option);

                                SurfPocket[ipocket] -= ang3 * (surfa + surfb);
                                VolPocket[ipocket] += ang3 * (vola + volb);

                                //fprintf(fp,"angle3 = %f\n",ang3);
                                //fprintf(fp,"edge3 volume = %f\n",vola + volb);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (delcx->DeluanayEdges[pair4].AlphaStatus == 1)
                        {
                                spm->Distance2(vertexList, j, k, &d2_4);
                                dist4 = sqrt(d2_4);

                                spm->TwoSphereVol(b, c, rb, rb2, rc, rc2, dist4, d2_4, &surfa, &surfb, &vola, &volb, dsurfa, dsurfb, dvola, dvolb, option);

                                SurfPocket[ipocket] -= ang4 * (surfa + surfb);
                                VolPocket[ipocket] += ang4 * (vola + volb);

                                //fprintf(fp,"angle4 = %f\n",ang4);
                                //fprintf(fp,"edge4 volume = %f\n",vola + volb);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (delcx->DeluanayEdges[pair5].AlphaStatus == 1)
                        {
                                spm->Distance2(vertexList, j, l, &d2_5);
                                dist5 = sqrt(d2_5);

                                spm->TwoSphereVol(b, d, rb, rb2, rd, rd2, dist5, d2_5, &surfa, &surfb, &vola, &volb, dsurfa, dsurfb, dvola, dvolb, option);

                                SurfPocket[ipocket] -= ang5 * (surfa + surfb);
                                VolPocket[ipocket] += ang5 * (vola + volb);

                                //fprintf(fp,"angle5 = %f\n",ang5);
                                //fprintf(fp,"edge5 volume = %f\n",vola + volb);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (delcx->DeluanayEdges[pair6].AlphaStatus == 1)
                        {
                                spm->Distance2(vertexList, k, l, &d2_6);
                                dist6 = sqrt(d2_6);

                                spm->TwoSphereVol(c, d, rc, rc2, rd, rd2, dist6, d2_6, &surfa, &surfb, &vola, &volb, dsurfa, dsurfb, dvola, dvolb, option);

                                SurfPocket[ipocket] -= ang6 * (surfa + surfb);
                                VolPocket[ipocket] += ang6 * (vola + volb);

                                //fprintf(fp,"angle6 = %f\n",ang6);
                                //fprintf(fp,"edge6 volume = %f\n",vola + volb);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        //Finally check faces
                        if (delcx->DeluanayTrigs[trig1].AlphaStatus == 1)
                        {
                                spm->ThreeSphereVol(a, b, c, ra, rb, rc, ra2, rb2, rc2, wa, wb, wc, dist1, dist2, dist4, d2_1, d2_2, d2_4,
                                                &surfa, &surfb, &surfc, &vola, &volb, &volc);

                                SurfPocket[ipocket] += 0.5 * (surfa + surfb + surfc);
                                VolPocket[ipocket] -= 0.5 * (vola + volb + volc);


                                //fprintf(fp,"trig1 volume = %f\n",vola + volb + volc);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (delcx->DeluanayTrigs[trig2].AlphaStatus == 1)
                        {
                                spm->ThreeSphereVol(a, b, d, ra, rb, rd, ra2, rb2, rd2, wa, wb, wd, dist1, dist3, dist5, d2_1, d2_3, d2_5,
                                                &surfa, &surfb, &surfc, &vola, &volb, &volc);

                                SurfPocket[ipocket] += 0.5 * (surfa + surfb + surfc);
                                VolPocket[ipocket] -= 0.5 * (vola + volb + volc);

                                //fprintf(fp,"trig2 volume = %f\n",vola + volb + volc);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (delcx->DeluanayTrigs[trig3].AlphaStatus == 1)
                        {
                                spm->ThreeSphereVol(a, c, d, ra, rc, rd, ra2, rc2, rd2, wa, wc, wd, dist2, dist3, dist6, d2_2, d2_3, d2_6,
                                                &surfa, &surfb, &surfc, &vola, &volb, &volc);

                                SurfPocket[ipocket] += 0.5 * (surfa + surfb + surfc);
                                VolPocket[ipocket] -= 0.5 * (vola + volb + volc);

                                //fprintf(fp,"trig3 volume = %f\n",vola + volb + volc);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }

                        if (delcx->DeluanayTrigs[trig4].AlphaStatus == 1)
                        {
                                spm->ThreeSphereVol(b, c, d, rb, rc, rd, rb2, rc2, rd2, wb, wc, wd, dist4, dist5, dist6, d2_4, d2_5, d2_6,
                                                &surfa, &surfb, &surfc, &vola, &volb, &volc);

                                SurfPocket[ipocket] += 0.5 * (surfa + surfb + surfc);
                                VolPocket[ipocket] -= 0.5 * (vola + volb + volc);

                                //fprintf(fp,"trig4 volume = %f\n",vola + volb + volc);
                                //fprintf(fp,"pocketvolume = %f\n",VolPocket[ipocket]);
                                //fprintf(fp,"pocarea = %f\n",SurfPocket[ipocket]);
                                //fprintf(fp,"\n");
                        }
                }
                //fprintf(fp,"----------------------------------\n----------------------------------\n");
        }
        //fclose(fp);
        for (i = 0; i < npockets; i++)
        {
                if (SurfPocket[i] < 0) SurfPocket[i] = 0;
                if (VolPocket[i] < 0) VolPocket[i] = 0;
        }
}

/*!
    \fn Pocket::FillTable(QTableWidget *table,std::vector<int> &PocketNMouths)
 */
void Pocket::FillTable(QTableWidget *table,std::vector<int> & PocketNMouths,QLineEdit *pocVol,QLineEdit *pocSurf)
{
        QString insertionString;
        QTableWidgetItem *insertionItem;
        TotVolume = 0;
        double TotPocVolume = 0.0;
        double TotVoidVolume = 0.0;
        TotSurfArea = 0;
        table->setRowCount(npockets);
        FILE *fp = fopen("pocvolgraph.txt","w");
        for(int i = 0;i<npockets;i++)
        {
                QString &tempstr = insertionString.setNum(i+1);
                insertionItem = new QTableWidgetItem(tempstr);
                table->setItem (i,0,insertionItem);

                insertionItem = new QTableWidgetItem();
                insertionItem->setCheckState(Qt::Checked);
                table->setItem(i,1,insertionItem);

                tempstr = insertionString.setNum(PocketNMouths[i]);
                insertionItem = new QTableWidgetItem(tempstr);
                table->setItem(i,2,insertionItem);

                tempstr = insertionString.setNum(VolPocket[i],'f',4);
                insertionItem = new QTableWidgetItem(tempstr);
                table->setItem(i,3,insertionItem);

                tempstr = insertionString.setNum(SurfPocket[i],'f',4);
                insertionItem = new QTableWidgetItem(tempstr);
                table->setItem(i,4,insertionItem);

                if(PocketNMouths[i] == 0)               //voids
                {
                    TotVoidVolume += VolPocket[i];
                }
                else                                    //pockets
                {
                    TotPocVolume += VolPocket[i];
                }

                TotVolume += VolPocket[i];
                TotSurfArea += SurfPocket[i];
        }

        QString &ran = insertionString.setNum(TotVolume,'f',5);
        pocVol->setText(ran);

        ran = insertionString.setNum(TotSurfArea,'f',5);
        pocSurf->setText(ran);
        //fprintf(fp,"%lf %lf %lf",TotVolume,TotPocVolume,TotVoidVolume);
        //fclose(fp);
}

/*!
    \fn Pocket::DrawTet(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int l)
 */
void Pocket::DrawTet(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,int l,bool mouth)
{
        int a,b,c,d;

        a = delcx->DeluanayTet[l].TetLink[1];
        b = delcx->DeluanayTet[l].TetLink[2];
        c = delcx->DeluanayTet[l].TetLink[3];
        d = delcx->DeluanayTet[l].TetLink[4];

        if(!mouth)
        {
            glBegin(GL_TRIANGLES);
            glNormal3d(delcx->DeluanayTrigs[a].Normal->X,delcx->DeluanayTrigs[a].Normal->Y,delcx->DeluanayTrigs[a].Normal->Z);
            glVertex3d(vertexList[delcx->DeluanayTrigs[a].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[a].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[a].Corners[1]].NormCoordinates[3]);
            glVertex3d(vertexList[delcx->DeluanayTrigs[a].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[a].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[a].Corners[2]].NormCoordinates[3]);
            glVertex3d(vertexList[delcx->DeluanayTrigs[a].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[a].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[a].Corners[3]].NormCoordinates[3]);
            glEnd();

            glBegin(GL_TRIANGLES);
            glNormal3d(delcx->DeluanayTrigs[b].Normal->X,delcx->DeluanayTrigs[b].Normal->Y,delcx->DeluanayTrigs[b].Normal->Z);
            glVertex3d(vertexList[delcx->DeluanayTrigs[b].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[b].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[b].Corners[1]].NormCoordinates[3]);
            glVertex3d(vertexList[delcx->DeluanayTrigs[b].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[b].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[b].Corners[2]].NormCoordinates[3]);
            glVertex3d(vertexList[delcx->DeluanayTrigs[b].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[b].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[b].Corners[3]].NormCoordinates[3]);
            glEnd();

            glBegin(GL_TRIANGLES);
            glNormal3d(delcx->DeluanayTrigs[c].Normal->X,delcx->DeluanayTrigs[c].Normal->Y,delcx->DeluanayTrigs[c].Normal->Z);
            glVertex3d(vertexList[delcx->DeluanayTrigs[c].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[c].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[c].Corners[1]].NormCoordinates[3]);
            glVertex3d(vertexList[delcx->DeluanayTrigs[c].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[c].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[c].Corners[2]].NormCoordinates[3]);
            glVertex3d(vertexList[delcx->DeluanayTrigs[c].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[c].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[c].Corners[3]].NormCoordinates[3]);
            glEnd();

            glBegin(GL_TRIANGLES);
            glNormal3d(delcx->DeluanayTrigs[d].Normal->X,delcx->DeluanayTrigs[d].Normal->Y,delcx->DeluanayTrigs[d].Normal->Z);
            glVertex3d(vertexList[delcx->DeluanayTrigs[d].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[d].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[d].Corners[1]].NormCoordinates[3]);
            glVertex3d(vertexList[delcx->DeluanayTrigs[d].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[d].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[d].Corners[2]].NormCoordinates[3]);
            glVertex3d(vertexList[delcx->DeluanayTrigs[d].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[d].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[d].Corners[3]].NormCoordinates[3]);
            glEnd();
        }
        else
        {
            if(!delcx->DeluanayTrigs[a].IsMouth)
            {
                    glBegin(GL_TRIANGLES);
                    glNormal3d(delcx->DeluanayTrigs[a].Normal->X,delcx->DeluanayTrigs[a].Normal->Y,delcx->DeluanayTrigs[a].Normal->Z);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[a].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[a].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[a].Corners[1]].NormCoordinates[3]);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[a].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[a].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[a].Corners[2]].NormCoordinates[3]);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[a].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[a].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[a].Corners[3]].NormCoordinates[3]);
                    glEnd();
            }

            if(!delcx->DeluanayTrigs[b].IsMouth)
            {
                    glBegin(GL_TRIANGLES);
                    glNormal3d(delcx->DeluanayTrigs[b].Normal->X,delcx->DeluanayTrigs[b].Normal->Y,delcx->DeluanayTrigs[b].Normal->Z);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[b].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[b].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[b].Corners[1]].NormCoordinates[3]);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[b].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[b].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[b].Corners[2]].NormCoordinates[3]);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[b].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[b].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[b].Corners[3]].NormCoordinates[3]);
                    glEnd();
            }

            if(!delcx->DeluanayTrigs[c].IsMouth)
            {
                    glBegin(GL_TRIANGLES);
                    glNormal3d(delcx->DeluanayTrigs[c].Normal->X,delcx->DeluanayTrigs[c].Normal->Y,delcx->DeluanayTrigs[c].Normal->Z);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[c].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[c].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[c].Corners[1]].NormCoordinates[3]);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[c].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[c].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[c].Corners[2]].NormCoordinates[3]);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[c].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[c].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[c].Corners[3]].NormCoordinates[3]);
                    glEnd();
            }

            if(!delcx->DeluanayTrigs[d].IsMouth)
            {
                    glBegin(GL_TRIANGLES);
                    glNormal3d(delcx->DeluanayTrigs[d].Normal->X,delcx->DeluanayTrigs[d].Normal->Y,delcx->DeluanayTrigs[d].Normal->Z);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[d].Corners[1]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[d].Corners[1]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[d].Corners[1]].NormCoordinates[3]);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[d].Corners[2]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[d].Corners[2]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[d].Corners[2]].NormCoordinates[3]);
                    glVertex3d(vertexList[delcx->DeluanayTrigs[d].Corners[3]].NormCoordinates[1],vertexList[delcx->DeluanayTrigs[d].Corners[3]].NormCoordinates[2],vertexList[delcx->DeluanayTrigs[d].Corners[3]].NormCoordinates[3]);
                    glEnd();
            }
        }
}

/*!
    \fn Pocket::InitMaterial(int i)
 */
void Pocket::InitMaterial(int i)
{
        if (i == 0)
        {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[i]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[i]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[i]);
                glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[i]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);
        }
        else
        {
                //int m = i % 9;
                int m = i % 15;
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[m]);
                glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);
        }
}

/*!
    \fn Pocket::RenderSingle(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,std::vector<int> & PocketNMouths,bool pc,bool vo,bool sk,int persistence,int rank,,int pocIndex)
 */
void Pocket::RenderSingle(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,std::vector<int> & PocketNMouths,bool allPockets,bool onlyPockets,bool onlyVoids,bool skinSurface,int persistence,int rank,int pocIndex,bool smoothShading,bool skinWireFrame,bool mouth,bool pocketWireFrame)
{
        int i = pocIndex-1;

        if(!skinSurface)
        {
            if(allPockets)
            {
                InitMaterial(i);
                for (int j = 0; j < AllPockets[i].size(); j++)
                {
                    if((pocketWireFrame) || ((delcx->DeluanayTet[AllPockets[i][j]].Persistence!=-1) && (delcx->DeluanayTet[AllPockets[i][j]].Persistence < persistence)))
                    {
                            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    }
                    else
                    {
                            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    }
                    DrawTet(delcx, vertexList, AllPockets[i][j],mouth);
                }
            }
            else
            {
                if (onlyVoids && !onlyPockets)
                {
                        if (PocketNMouths[i] != 0) return;
                }
                else if (onlyPockets && !onlyVoids)
                {
                        if (PocketNMouths[i] == 0) return;
                }
                InitMaterial(i);
                for (int j = 0; j < AllPockets[i].size(); j++)
                {
                    if((pocketWireFrame) || ((delcx->DeluanayTet[AllPockets[i][j]].Persistence!=-1) && (delcx->DeluanayTet[AllPockets[i][j]].Persistence < persistence)))
                    {
                            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    }
                    else
                    {
                            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    }
                    DrawTet(delcx, vertexList, AllPockets[i][j],mouth);
                }
            }
        }
        else
        {
            skin->Draw(smoothShading,skinWireFrame,i);
        }
}

/*!
    \fn Pocket::Render(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,std::vector<int> & PocketNMouths,bool pc,bool vo,bool sk,int persistence,int rank)
 */
void Pocket::Render(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,std::vector<int> & PocketNMouths,bool allPockets,bool onlyPockets,bool onlyVoids,bool skinSurface,int persistence,int rank,bool smoothShading,bool skinWireFrame,bool mouth,bool pocketWireFrame)
{
        std::vector<Vertex> skinList;
        FILE *fp;
        int  j;


        if(!skinSurface)
        {
                if(allPockets)
                {
                    for (int i = 0; i < npockets; i++)
                    {
                        InitMaterial(i);
                        for (j = 0; j < AllPockets[i].size(); j++)
                        {
                            if((pocketWireFrame) || ((delcx->DeluanayTet[AllPockets[i][j]].Persistence!=-1) && (delcx->DeluanayTet[AllPockets[i][j]].Persistence < persistence)))
                            {
                                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                            }
                            else
                            {
                                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                            }
                            DrawTet(delcx, vertexList, AllPockets[i][j],mouth);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < npockets; i++)
                    {
                        if (onlyVoids && !onlyPockets)
                        {
                                if (PocketNMouths[i] != 0) continue;
                        }
                        else if (onlyPockets && !onlyVoids)
                        {
                                if (PocketNMouths[i] == 0) continue;
                        }
                        InitMaterial(i);
                        for (j = 0; j < AllPockets[i].size(); j++)
                        {
                                if((pocketWireFrame) || ((delcx->DeluanayTet[AllPockets[i][j]].Persistence!=-1) && (delcx->DeluanayTet[AllPockets[i][j]].Persistence < persistence)))
                                {
                                        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                                }
                                else
                                {
                                        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                                }
                                DrawTet(delcx, vertexList, AllPockets[i][j],mouth);
                        }
                    }
                }
        }
        else
        {
                if(rank != Rank)
                {
                        system("rm -f skin_lev0.off");
                        system("rm -f skin");
                        skinList.clear();

                        for (int i = 0; i < npockets; i++)
                        {
                                for (j = 0; j < AllPockets[i].size(); j++)
                                {
                                        if((delcx->DeluanayTet[AllPockets[i][j]].Persistence!=-1) && (delcx->DeluanayTet[AllPockets[i][j]].Persistence >= persistence))
                                        {
                                                delcx->DeluanayTet[AllPockets[i][j]].center4(vertexList,skinList);
                                        }
                                }
                                int temp = skinList.size();
                                if(temp < 4)
                                {
                                        switch(temp)
                                        {
                                                case 3:
                                                {
                                                        Vertex vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]-0.723;
                                                        vert.Coordinates[2] = vert.Coordinates[2]-1.474;
                                                        vert.Coordinates[3] = vert.Coordinates[3]+0.148;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);
                                                }
                                                break;
                                                case 2:
                                                {
                                                        Vertex vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]-0.723;
                                                        vert.Coordinates[2] = vert.Coordinates[2]-1.474;
                                                        vert.Coordinates[3] = vert.Coordinates[3]+0.148;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);

                                                        vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]+0.834;
                                                        vert.Coordinates[2] = vert.Coordinates[2]+0.487;
                                                        vert.Coordinates[3] = vert.Coordinates[3]-0.261;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);
                                                }
                                                break;
                                                case 1:
                                                {
                                                        Vertex vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]-0.723;
                                                        vert.Coordinates[2] = vert.Coordinates[2]-1.474;
                                                        vert.Coordinates[3] = vert.Coordinates[3]+0.148;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);

                                                        vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]+0.834;
                                                        vert.Coordinates[2] = vert.Coordinates[2]+0.487;
                                                        vert.Coordinates[3] = vert.Coordinates[3]-0.261;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);

                                                        vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]+0.078;
                                                        vert.Coordinates[2] = vert.Coordinates[2]-0.243;
                                                        vert.Coordinates[3] = vert.Coordinates[3]-0.967;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);
                                                }
                                                break;
                                                case 0:
                                                {
                                                        Vertex vert(0.0,0.0,0.0,-1,-1,0.0);
                                                        skinList.push_back(vert);

                                                        vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]-0.723;
                                                        vert.Coordinates[2] = vert.Coordinates[2]-1.474;
                                                        vert.Coordinates[3] = vert.Coordinates[3]+0.148;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);

                                                        vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]+0.834;
                                                        vert.Coordinates[2] = vert.Coordinates[2]+0.487;
                                                        vert.Coordinates[3] = vert.Coordinates[3]-0.261;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);

                                                        vert = skinList[skinList.size()-1];
                                                        vert.Coordinates[1] = vert.Coordinates[1]+0.078;
                                                        vert.Coordinates[2] = vert.Coordinates[2]-0.243;
                                                        vert.Coordinates[3] = vert.Coordinates[3]-0.967;
                                                        vert.Radius = -1.0;
                                                        skinList.push_back(vert);
                                                }
                                                break;
                                        }
                                }
                                assert(skinList.size()>=4);
                                fp = fopen("skin","w");
                                fprintf(fp,"%d\n",skinList.size());
                                fprintf(fp,"#junk\n");
                                for(uint ll = 0 ; ll < skinList.size() ; ll++)
                                {
                                        fprintf(fp,"%d %f %f %f %f\n",ll+1,skinList[ll].Coordinates[1],skinList[ll].Coordinates[2],skinList[ll].Coordinates[3],skinList[ll].Radius);
                                }
                                fclose(fp);
                                skinList.clear();
                                system("./smesh skin -s skin.off -t skin.tet");
                                skin->Read("skin_lev0.off",i);
                        }
                        skin->Process();
                }
                skin->Draw(smoothShading,skinWireFrame);
                Rank = rank;
        }
}
