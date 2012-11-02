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

#include "volume.h"

Volume::Volume()
{
        spm = new SpaceFillMeasure();
        test = 0;
        PI = 4 * atan(1);
}

/*!
    \fn Volume::Volume(int vcount,int trcount,int ecount)
 */
Volume::Volume(std::vector<Vertex> &vertexList,int vcount,int trcount,int ecount)
{
        spm = new SpaceFillMeasure();
        met = new Metric(vertexList,vcount);
        test = 0;
        PI = 4 * atan(1);
        IndVolume.reserve(vcount);
        IndSurfArea.reserve(vcount);

        for(int i = 0;i< vcount;i++)
        {
            IndVolume[i] = 0.0;
            IndSurfArea[i] = 0.0;
        }

        DerivSurf.reserve(vcount);
        DerivVol.reserve(vcount);

        trig_eps.reserve(trcount);
        trig_sh.reserve(trcount);
        trig_dual1.reserve(trcount);
        trig_dual2.reserve(trcount);
        trig_ang.reserve(trcount);

        nlink_trig.reserve(trcount);
        link_trig.reserve(trcount);

        distpair.reserve(ecount);
        distpair2.reserve(ecount);

        edgeCoef.reserve(ecount);

        Radius2.reserve(vcount);

        CoefAsp.reserve(vcount);
}

Volume::~Volume()
{
        IndVolume.clear();
        IndSurfArea.clear();
        DerivSurf.clear();
        DerivVol.clear();

        trig_eps.clear();
        trig_sh.clear();
        trig_dual1.clear();
        trig_dual2.clear();
        trig_ang.clear();

        nlink_trig.clear();
        link_trig.clear();

        distpair.clear();
        distpair2.clear();

        edgeCoef.clear();

        Radius2.clear();

        CoefAsp.clear();
}

void Volume::VertexProperties(std::vector<Vertex> &vertexList,int idx,FILE *fp1)
{
    double ball_volume, ball_area;

    ball_volume = met->BallVolume (vertexList,idx);
    ball_area = met->BallArea (idx);

    IndSurfArea[idx] = ball_area;
    IndVolume[idx] = ball_volume;

    fprintf(fp1,"point = %d\n",idx);
    fprintf(fp1,"radius2 = %lf\n",Radius2[idx]);
    fprintf(fp1,"surf =  %lf\n",IndSurfArea[idx]);
    fprintf(fp1,"vol =  %lf\n",IndVolume[idx]);
    fprintf(fp1,"\n");
}

void Volume::EdgeProperties(std::vector<Vertex> &vertexList,int i,int j,FILE *fp1)
{
    double aux_i=0.0,aux_j=0.0,vol_i=0.0,vol_j=0.0;

    aux_i = met->CapArea (vertexList,i,j);
    aux_j = met->CapArea (vertexList,j,i);

    vol_i = met->CapVolume (vertexList,i,j);
    vol_j = met->CapVolume (vertexList,j,i);

    fprintf(fp1,"surfi before = %lf\n",IndSurfArea[i]);
    fprintf(fp1,"surfj before = %lf\n",IndSurfArea[j]);
    fprintf(fp1,"voli before = %lf\n",IndVolume[i]);
    fprintf(fp1,"volj before = %lf\n",IndVolume[j]);

    IndSurfArea[i] -= aux_i;
    IndSurfArea[j] -= aux_j;

    IndVolume[i] -= vol_i;
    IndVolume[j] -= vol_j;

    fprintf(fp1,"------------\n");

    fprintf(fp1,"surfi after = %lf\n",IndSurfArea[i]);
    fprintf(fp1,"surfj after = %lf\n",IndSurfArea[j]);
    fprintf(fp1,"voli after = %lf\n",IndVolume[i]);
    fprintf(fp1,"volj after = %lf\n",IndVolume[j]);
    fprintf(fp1,"\n");

}

void Volume::TriangleProperties(std::vector<Vertex> &vertexList,int i,int j,int k,FILE *fp1)
{
    double aux_i,aux_j,aux_k,vol_i,vol_j,vol_k;

    aux_i = met->Cap2Area (vertexList,i,j,k);
    aux_j = met->Cap2Area (vertexList,j,i,k);
    aux_k = met->Cap2Area (vertexList,k,i,j);

    vol_i = met->Cap2Volume (vertexList,i,j,k);
    vol_j = met->Cap2Volume (vertexList,j,i,k);
    vol_k = met->Cap2Volume (vertexList,k,i,j);

    fprintf(fp1,"surfi before = %lf\n",IndSurfArea[i]);
    fprintf(fp1,"surfj before = %lf\n",IndSurfArea[j]);
    fprintf(fp1,"surfk before = %lf\n",IndSurfArea[k]);
    fprintf(fp1,"voli before = %lf\n",IndVolume[i]);
    fprintf(fp1,"volj before = %lf\n",IndVolume[j]);
    fprintf(fp1,"volk before = %lf\n",IndVolume[k]);

    IndSurfArea[i] += aux_i;
    IndSurfArea[j] += aux_j;
    IndSurfArea[k] += aux_k;

    IndVolume[i] += vol_i;
    IndVolume[j] += vol_j;
    IndVolume[k] += vol_k;

    fprintf(fp1,"------------\n");

    fprintf(fp1,"surfi after = %lf\n",IndSurfArea[i]);
    fprintf(fp1,"surfj after = %lf\n",IndSurfArea[j]);
    fprintf(fp1,"surfk after = %lf\n",IndSurfArea[k]);
    fprintf(fp1,"voli after = %lf\n",IndVolume[i]);
    fprintf(fp1,"volj after = %lf\n",IndVolume[j]);
    fprintf(fp1,"volk after = %lf\n",IndVolume[k]);
    fprintf(fp1,"\n");
}

void Volume::TetProperties (std::vector<Vertex> &vertexList, int i, int j, int k, int l,FILE *fp1)
{
    double aux_i,aux_j,aux_k,aux_l,vol_i,vol_j,vol_k,vol_l;

    aux_i = met->Cap3Area (vertexList,i,j,k,l);
    aux_j = met->Cap3Area (vertexList,j,i,k,l);
    aux_k = met->Cap3Area (vertexList,k,i,j,l);
    aux_l = met->Cap3Area (vertexList,l,i,j,k);

    vol_i = met->Cap3Volume (vertexList,i,j,k,l);
    vol_j = met->Cap3Volume (vertexList,j,i,k,l);
    vol_k = met->Cap3Volume (vertexList,k,i,j,l);
    vol_l = met->Cap3Volume (vertexList,l,i,j,k);

    fprintf(fp1,"surfi before = %lf\n",IndSurfArea[i]);
    fprintf(fp1,"surfj before = %lf\n",IndSurfArea[j]);
    fprintf(fp1,"surfk before = %lf\n",IndSurfArea[k]);
    fprintf(fp1,"surfl before = %lf\n",IndSurfArea[l]);
    fprintf(fp1,"voli before = %lf\n",IndVolume[i]);
    fprintf(fp1,"volj before = %lf\n",IndVolume[j]);
    fprintf(fp1,"volk before = %lf\n",IndVolume[k]);
    fprintf(fp1,"voll before = %lf\n",IndVolume[l]);

    IndSurfArea[i] -= aux_i;
    IndSurfArea[j] -= aux_j;
    IndSurfArea[k] -= aux_k;
    IndSurfArea[l] -= aux_l;

    IndVolume[i] -= vol_i;
    IndVolume[j] -= vol_j;
    IndVolume[k] -= vol_k;
    IndVolume[l] -= vol_l;

    fprintf(fp1,"------------\n");

    fprintf(fp1,"surfi after = %lf\n",IndSurfArea[i]);
    fprintf(fp1,"surfj after = %lf\n",IndSurfArea[j]);
    fprintf(fp1,"surfk after = %lf\n",IndSurfArea[k]);
    fprintf(fp1,"surfl after = %lf\n",IndSurfArea[l]);
    fprintf(fp1,"voli after = %lf\n",IndVolume[i]);
    fprintf(fp1,"volj after = %lf\n",IndVolume[j]);
    fprintf(fp1,"volk after = %lf\n",IndVolume[k]);
    fprintf(fp1,"voll after = %lf\n",IndVolume[l]);
    fprintf(fp1,"\n");
}

void Volume::FindTotal(std::vector<Vertex> &vertexList,double *TotVol, double *TotSurf)
{
    FILE * fp = fopen("IndivVolume.txt","w");
    for (unsigned int i = 1; i < vertexList.size(); i++)
    {
        fprintf(fp,"IndSurfArea[%d] = %lf\n",i,IndSurfArea[i]);
        fprintf(fp,"IndVolume[%d] = %lf\n",i,IndVolume[i]);
        *TotSurf = *TotSurf + IndSurfArea[i];
        *TotVol = *TotVol + IndVolume[i];
     }
    fclose(fp);
}

void Volume::DoTetraVertex (std::vector<Vertex> &vertexList, int i, int j, int k, int l, int m,double *surf,double *vol)
{
    double sector_a, sector_v;

    sector_a = met->SectorArea (vertexList,i, j, k, l);
    sector_v = (1.0/3.0) * sector_a * vertexList[i].Radius;

    *surf = sector_a;
    *vol = sector_v;
}

void Volume::DoTetraEdge (std::vector<Vertex> &vertexList, int i, int j, int k, int l, int m,double *surf,double *vol)
{

    double angle_dih, auxi, auxj;
    double wedge_vol, cap_vol_ij, cap_vol_ji;

    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);
    Vector3 v(vertexList[l].Coordinates[1],vertexList[l].Coordinates[2],vertexList[l].Coordinates[3]);

    angle_dih = met->AngleDihedral (s,t,u,v);
    cap_vol_ij = met->CapVolume (vertexList,i,j);
    cap_vol_ji = met->CapVolume (vertexList,j,i);
    wedge_vol = angle_dih * (cap_vol_ij + cap_vol_ji);
    auxi = angle_dih * met->CapArea (vertexList,i,j);
    auxj = angle_dih * met->CapArea (vertexList,j,i);

    *surf = auxi + auxj;
    *vol = wedge_vol;
}

void Volume::DoTetraTriangle (std::vector<Vertex> &vertexList, int i, int j, int k, int l, int m,double *surf,double *vol)
{
    double auxi, auxj, auxk;
    double pawn_vol;
    double cap2_vol_i_jk, cap2_vol_j_ik, cap2_vol_k_ij;

    auxi = (1.0/2.0) * met->Cap2Area (vertexList,i,j,k);
    auxj = (1.0/2.0) * met->Cap2Area (vertexList,j,i,k);
    auxk = (1.0/2.0) * met->Cap2Area (vertexList,k,i,j);

    cap2_vol_i_jk = met->Cap2Volume (vertexList,i,j,k);
    cap2_vol_j_ik = met->Cap2Volume (vertexList,j,i,k);
    cap2_vol_k_ij = met->Cap2Volume (vertexList,k,i,j);

    pawn_vol = (cap2_vol_i_jk + cap2_vol_j_ik + cap2_vol_k_ij)/2.0;

    *surf = auxi+auxj+auxk;
    *vol = pawn_vol;
}

double Volume::DoTetraVolume (std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    return met->TetraVolume (vertexList[i].Coordinates,vertexList[j].Coordinates,vertexList[k].Coordinates,vertexList[l].Coordinates);
}

/*!
    \fn Volume::MeasureVolume(DelaunayComplex *delcx,std::vector<Vertex> & vertexList,int option)
 */
void Volume::MeasureVolume(DeluanayComplex *delcx,std::vector<Vertex> & vertexList,QLineEdit *totVol,QLineEdit *totSurf,int option)
{
        uint idx, i;
        int j, k, l, m;
        int pair1, pair2, pair3, pair4, pair5, pair6;
        int trig1, trig2, trig3, trig4;
        int flaga, flagb, flagc;
        int nlink, l1, l2;

        double distij = 0.0, distij2 = 0.0;
        double eps = 0.0, eps1 = 0.0, eps2 = 0.0, eps3 = 0.0, eps4 = 0.0;
        double ra, rb, rc, rd, ra2, rb2, rc2, rd2;
        double wa, wb, wc, wd, wl1, wl2;
        double surfa = 0.0, surfb = 0.0, surfc = 0.0, surfd = 0.0, vola = 0.0, volb = 0.0, volc = 0.0, vold = 0.0;
        double dist1 = 0.0, dist2 = 0.0, dist3 = 0.0, dist4 = 0.0, dist5 = 0.0, dist6 = 0.0;
        double d2_1 = 0.0, d2_2 = 0.0, d2_3 = 0.0, d2_4 = 0.0, d2_5 = 0.0, d2_6 = 0.0;
        double coefa, coefb, coefc, coefd;
        double sh_abc = 0.0, sh_acb = 0.0, sh_bca = 0.0, sh_abd = 0.0, sh_adb = 0.0, sh_bda = 0.0;
        double sh_acd = 0.0, sh_adc = 0.0, sh_cda = 0.0, sh_bcd = 0.0, sh_bdc = 0.0, sh_cdb = 0.0;

        double a[4];
        double b[4];
        double c[4];
        double d[4];

        double pabc[4];
        double pacb[4];
        double pabd[4];
        double padb[4];
        double padc[4];
        double pacd[4];
        double pbcd[4];
        double pbdc[4];

        double x[4];
        double y[4];

        double angles[4];

        double dsurfa3[4][3];
        double dsurfb3[4][3];

        double dsurfa4[4][4];
        double dsurfb4[4][4];
        double dsurfc4[4][4];

        double dsurfa5[4][5];
        double dsurfb5[4][5];
        double dsurfc5[4][5];
        double dsurfd5[4][5];

        double dvola3[4][3];
        double dvolb3[4][3];

        double dvola4[4][4];
        double dvolb4[4][4];
        double dvolc4[4][4];

        double dvolda4[4][4];
        double dvoldb4[4][4];
        double dvoldc4[4][4];

        double dvola5[4][5];
        double dvolb5[4][5];
        double dvolc5[4][5];
        double dvold5[4][5];

        double ang_abc[4];
        double ang_abd[4];
        double ang_acd[4];
        double ang_bcd[4];

        for (i = 1; i < vertexList.size(); i++)
        {
                for (j = 1; j <= 3; j++)
                {
                        DerivSurf[i].a[j] = 0;
                        DerivVol[i].a[j] = 0;
                }
                Radius2[i] = pow(vertexList[i].Radius, 2.0);
                CoefAsp[i] = 1.0;
        }

        spm->PrepareDeriv(delcx,nlink_trig,link_trig);


        for (i = 1; i < delcx->DeluanayEdges.size(); i++)
        {
                edgeCoef[i] = delcx->DeluanayEdges[i].AlphaStatus;
        }

        for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
        {
                if (delcx->DeluanayTrigs[i].AlphaStatus == 1)
                {
                        pair1 = delcx->DeluanayTrigs[i].TrigLink[1];
                        pair2 = delcx->DeluanayTrigs[i].TrigLink[2];
                        pair3 = delcx->DeluanayTrigs[i].TrigLink[3];
                        edgeCoef[pair1] = 0;
                        edgeCoef[pair2] = 0;
                        edgeCoef[pair3] = 0;
                }
        }

        for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
        {
                if (delcx->DeluanayTrigs[i].trigCoef != 0)
                {
                        pair1 = delcx->DeluanayTrigs[i].TrigLink[1];
                        pair2 = delcx->DeluanayTrigs[i].TrigLink[2];
                        pair3 = delcx->DeluanayTrigs[i].TrigLink[3];
                        edgeCoef[pair1] = 1;
                        edgeCoef[pair2] = 1;
                        edgeCoef[pair3] = 1;
                }
        }

        FILE *fp1 = fopen("status.txt","w");

        for (i = 1; i < vertexList.size(); i++)
        {
                fprintf(fp1,"%d %d\n",i,vertexList[i].AlphaStatus);
        }
        fprintf(fp1,"\n");
        for (i = 1; i < delcx->DeluanayEdges.size(); i++)
        {
                fprintf(fp1,"%d %d\n",i,delcx->DeluanayEdges[i].AlphaStatus);
        }
        fprintf(fp1,"\n");
        for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
        {
                fprintf(fp1,"%d %d\n",i,delcx->DeluanayTrigs[i].AlphaStatus);
        }
        fprintf(fp1,"\n");
        for (i = 1; i < delcx->DeluanayTet.size(); i++)
        {
                fprintf(fp1,"%d %d\n",i,delcx->DeluanayTet[i].AlphaStatus);
        }
        fclose(fp1);


        FILE *fp;
        fp = fopen("volume.txt","w");

        for (idx = 1; idx < vertexList.size(); idx++)
        {
                if (vertexList[idx].AlphaStatus != 0) continue;

                IndSurfArea[idx] = 4.0 * PI * Radius2[idx];
                IndVolume[idx] = IndSurfArea[idx] * vertexList[idx].Radius / 3.0;

                fprintf(fp,"point = %d\n",idx);
                fprintf(fp,"radius2 = %lf\n",Radius2[idx]);
                fprintf(fp,"surf =  %lf\n",IndSurfArea[idx]);
                fprintf(fp,"vol =  %lf\n",IndVolume[idx]);
                fprintf(fp,"\n");
        }

        for (idx = 1; idx < delcx->DeluanayEdges.size(); idx++)
        {
                if (delcx->DeluanayEdges[idx].AlphaStatus == 0) continue;

                fprintf(fp,"edge = %d\n",idx);

                i = delcx->DeluanayEdges[idx].Corners[1];
                j = delcx->DeluanayEdges[idx].Corners[2];

                fprintf(fp,"i = %d  | j = %d\n",i,j);

                if(i == 534 && j == 541)
                {
                    printf("found\n");
                }

                coefa = CoefAsp[i];
                coefb = CoefAsp[j];

                spm->Distance2(vertexList, i, j, &distij2);
                distij = sqrt(distij2);

                distpair[idx] = distij;
                distpair2[idx] = distij2;

                for (k = 1; k <= 3; k++)
                {
                        a[k] = vertexList[i].Coordinates[k];
                        b[k] = vertexList[j].Coordinates[k];
                }

                ra = vertexList[i].Radius;
                rb = vertexList[j].Radius;

                ra2 = Radius2[i];
                rb2 = Radius2[j];

                spm->TwoSphereVol(a, b, ra, ra2, rb, rb2, distij, distij2, &surfa, &surfb, &vola, &volb, dsurfa3, dsurfb3, dvola3, dvolb3, option);


                fprintf(fp,"surfi before = %lf\n",IndSurfArea[i]);
                fprintf(fp,"surfj before = %lf\n",IndSurfArea[j]);
                fprintf(fp,"voli before = %lf\n",IndVolume[i]);
                fprintf(fp,"volj before = %lf\n",IndVolume[j]);
                IndSurfArea[i] = IndSurfArea[i] - surfa;
                IndSurfArea[j] = IndSurfArea[j] - surfb;
                fprintf(fp,"------------\n");
                IndVolume[i] = IndVolume[i] - vola;
                IndVolume[j] = IndVolume[j] - volb;
                fprintf(fp,"surfi after = %lf\n",IndSurfArea[i]);
                fprintf(fp,"surfj after = %lf\n",IndSurfArea[j]);
                fprintf(fp,"voli after = %lf\n",IndVolume[i]);
                fprintf(fp,"volj after = %lf\n",IndVolume[j]);
                fprintf(fp,"\n");

                if (option == 0) continue;

                for (k = 1; k <= 3; k++)
                {
                        DerivSurf[i].a[k] = DerivSurf[i].a[k] - coefa * dsurfa3[k][1] - coefb * dsurfb3[k][1];
                        DerivSurf[j].a[k] = DerivSurf[j].a[k] - coefa * dsurfa3[k][2] - coefb * dsurfb3[k][2];

                        DerivVol[i].a[k] = DerivVol[i].a[k] - coefa * dvola3[k][1] - coefb * dvolb3[k][1];
                        DerivVol[j].a[k] = DerivVol[j].a[k] - coefa * dvola3[k][2] - coefb * dvolb3[k][2];
                }
        }

        for (idx = 1; idx < delcx->DeluanayTrigs.size(); idx++)
        {
                if (delcx->DeluanayTrigs[idx].AlphaStatus == 0) continue;

                if(idx == 2577)//53
                {
                        idx = 2577;
                }

                fprintf(fp,"trig = %d\n",idx);

                i = delcx->DeluanayTrigs[idx].Corners[1];
                j = delcx->DeluanayTrigs[idx].Corners[2];
                k = delcx->DeluanayTrigs[idx].Corners[3];

                fprintf(fp,"i = %d  | j = %d | k = %d\n",i,j,k);

                coefa = CoefAsp[i];
                coefb = CoefAsp[j];
                coefc = CoefAsp[k];

                for (m = 1; m <= 3; m++)
                {
                        a[m] = vertexList[i].Coordinates[m];
                        b[m] = vertexList[j].Coordinates[m];
                        c[m] = vertexList[k].Coordinates[m];
                }

                ra = vertexList[i].Radius;
                rb = vertexList[j].Radius;
                rc = vertexList[k].Radius;

                ra2 = Radius2[i];
                rb2 = Radius2[j];
                rc2 = Radius2[k];

                wa = 0.5 * vertexList[i].Weight;
                wb = 0.5 * vertexList[j].Weight;
                wc = 0.5 * vertexList[k].Weight;

                pair1 = delcx->DeluanayTrigs[idx].TrigLink[3];
                pair2 = delcx->DeluanayTrigs[idx].TrigLink[2];
                pair3 = delcx->DeluanayTrigs[idx].TrigLink[1];

                dist1 = distpair[pair1];
                dist2 = distpair[pair2];
                dist3 = distpair[pair3];

                d2_1 = distpair2[pair1];
                d2_2 = distpair2[pair2];
                d2_3 = distpair2[pair3];

                //ivol = 0;
                if (delcx->DeluanayTrigs[idx].trigCoef == 0)
                {
                        spm->ThreeVolDir(a, b, c, ra, rb, rc, ra2, rb2, rc2, wa, wb, wc, dist1, dist2, dist3, d2_1, d2_2, d2_3, dvola4, dvolb4, dvolc4, angles, pabc, pacb, &eps, &sh_abc, &sh_acb, &sh_bca, option);

                        fprintf(fp,"angles = %lf\n%lfn%lf\n",angles[1],angles[2],angles[3]);
                        fprintf(fp,"sh = %lf\n%lf\n%lf\n",sh_abc,sh_acb,sh_bca);

                        fprintf(fp,"pabc[1] = %lf\n",pabc[1]);
                        fprintf(fp,"pabc[2] = %lf\n",pabc[2]);
                        fprintf(fp,"pabc[3] = %lf\n",pabc[3]);

                        fprintf(fp,"pacb[1] = %lf\n",pacb[1]);
                        fprintf(fp,"pacb[2] = %lf\n",pacb[2]);
                        fprintf(fp,"pacb[3] = %lf\n",pacb[3]);
                        fprintf(fp,"\n");

                        trig_ang[idx].a[1] = angles[1];
                        trig_ang[idx].a[2] = angles[2];
                        trig_ang[idx].a[3] = angles[3];

                        trig_eps[idx] = eps;
                        trig_sh[idx].a[1] = sh_abc;
                        trig_sh[idx].a[2] = sh_acb;
                        trig_sh[idx].a[3] = sh_bca;

                        for (l = 1; l <= 3; l++)
                        {
                                trig_dual1[idx].a[l] = pabc[l];
                                trig_dual2[idx].a[l] = pacb[l];
                        }

                        if (option == 0) continue;

                        for (l = 1; l <= 3; l++)
                        {
                                DerivVol[i].a[l] = DerivVol[i].a[l] + coefa * dvola4[l][1] + coefb * dvolb4[l][1] + coefc * dvolc4[l][1];
                                DerivVol[j].a[l] = DerivVol[j].a[l] + coefa * dvola4[l][2] + coefb * dvolb4[l][2] + coefc * dvolc4[l][2];
                                DerivVol[k].a[l] = DerivVol[k].a[l] + coefa * dvola4[l][3] + coefb * dvolb4[l][3] + coefc * dvolc4[l][3];
                        }
                }
                else
                {
                        //ivol = 1;
                        spm->ThreeVolDist(a, b, c, ra, rb, rc, ra2, rb2, rc2, wa, wb, wc, dist1, dist2, dist3, d2_1, d2_2, d2_3, &surfa, &surfb, &surfc, &vola, &volb, &volc, dsurfa4, dsurfb4, dsurfc4, dvola4, dvolb4, dvolc4, dvolda4, dvoldb4, dvoldc4, &eps, pabc, pacb, angles, &sh_abc, &sh_acb, &sh_bca, option);

                        fprintf(fp,"angles = %lf\n%lf\n%lf\n",angles[1],angles[2],angles[3]);
                        fprintf(fp,"sh = %lf\n%lf\n%lf\n",sh_abc,sh_acb,sh_bca);

                        fprintf(fp,"pabc[1] = %lf\n",pabc[1]);
                        fprintf(fp,"pabc[2] = %lf\n",pabc[2]);
                        fprintf(fp,"pabc[3] = %lf\n",pabc[3]);

                        fprintf(fp,"pacb[1] = %lf\n",pacb[1]);
                        fprintf(fp,"pacb[2] = %lf\n",pacb[2]);
                        fprintf(fp,"pacb[3] = %lf\n",pacb[3]);

                        fprintf(fp,"surfi before = %lf\n",IndSurfArea[i]);
                        fprintf(fp,"surfj before = %lf\n",IndSurfArea[j]);
                        fprintf(fp,"surfk before = %lf\n",IndSurfArea[k]);
                        fprintf(fp,"voli before = %lf\n",IndVolume[i]);
                        fprintf(fp,"volj before = %lf\n",IndVolume[j]);
                        fprintf(fp,"volk before = %lf\n",IndVolume[k]);
                        IndSurfArea[i] = IndSurfArea[i] + 0.5 * surfa * delcx->DeluanayTrigs[idx].trigCoef;
                        IndSurfArea[j] = IndSurfArea[j] + 0.5 * surfb * delcx->DeluanayTrigs[idx].trigCoef;
                        IndSurfArea[k] = IndSurfArea[k] + 0.5 * surfc * delcx->DeluanayTrigs[idx].trigCoef;
                        fprintf(fp,"------------\n");
                        IndVolume[i] = IndVolume[i] + 0.5 * vola * delcx->DeluanayTrigs[idx].trigCoef;
                        IndVolume[j] = IndVolume[j] + 0.5 * volb * delcx->DeluanayTrigs[idx].trigCoef;
                        IndVolume[k] = IndVolume[k] + 0.5 * volc * delcx->DeluanayTrigs[idx].trigCoef;
                        fprintf(fp,"surfi after = %lf\n",IndSurfArea[i]);
                        fprintf(fp,"surfj after = %lf\n",IndSurfArea[j]);
                        fprintf(fp,"surfk after = %lf\n",IndSurfArea[k]);
                        fprintf(fp,"voli after = %lf\n",IndVolume[i]);
                        fprintf(fp,"volj after = %lf\n",IndVolume[j]);
                        fprintf(fp,"volk after = %lf\n",IndVolume[k]);
                        fprintf(fp,"\n");


                        trig_ang[idx].a[1] = angles[1];
                        trig_ang[idx].a[2] = angles[2];
                        trig_ang[idx].a[3] = angles[3];

                        trig_eps[idx] = eps;
                        trig_sh[idx].a[1] = sh_abc;
                        trig_sh[idx].a[2] = sh_acb;
                        trig_sh[idx].a[3] = sh_bca;


                        for (l = 1; l <= 3; l++)
                        {
                                trig_dual1[idx].a[l] = pabc[l];
                                trig_dual2[idx].a[l] = pacb[l];
                        }

                        if (option == 0) continue;

                        for (l = 1; l <= 3; l++)
                        {
                                DerivSurf[i].a[l] = DerivSurf[i].a[l] + 0.5 * delcx->DeluanayTrigs[idx].trigCoef * (coefa * dsurfa4[l][1] + coefb * dsurfb4[l][1] + coefc * dsurfc4[l][1]);

                                DerivSurf[j].a[l] = DerivSurf[j].a[l] + 0.5 * delcx->DeluanayTrigs[idx].trigCoef * (coefa * dsurfa4[l][2] + coefb * dsurfb4[l][2] + coefc * dsurfc4[l][2]);

                                DerivSurf[k].a[l] = DerivSurf[k].a[l] + 0.5 * delcx->DeluanayTrigs[idx].trigCoef * (coefa * dsurfa4[l][3] + coefb * dsurfb4[l][3] + coefc * dsurfc4[l][3]);

                                DerivVol[i].a[l] = DerivVol[i].a[l] + 0.5 * delcx->DeluanayTrigs[idx].trigCoef * (coefa * dvola4[l][1] + coefb * dvolb4[l][1] + coefc * dvolc4[l][1]);

                                DerivVol[j].a[l] = DerivVol[j].a[l] + 0.5 * delcx->DeluanayTrigs[idx].trigCoef * (coefa * dvola4[l][2] + coefb * dvolb4[l][2] + coefc * dvolc4[l][2]);

                                DerivVol[k].a[l] = DerivVol[k].a[l] + 0.5 * delcx->DeluanayTrigs[idx].trigCoef * (coefa * dvola4[l][3] + coefb * dvolb4[l][3] + coefc * dvolc4[l][3]);

                                DerivVol[i].a[l] = DerivVol[i].a[l] + coefa * dvolda4[l][1] + coefb * dvoldb4[l][1] + coefc * dvoldc4[l][1];

                                DerivVol[j].a[l] = DerivVol[j].a[l] + coefa * dvolda4[l][2] + coefb * dvoldb4[l][2] + coefc * dvoldc4[l][2];

                                DerivVol[k].a[l] = DerivVol[k].a[l] + coefa * dvolda4[l][3] + coefb * dvoldb4[l][3] + coefc * dvoldc4[l][3];
                        }
                }

                if ((edgeCoef[pair1] == 0) && (edgeCoef[pair2] == 0) && (edgeCoef[pair3] == 0)) continue;

                nlink = nlink_trig[idx];//delcx->DeluanayTrigs[idx].nLink;

                if (nlink == 0)
                {
                        l1 = i;
                        l2 = i;
                }
                else if (nlink == 1)
                {
                        l1 = link_trig[idx].l1;//delcx->DeluanayTrigs[idx].ReverseLink1;
                        l2 = i;
                }
                else
                {
                        l1 = link_trig[idx].l1;//delcx->DeluanayTrigs[idx].ReverseLink1;
                        l2 = link_trig[idx].l2;//delcx->DeluanayTrigs[idx].ReverseLink2;
                }

                wl1 = 0.5 * vertexList[l1].Weight;
                wl2 = 0.5 * vertexList[l2].Weight;

                flaga = edgeCoef[pair1];
                flagb = edgeCoef[pair2];
                flagc = edgeCoef[pair3];

                for (l = 1; l <= 3; l++)
                {
                        x[l] = vertexList[l1].Coordinates[l];
                        y[l] = vertexList[l2].Coordinates[l];
                }

                spm->ThreeSurfDir(a, b, c, x, y, nlink, pabc, pacb, eps, ra, rb, rc, ra2, wa, wb, wc, wl1, wl2, dist1, dist2, dist3, flaga, flagb, flagc, dsurfa4, dsurfb4, dsurfc4);

                for (l = 1; l <= 3; l++)
                {
                        DerivSurf[i].a[l] = DerivSurf[i].a[l] + coefa * dsurfa4[l][1] + coefb * dsurfb4[l][1] + coefc * dsurfc4[l][1];
                        DerivSurf[j].a[l] = DerivSurf[j].a[l] + coefa * dsurfa4[l][2] + coefb * dsurfb4[l][2] + coefc * dsurfc4[l][2];
                        DerivSurf[k].a[l] = DerivSurf[k].a[l] + coefa * dsurfa4[l][3] + coefb * dsurfb4[l][3] + coefc * dsurfc4[l][3];
                }
        }

        for (idx = 1; idx < delcx->DeluanayTet.size(); idx++)
        {
                fprintf(fp,"tet = %d\n",idx);
                fprintf(fp,"orient = %d\n",delcx->DeluanayTet[idx].Orient);

                if (delcx->DeluanayTet[idx].AlphaStatus != 1) continue;



                i = delcx->DeluanayTet[idx].Corners[1];
                j = delcx->DeluanayTet[idx].Corners[2];
                k = delcx->DeluanayTet[idx].Corners[3];
                l = delcx->DeluanayTet[idx].Corners[4];

                fprintf(fp,"i = %d  | j = %d | k = %d | l = %d\n",i,j,k,l);

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

                coefa = CoefAsp[i];
                coefb = CoefAsp[j];
                coefc = CoefAsp[k];
                coefd = CoefAsp[l];

                for (m = 1; m <= 3; m++)
                {
                        a[m] = vertexList[i].Coordinates[m];
                        b[m] = vertexList[j].Coordinates[m];
                        c[m] = vertexList[k].Coordinates[m];
                        d[m] = vertexList[l].Coordinates[m];
                }

                trig1 = delcx->DeluanayTet[idx].TetLink[4];
                trig2 = delcx->DeluanayTet[idx].TetLink[3];
                trig3 = delcx->DeluanayTet[idx].TetLink[2];
                trig4 = delcx->DeluanayTet[idx].TetLink[1];

                fprintf(fp,"trig1 = %d\n",trig1);
                fprintf(fp,"trig2 = %d\n",trig2);
                fprintf(fp,"trig3 = %d\n",trig3);
                fprintf(fp,"trig4 = %d\n",trig4);

                pair1 = delcx->DeluanayTrigs[trig1].TrigLink[3];
                pair2 = delcx->DeluanayTrigs[trig1].TrigLink[2];
                pair4 = delcx->DeluanayTrigs[trig1].TrigLink[1];
                pair3 = delcx->DeluanayTrigs[trig2].TrigLink[2];
                pair5 = delcx->DeluanayTrigs[trig2].TrigLink[1];
                pair6 = delcx->DeluanayTrigs[trig3].TrigLink[1];

                dist1 = distpair[pair1];
                dist2 = distpair[pair2];
                dist3 = distpair[pair3];
                dist4 = distpair[pair4];
                dist5 = distpair[pair5];
                dist6 = distpair[pair6];

                d2_1 = distpair2[pair1];
                d2_2 = distpair2[pair2];
                d2_3 = distpair2[pair3];
                d2_4 = distpair2[pair4];
                d2_5 = distpair2[pair5];
                d2_6 = distpair2[pair6];

                eps1 = trig_eps[trig1];
                eps2 = trig_eps[trig2];
                eps3 = trig_eps[trig3];
                eps4 = trig_eps[trig4];

                sh_abc = trig_sh[trig1].a[1];
                sh_acb = trig_sh[trig1].a[2];
                sh_bca = trig_sh[trig1].a[3];
                sh_abd = trig_sh[trig2].a[1];
                sh_adb = trig_sh[trig2].a[2];
                sh_bda = trig_sh[trig2].a[3];
                sh_acd = trig_sh[trig3].a[1];
                sh_adc = trig_sh[trig3].a[2];
                sh_cda = trig_sh[trig3].a[3];
                sh_bcd = trig_sh[trig4].a[1];
                sh_bdc = trig_sh[trig4].a[2];
                sh_cdb = trig_sh[trig4].a[3];

                for (m = 1; m <= 3; m++)
                {
                        pabc[m] = trig_dual1[trig1].a[m];
                        pacb[m] = trig_dual2[trig1].a[m];
                        pabd[m] = trig_dual1[trig2].a[m];
                        padb[m] = trig_dual2[trig2].a[m];
                        pacd[m] = trig_dual1[trig3].a[m];
                        padc[m] = trig_dual2[trig3].a[m];
                        pbcd[m] = trig_dual1[trig4].a[m];
                        pbdc[m] = trig_dual2[trig4].a[m];
                }

                ang_abc[1] = trig_ang[trig1].a[1];
                ang_abc[2] = trig_ang[trig1].a[2];
                ang_abc[3] = trig_ang[trig1].a[3];

                ang_abd[1] = trig_ang[trig2].a[1];
                ang_abd[2] = trig_ang[trig2].a[2];
                ang_abd[3] = trig_ang[trig2].a[3];

                ang_acd[1] = trig_ang[trig3].a[1];
                ang_acd[2] = trig_ang[trig3].a[2];
                ang_acd[3] = trig_ang[trig3].a[3];

                ang_bcd[1] = trig_ang[trig4].a[1];
                ang_bcd[2] = trig_ang[trig4].a[2];
                ang_bcd[3] = trig_ang[trig4].a[3];

                test = delcx->DeluanayTet[idx].Orient;

                if (test == 1)
                {
                        spm->FourSphereVol(i,j,k,l,a, b, c, d, ra, rb, rc, rd, ra2, rb2, rc2, rd2, dist1, dist2, dist3, dist4, dist5, dist6, d2_1, d2_2, d2_3, d2_4, d2_5, d2_6, wa, wb, wc, wd, eps1, eps2, eps3, eps4, sh_abc, sh_acb, sh_bca, sh_abd, sh_adb, sh_bda, sh_acd, sh_adc, sh_cda, sh_bcd, sh_bdc, sh_cdb, pacb, pabd, padc, pbcd, ang_abc, ang_abd, ang_acd, ang_bcd, &surfa, &surfb, &surfc, &surfd, &vola, &volb, &volc, &vold, dsurfa5, dsurfb5, dsurfc5, dsurfd5, dvola5, dvolb5, dvolc5, dvold5, option);
                }
                else
                {
                        ang_acd[1] = trig_ang[trig3].a[2];
                        ang_acd[2] = trig_ang[trig3].a[1];
                        ang_acd[3] = trig_ang[trig3].a[3];
                        ang_bcd[1] = trig_ang[trig4].a[2];
                        ang_bcd[2] = trig_ang[trig4].a[1];
                        ang_bcd[3] = trig_ang[trig4].a[3];

                        spm->FourSphereVol(i,j,k,l,a, b, d, c, ra, rb, rd, rc, ra2, rb2, rd2, rc2, dist1, dist3, dist2, dist5, dist4, dist6, d2_1, d2_3, d2_2, d2_5, d2_4, d2_6, wa, wb, wd, wc, eps2, eps1, eps3, eps4, sh_abd, sh_adb, sh_bda, sh_abc, sh_acb, sh_bca,sh_adc, sh_acd, sh_cda, sh_bdc, sh_bcd, sh_cdb, padb, pabc, pacd, pbdc, ang_abd, ang_abc, ang_acd, ang_bcd, &surfa, &surfb, &surfd, &surfc, &vola, &volb, &vold, &volc, dsurfa5, dsurfb5, dsurfd5, dsurfc5, dvola5, dvolb5, dvold5, dvolc5, option);
                }

                fprintf(fp,"test = %d\n",test);
                fprintf(fp,"surfi before = %lf\n",IndSurfArea[i]);
                fprintf(fp,"surfj before = %lf\n",IndSurfArea[j]);
                fprintf(fp,"surfk before = %lf\n",IndSurfArea[k]);
                fprintf(fp,"surfl before = %lf\n",IndSurfArea[l]);
                fprintf(fp,"voli before = %lf\n",IndVolume[i]);
                fprintf(fp,"volj before = %lf\n",IndVolume[j]);
                fprintf(fp,"volk before = %lf\n",IndVolume[k]);
                fprintf(fp,"voll before = %lf\n",IndVolume[l]);
                IndSurfArea[i] = IndSurfArea[i] - surfa;
                IndSurfArea[j] = IndSurfArea[j] - surfb;
                IndSurfArea[k] = IndSurfArea[k] - surfc;
                IndSurfArea[l] = IndSurfArea[l] - surfd;
                fprintf(fp,"------------\n");
                IndVolume[i] = IndVolume[i] - vola;
                IndVolume[j] = IndVolume[j] - volb;
                IndVolume[k] = IndVolume[k] - volc;
                IndVolume[l] = IndVolume[l] - vold;
                fprintf(fp,"surfi after = %lf\n",IndSurfArea[i]);
                fprintf(fp,"surfj after = %lf\n",IndSurfArea[j]);
                fprintf(fp,"surfk after = %lf\n",IndSurfArea[k]);
                fprintf(fp,"surfl after = %lf\n",IndSurfArea[l]);
                fprintf(fp,"voli after = %lf\n",IndVolume[i]);
                fprintf(fp,"volj after = %lf\n",IndVolume[j]);
                fprintf(fp,"volk after = %lf\n",IndVolume[k]);
                fprintf(fp,"voll after = %lf\n",IndVolume[l]);
                fprintf(fp,"\n");

                if (option == 0) continue;

                for (m = 1; m <= 3; m++)
                {
                        DerivSurf[i].a[m] = DerivSurf[i].a[m] - coefa * dsurfa5[m][1] - coefb * dsurfb5[m][1] - coefc * dsurfc5[m][1] - coefd * dsurfd5[m][1];
                        DerivSurf[j].a[m] = DerivSurf[j].a[m] - coefa * dsurfa5[m][2] - coefb * dsurfb5[m][2] - coefc * dsurfc5[m][2] - coefd * dsurfd5[m][2];

                        DerivVol[i].a[m] = DerivVol[i].a[m] - coefa * dvola5[m][1] - coefb * dvolb5[m][1] - coefc * dvolc5[m][1] - coefd * dvold5[m][1];
                        DerivVol[j].a[m] = DerivVol[j].a[m] - coefa * dvola5[m][2] - coefb * dvolb5[m][2] - coefc * dvolc5[m][2] - coefd * dvold5[m][2];

                        if (test == 1)
                        {
                                DerivSurf[k].a[m] = DerivSurf[k].a[m] - coefa * dsurfa5[m][3] - coefb * dsurfb5[m][3] - coefc * dsurfc5[m][3] - coefd * dsurfd5[m][3];
                                DerivSurf[l].a[m] = DerivSurf[l].a[m] - coefa * dsurfa5[m][4] - coefb * dsurfb5[m][4] - coefc * dsurfc5[m][4] - coefd * dsurfd5[m][4];

                                DerivVol[k].a[m] = DerivVol[k].a[m] - coefa * dvola5[m][3] - coefb * dvolb5[m][3] - coefc * dvolc5[m][3] - coefd * dvold5[m][3];
                                DerivVol[l].a[m] = DerivVol[l].a[m] - coefa * dvola5[m][4] - coefb * dvolb5[m][4] - coefc * dvolc5[m][4] - coefd * dvold5[m][4];
                        }
                        else
                        {
                                DerivSurf[k].a[m] = DerivSurf[k].a[m] - coefa * dsurfa5[m][4] - coefb * dsurfb5[m][4] - coefc * dsurfc5[m][4] - coefd * dsurfd5[m][4];
                                DerivSurf[l].a[m] = DerivSurf[l].a[m] - coefa * dsurfa5[m][3] - coefb * dsurfb5[m][3] - coefc * dsurfc5[m][3] - coefd * dsurfd5[m][3];

                                DerivVol[k].a[m] = DerivVol[k].a[m] - coefa * dvola5[m][4] - coefb * dvolb5[m][4] - coefc * dvolc5[m][4] - coefd * dvold5[m][4];
                                DerivVol[l].a[m] = DerivVol[l].a[m] - coefa * dvola5[m][3] - coefb * dvolb5[m][3] - coefc * dvolc5[m][3] - coefd * dvold5[m][3];
                        }
                }
        }
        fclose(fp);

        TotSurfaceArea = 0.0;
        TotVolume = 0.0;
        Ssolv = 0.0;
        Vsolv = 0.0;

        FILE * fp2 = fopen("anothervolume.txt","w");

        for (i = 1; i < vertexList.size(); i++)
        {
                fprintf(fp2,"surf[%d] = %lf\n",i,IndSurfArea[i]);
                fprintf(fp2,"vol[%d] = %lf\n",i,IndVolume[i]);
                fprintf(fp2,"\n");
                TotSurfaceArea = TotSurfaceArea + IndSurfArea[i];
                TotVolume = TotVolume + IndVolume[i];
                Ssolv = Ssolv + CoefAsp[i] * IndSurfArea[i];
                Vsolv = Vsolv + CoefAsp[i] * IndVolume[i];
        }

        fprintf(fp2,"total surf = %lf\n",TotSurfaceArea);
        fprintf(fp2,"total vol = %lf\n",TotVolume);

        fclose(fp2);
        QString insertionString;
        QString &ran = insertionString.setNum(TotVolume,'f',5);
        totVol->setText(ran);

        ran = insertionString.setNum(TotSurfaceArea,'f',5);
        totSurf->setText(ran);
}
