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

#include "deluanaycomplex.h"

int DeluanayComplex::inf4_1[5] =	{ 0, 2, 2, 1, 1 };

int DeluanayComplex::sign4_1[5] =	{ 0, -1, 1, 1, -1 };

int DeluanayComplex::inf4_2[5][5] =	{
                                        { 0, 0, 0, 0, 0 },
                                        { 0, 0, 2, 3, 3 },
                                        { 0, 2, 0, 3, 3 },
                                        { 0, 3, 3, 0, 1 },
                                        { 0, 3, 3, 1, 0 }
                                        };

int DeluanayComplex::sign4_2[5][5] =	{
                                        { 0, 0, 0,  0, 0 },
                                        { 0, 0, 1, -1, 1 },
                                        { 0,-1, 0, 1, -1 },
                                        { 0, 1, -1, 0, 1 },
                                        { 0,-1, 1, -1, 0 }
                                        };

int DeluanayComplex::sign4_3[5] =	{ 0, -1, 1, -1, 1 };

int DeluanayComplex::inf5_2[5][5] =	{
                                        { 0, 0, 0, 0, 0 },
                                        { 0, 0, 2, 1, 1 },
                                        { 0, 2, 0, 1, 1 },
                                        { 0, 1, 1, 0, 1 },
                                        { 0, 1, 1, 1, 0 }
                                        };

int DeluanayComplex::sign5_2[5][5] =	{
                                        { 0, 0, 0, 0, 0 },
                                        { 0, 0, -1, -1, 1 },
                                        { 0, 1, 0, -1, 1 },
                                        { 0, 1, 1, 0, 1 },
                                        { 0, -1, -1, -1, 0 }
                                        };

int DeluanayComplex::inf5_3[5] =	{ 0, 1, 1, 3, 3 };

int DeluanayComplex::sign5_3[5] =	{ 0, 1, 1, -1, 1 };

int DeluanayComplex::order1[5][4] =	{
                                        { 0, 0, 0, 0 },
                                        { 0, 3, 2, 4 },
                                        { 0, 1, 3, 4 },
                                        { 0, 2, 1, 4 },
                                        { 0, 1, 2, 3 }
                                        };

int DeluanayComplex::order2[7][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 3, 4 },
                                        { 0, 4, 2 },
                                        { 0, 2, 3 },
                                        { 0, 1, 4 },
                                        { 0, 3, 1 },
                                        { 0, 1, 2 }
                                        };

int DeluanayComplex::order3[7][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 1, 2 },
                                        { 0, 1, 3 },
                                        { 0, 1, 4 },
                                        { 0, 2, 3 },
                                        { 0, 2, 4 },
                                        { 0, 3, 4 }
                                        };

int DeluanayComplex::idxList[5][4] =	{
                                        { 0, 0, 0, 0 },
                                        { 0, 1, 1, 1 },
                                        { 0, 1, 2, 2 },
                                        { 0, 2, 2, 3 },
                                        { 0, 3, 3, 3 }
                                        };

int DeluanayComplex::table32[4][4] =	{
                                        { 0, 0, 0, 0 },
                                        { 0, 1, 2, 3 },
                                        { 0, 1, 3, 2 },
                                        { 0, 3, 1, 2 }
                                        };

int DeluanayComplex::table32_2[4][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 1, 2 },
                                        { 0, 1, 3 },
                                        { 0, 2, 3 }
                                        };

int DeluanayComplex::table41[4][4] =	{
                                        { 0, 0, 0, 0 },
                                        { 0, 2, 1, 3 },
                                        { 0, 1, 2, 3 },
                                        { 0, 1, 3, 2 }
                                        };

int DeluanayComplex::table41_2[4][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 1, 1 },
                                        { 0, 2, 1 },
                                        { 0, 2, 2 }
                                        };

int DeluanayComplex::order[4][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 2, 3 },
                                        { 0, 3, 1 },
                                        { 0, 1, 2 }
                                        };

int DeluanayComplex::order_1[4][4] =	{
                                        { 0, 0, 0, 0 },
                                        { 0, 1, 2, 3 },
                                        { 0, 3, 1, 2 },
                                        { 0, 2, 3, 1 }
                                        };

int DeluanayComplex::other[5][4] =	{
                                        { 0, 0, 0, 0 },
                                        { 0, 2, 3, 4 },
                                        { 0, 1, 3, 4 },
                                        { 0, 1, 2, 4 },
                                        { 0, 1, 2, 3 }
                                        };

int DeluanayComplex::idx_list[4][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 1, 1 },
                                        { 0, 1, 2 },
                                        { 0, 2, 2 }
                                        };

int DeluanayComplex::other2[5][5][3] =	{
                                                {
                                                        { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }
                                                },
                                                {
                                                        { 0, 0, 0 }, { 0, 0, 0 }, { 0, 3, 4 }, { 0, 2, 4 }, { 0, 2, 3 }
                                                },
                                                {
                                                        { 0, 0, 0 }, { 0, 3, 4 }, { 0, 0, 0 }, { 0, 1, 4 }, { 0, 1, 3 }
                                                },
                                                {
                                                        { 0, 0, 0 }, { 0, 2, 4 }, { 0, 1, 4 }, { 0, 0, 0 }, { 0, 1, 2 }
                                                },
                                                {
                                                        { 0, 0, 0 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 }, { 0, 0, 0 }
                                                },
                                        };

int DeluanayComplex::mask[7] =		{ 0, 2, 4, 8, 16, 32, 64 };

int DeluanayComplex::pair[7][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 3, 4 },
                                        { 0, 2, 4 },
                                        { 0, 2, 3 },
                                        { 0, 1, 4 },
                                        { 0, 1, 3 },
                                        { 0, 1, 2 }
                                        };

int DeluanayComplex::trig_info[7][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 1, 2 },
                                        { 0, 1, 3 },
                                        { 0, 1, 4 },
                                        { 0, 2, 3 },
                                        { 0, 2, 4 },
                                        { 0, 3, 4 }
                                        };

int DeluanayComplex::trig_pos[7][3] =	{
                                        { 0, 0, 0 },
                                        { 0, 1, 1 },
                                        { 0, 2, 1 },
                                        { 0, 3, 1 },
                                        { 0, 2, 2 },
                                        { 0, 3, 2 },
                                        { 0, 3, 3 }
                                        };

double DeluanayComplex::eps = 0.01;
int DeluanayComplex::ir[98];
int DeluanayComplex::iy;

int Nvert = 0;
int NTri = 0;
std::vector<int> mark_buf;
std::vector<int> tri_buf;
int tcount = 0;

static double myMax(double a, double b)
{
        if (b > a) return b;
        return a;
}

static double myMin(double a, double b)
{
        if (b < a) return b;
        return a;
}

DeluanayComplex::DeluanayComplex()
{
        tetraLoc = -1;
        iredundant = 0;
        redundantCount = 0;
        nfree = 0;
        iff = 0;
        sos = new Sos();
        DeluanayTet.clear();
}

DeluanayComplex::~DeluanayComplex()
{
        DeluanayTet.clear();
        delete sos;
}

/*!
    \fn DeluanayComplex::ValSort2(int a,int b,int *ia,int *ib,int *iswap)
 */
void DeluanayComplex::ValSort2(int a,int b,int *ia,int *ib,int *iswap)
{
        *iswap = 1;
        if (a > b)
        {
                *ia = b;
                *ib = a;
                *iswap = -*iswap;
        }
        else
        {
                *ia = a;
                *ib = b;
        }
}

/*!
    \fn DeluanayComplex::ValSort3(int a,int b,int c,int *ia,int *ib,int *ic,int *iswap)
 */
void DeluanayComplex::ValSort3(int a,int b,int c,int *ia,int *ib,int *ic,int *iswap)
{
        int temp;

        ValSort2(a, b, ia, ib, iswap);

        *ic = c;

        if (*ib > *ic)
        {
                temp = *ib;
                *ib = *ic;
                *ic = temp;
                *iswap = -*iswap;
                if (*ia > *ib)
                {
                        temp = *ia;
                        *ia = *ib;
                        *ib = temp;
                        *iswap = -*iswap;
                }
        }
}

/*!
    \fn DeluanayComplex::ValSort4(int a,int b,int c,int d,int *ia,int *ib,int *ic,int *id,int *iswap)
 */
void DeluanayComplex::ValSort4(int a,int b,int c,int d,int *ia,int *ib,int *ic,int *id,int *iswap)
{
        int temp;

        ValSort3(a, b, c, ia, ib, ic, iswap);

        *id = d;

        if(*ic > *id)
        {
                temp = *ic;
                *ic = *id;
                *id = temp;
                *iswap = -*iswap;
                if(*ib > *ic)
                {
                        temp = *ib;
                        *ib = *ic;
                        *ic = temp;
                        *iswap = -*iswap;
                        if(*ia > *ib)
                        {
                                temp = *ia;
                                *ia = *ib;
                                *ib = temp;
                                *iswap = -*iswap;
                        }
                }
        }
}

/*!
    \fn DeluanayComplex::ValSort5(int a,int b,int c,int d,int e,int *ia,int *ib,int *ic,int *id,int *ie,int *iswap)
 */
void DeluanayComplex::ValSort5(int a,int b,int c,int d,int e,int *ia,int *ib,int *ic,int *id,int *ie,int *iswap)
{
        int temp;

        ValSort4(a, b, c, d, ia, ib, ic, id, iswap);

        *ie = e;

        if (*id > *ie)
        {
                temp = *id;
                *id = *ie;
                *ie = temp;
                *iswap = -*iswap;
                if (*ic > *id)
                {
                        temp = *ic;
                        *ic = *id;
                        *id = temp;
                        *iswap = -*iswap;
                        if (*ib > *ic)
                        {
                                temp = *ib;
                                *ib = *ic;
                                *ic = temp;
                                *iswap = -*iswap;
                                if (*ia > *ib)
                                {
                                        temp = *ia;
                                        *ia = *ib;
                                        *ib = temp;
                                        *iswap = -*iswap;
                                }
                        }
                }
        }
}

/*!
    \fn DeluanayComplex::EliminateDup(std::vector <Vertex> &sortedVertexList)
 */
void DeluanayComplex::EliminateDup(std::vector <Vertex> &sortedVertexList)
{
        uint ip, jp;
        double x, y, z, w;
        double xi, yi, zi, wi;

        jp = 1;
        x = sortedVertexList[jp].Coordinates[0];
        y = sortedVertexList[jp].Coordinates[1];
        z = sortedVertexList[jp].Coordinates[2];
        w = sortedVertexList[jp].Radius;

        for (uint i = 2; i < sortedVertexList.size(); i++)
        {
                ip = i;
                xi = sortedVertexList[ip].Coordinates[0];
                yi = sortedVertexList[ip].Coordinates[1];
                zi = sortedVertexList[ip].Coordinates[2];
                wi = sortedVertexList[ip].Radius;

                if ((xi == x) && (yi == y) && (zi == z))
                {
                        if (wi <= w)
                        {
                                sortedVertexList[ip].Redinfo = 1;
                        }
                        else
                        {
                                sortedVertexList[jp].Redinfo = 1;
                                jp = ip;
                                w = wi;
                        }
                }
                else
                {
                        x = xi;
                        y = yi;
                        z = zi;
                        w = wi;
                }
        }
}

/*!
    \fn DeluanayComplex::AddDummyLF()
 */
void DeluanayComplex::AddDummyLF()
{
        assert((linkFacet.size() == 0) && (linkIndex.size() == 0));
        LinkFacet lf;
        LinkIndex li;
        lf.f1 = 0;
        lf.f2 = 0;
        li.i1 = 0;
        li.i2 = 0;
        linkFacet.push_back(lf);
        linkIndex.push_back(li);
}

/*!
    \fn DeluanayComplex::ran2(int *idum)
 */
double DeluanayComplex::ran2(int *idum)
{
        int m, ia, ic, j;
        double rm;
        double returnVal;
        m = 714025;
        ia = 1366;
        ic = 150889;
        rm = 1.0 / m;

        if ((*idum < 0) || (iff == 0))
        {
                iff = 1;
                *idum = (ic - *idum) % m;
                //init shuffle table
                for (j = 1; j <= 97; j++)
                {
                        *idum = (ia * *idum + ic) % m;
                        ir[j] = *idum;
                }
                *idum = (ia * (*idum) + ic) % m;
                iy = *idum;
        }
        j = 1 + (97 * iy) / m;
        assert((j>=1) && (j<=97));
        iy = ir[j];

        returnVal = iy * rm;
        *idum = (ia * (*idum) + ic) % m;
        ir[j] = *idum;

        return (returnVal);
}

/*!
    \fn DeluanayComplex::compare(Vertex *a,Vertex *b)
 */
int DeluanayComplex::compare(Vertex *a,Vertex *b)
{
        int ret;
        if(a->ranValue > b->ranValue)
        {
                ret = 1;
        }
        else if(a->ranValue < b->ranValue)
        {
                ret = -1;
        }
        else
        {
                ret = 0;
        }
        return(ret);
}

/*!
    \fn DeluanayComplex::ISortIdx(int vertice[],int idx[],int *nswap,int *n)
 */
void DeluanayComplex::ISortIdx(int vertice[],int idx[],int *nswap,int *n)
{
        int i, j, a;

        for (i = 1; i <= *n; i++)
        {
                idx[i] = i;
        }

        *nswap = 0;

        for (i = 1; i < *n; i++)
        {
                for (j = i + 1; j <= *n; j++)
                {
                        if (vertice[i] > vertice[j])
                        {
                                a = vertice[i];
                                vertice[i] = vertice[j];
                                vertice[j] = a;
                                a = idx[i];
                                idx[i] = idx[j];
                                idx[j] = a;
                                (*nswap)++;         //the most stupid mistake committed was to do *nswap++ :(
                        }
                }
        }
}

/*!
    \fn DeluanayComplex::MissInfSign(int i,int j,int k,int *l,int *sign)
 */
void DeluanayComplex::MissInfSign(int i,int j,int k,int *l,int *sign)
{
        int a, b, c, d;

        *l = 10 - i - j - k;

        a = i;
        b = j;
        c = k;

        *sign = 1;

        if (a > b)
        {
                d = a;
                a = b;
                b = d;
                *sign = -*sign;
        }

        if (a > c)
        {
                d = a;
                a = c;
                c = d;
                *sign = -*sign;
        }

        if (b > c)
        {
                *sign = -*sign;
        }
}

/*!
    \fn DeluanayComplex::TetraVol(std::vector<Vertex>& vertexList, int ia, int ib, int ic, int id, double *vol)
 */
void DeluanayComplex::TetraVol(std::vector<Vertex>& vertexList, int ia, int ib, int ic, int id, double *vol)
{
        int i;
        double ad[4];
        double bd[4];
        double cd[4];
        double Sbcd[4];
        /*
        The volume of the tetrahedron is proportional to:

        vol =det| a(1)  a(2)  a(3)  1|
                | b(1)  b(2)  b(3)  1|
                | c(1)  c(2)  c(3)  1|
                | d(1)  d(2)  d(3)  1|

        substract last row from the first 3 rows, and expand w.r.t last column

        vol =det| ad(1)  ad(2)  ad(3) |
                | bd(1)  bd(2)  bd(3) |
                | cd(1)  cd(2)  cd(3) |

        where ad(i) = a(i) - d(i), ...
        */
        for (i = 1; i <= 3; i++)
        {
                ad[i] = vertexList[ia].Coordinates[i] - vertexList[id].Coordinates[i];
                bd[i] = vertexList[ib].Coordinates[i] - vertexList[id].Coordinates[i];
                cd[i] = vertexList[ic].Coordinates[i] - vertexList[id].Coordinates[i];
        }

        Sbcd[3] = bd[1] * cd[2] - cd[1] * bd[2];
        Sbcd[2] = bd[1] * cd[3] - cd[1] * bd[3];
        Sbcd[1] = bd[2] * cd[3] - cd[2] * bd[3];

        *vol = ad[1] * Sbcd[1] - ad[2] * Sbcd[2] + ad[3] * Sbcd[3];
}

/*!
    \fn DeluanayComplex:: Peel(std::vector<Vertex> & vertexList)
 */
void DeluanayComplex::Peel(std::vector<Vertex> & vertexList)
{
        uint i;
        int j, k, l;
        int ia, ib, ic, id, val = 0;
        double vol = 0.0;
        bool flag = false;

        for (i = 1; i < DeluanayTet.size(); i++)
        {
                //ignore dead tets
                if (DeluanayTet[i].Status == 0) continue;

                for (j = 1; j <= 4; j++)
                {
                        if (DeluanayTet[i].Neighbours[j] == 0)
                        {
                                flag = true;
                                break;
                        }
                }

                if (!flag) continue;       //tet is interior hence cant be flat
                else
                {
                        //tet is at boundary test if its flat be see if vol == 0
                        ia = DeluanayTet[i].Corners[1];
                        ib = DeluanayTet[i].Corners[2];
                        ic = DeluanayTet[i].Corners[3];
                        id = DeluanayTet[i].Corners[4];

                        TetraVol(vertexList, ia, ib, ic, id, &vol);

                        if (fabs(vol) < eps)
                        {
                                sos->Minor4(vertexList, &ia, &ib, &ic, &id, &val);
                                if (val == 0) DeluanayTet[i].Status = -1;
                        }
                }
        }

        //now remove flat tets and update its neighbours
        for (i = 1; i < DeluanayTet.size(); i++)
        {
                if (DeluanayTet[i].Status == -1)
                {
                        DeluanayTet[i].Status = 0;
                        for (j = 1; j <= 4; j++)
                        {
                                k = DeluanayTet[i].Neighbours[j];
                                if (k != 0)
                                {
                                        l = DeluanayTet[i].Nindex[j];
                                        DeluanayTet[k].Neighbours[l] = 0;
                                }
                        }
                }
        }
}

/*!
    \fn DeluanayComplex::MarkZero(int itetra, int ivertex)
 */
void DeluanayComplex::MarkZero(int itetra, int ivertex)
{
        int jtetra, jvertex;

        jtetra = DeluanayTet[itetra].Neighbours[ivertex];
        jvertex = DeluanayTet[itetra].Nindex[ivertex];

        if (jtetra != 0)
        {
                DeluanayTet[jtetra].Neighbours[jvertex] = 0;
        }
}

/*!
    \fn DeluanayComplex::RemoveInf()
 */
void DeluanayComplex::RemoveInf()
{
        uint i;
        int a, b, c, d;
        int ninf;

        for (i = 1; i < DeluanayTet.size(); i++)
        {
                //ignore dead tet
                if (DeluanayTet[i].Status == 0) continue;

                a = DeluanayTet[i].Corners[1];
                b = DeluanayTet[i].Corners[2];
                c = DeluanayTet[i].Corners[3];
                d = DeluanayTet[i].Corners[4];

                ninf = infPoint[a] + infPoint[b] + infPoint[c] + infPoint[d];

                if (ninf != 0)
                {
                        DeluanayTet[i].Status = 0;
                        if (infPoint[a] == 1) MarkZero(i, 1);
                        if (infPoint[b] == 1) MarkZero(i, 2);
                        if (infPoint[c] == 1) MarkZero(i, 3);
                        if (infPoint[d] == 1) MarkZero(i, 4);
                }
        }
}

/*!
    \fn DeluanayComplex::RegularConvex(std::vector<Vertex> & vertexList,int a, int b, int c, int p, int o, int itest_abcp, bool *regular, bool *convex, bool *test_abpo, bool *test_bcpo, bool *test_capo)
 */
void DeluanayComplex::RegularConvex(std::vector<Vertex> & vertexList,int a, int b, int c, int p, int o, int itest_abcp, bool *regular, bool *convex, bool *test_abpo, bool *test_bcpo, bool *test_capo)
{
        int i,j,k,l,m;
        int ia = 0, ib = 0, ic = 0, id = 0, ie = 0;
        int npoints = vertexList.size() - 1;
        int ninf, infp, info, iswap = 0, iswap2 = 0, idx, val = 0;
        int icol1, sign1, icol2, sign2, icol4 = 0, sign4, icol5, sign5;
        int list[4];

        double det_abpo,det_bcpo,det_capo,det_abcpo,det_abpc;
        double a_p[5];
        double b_p[5];
        double c_p[5];
        double o_p[5];
        double i_p[4];
        double j_p[4];
        double Mbo[4];
        double Mca[4];
        double Mjo[4];
        double Mio[4];
        double coordp[4];

        bool testc[4];
        /*
        To test if the union of the two tetrahedron is convex, check the position of o w.r.t
        three faces (a,b,p), (b,c,p) and (c,a,p) of (a,b,c,p).evaluate the three determinants:
        det(a,b,p,o)
        det(b,c,p,o)
        det(c,a,p,o)
        If the three determinants are positive, & det(a,b,c,p) is negative,the union is convex
        If the three determinants are negative, & det(a,b,c,p) is positive,the union is convex
        In all other cases, the union is non convex
        The regularity is tested by computing det(a,b,c,p,o)

        first count how many infinite points (except o)
        only a and/or b and/or c can be infinite:
        */
        ninf = infPoint[a] + infPoint[b] + infPoint[c];

        list[1] = a;
        list[2] = b;
        list[3] = c;

        for (m = 1; m <= 3; m++)
        {
                coordp[m] = vertexList[p].Coordinates[m];
        }
        //general case:no inf points
        if (ninf == 0)
        {
                //First, simple case:
                //if o is infinite, then det(a,b,c,p,o) = -det(a,b,c,p)
                //and consequently (a,b,c,p,o) is regular:nothing to do!
                if(infPoint[o] == 1)
                {
                        *regular = true;
                        return;
                }
                /*
                The three determinants det(a,b,p,o), det(b,c,p,o), and det(c,a,p,o) are "real" 4x4 determinants.
                Subtract the row corresponding to p from the other row, and develop with respect to p.
                The determinants become:
                det(a,b,p,o)= - | ap(1) ap(2) ap(3) |
                                | bp(1) bp(2) bp(3) |
                                | op(1) op(2) op(3) |

                det(b,c,p,o)=  -| bp(1) bp(2) bp(3) |
                                | cp(1) cp(2) cp(3) |
                                | op(1) op(2) op(3) |

                det(c,a,p,o)=  -| cp(1) cp(2) cp(3) |
                                | ap(1) ap(2) ap(3) |
                                | op(1) op(2) op(3) |

                Where ip(j) = i(j) - p(j) for all i in {a,b,c,o} and j in {1,2,3}
                        Compute two types of minors:

                Mbo_ij = bp(i)op(j) - bp(j)op(i) and Mca_ij = cp(i)ap(j) - cp(j)op(i)

                Store Mbo_12 in Mbo(3), Mbo_13 in Mbo(2),...
                */
                //get the coordinates
                for (m = 1; m <= 3; m++)
                {
                        a_p[m] = vertexList[a].Coordinates[m] - coordp[m];
                        b_p[m] = vertexList[b].Coordinates[m] - coordp[m];
                        c_p[m] = vertexList[c].Coordinates[m] - coordp[m];
                        o_p[m] = vertexList[o].Coordinates[m] - coordp[m];
                }
                a_p[4] = vertexList[a].Weight - vertexList[p].Weight;
                b_p[4] = vertexList[b].Weight - vertexList[p].Weight;
                c_p[4] = vertexList[c].Weight - vertexList[p].Weight;
                o_p[4] = vertexList[o].Weight - vertexList[p].Weight;

                //compute 2x2 determinants Mbo and Mca
                Mbo[1] = b_p[2] * o_p[3] - b_p[3] * o_p[2];
                Mbo[2] = b_p[1] * o_p[3] - b_p[3] * o_p[1];
                Mbo[3] = b_p[1] * o_p[2] - b_p[2] * o_p[1];

                Mca[1] = c_p[2] * a_p[3] - c_p[3] * a_p[2];
                Mca[2] = c_p[1] * a_p[3] - c_p[3] * a_p[1];
                Mca[3] = c_p[1] * a_p[2] - c_p[2] * a_p[1];

                det_abpo = -a_p[1] * Mbo[1] + a_p[2] * Mbo[2] - a_p[3] * Mbo[3];
                det_bcpo =  c_p[1] * Mbo[1] - c_p[2] * Mbo[2] + c_p[3] * Mbo[3];
                det_capo = -o_p[1] * Mca[1] + o_p[2] * Mca[2] - o_p[3] * Mca[3];
                det_abpc = -b_p[1] * Mca[1] + b_p[2] * Mca[2] - b_p[3] * Mca[3];
                /*
                To compute:
                det(a,b,c,p,o) =	| a(1) a(2) a(3) a(4) 1 |
                                        | b(1) b(2) b(3) b(4) 1 |
                                        | c(1) c(2) c(3) c(4) 1 |
                                        | p(1) p(2) p(3) p(4) 1 |
                                        | o(1) o(2) o(3) o(4) 1 |

                First substract row p :
                det(a,b,c,p,o) =      - | ap(1) ap(2) ap(3) ap(4) |
                                        | bp(1) bp(2) bp(3) bp(4) |
                                        | cp(1) cp(2) cp(3) cp(4) |
                                        | op(1) op(2) op(3) op(4) |

                Expand w.r.t last column
                */
                det_abcpo = -a_p[4] * det_bcpo - b_p[4] * det_capo - c_p[4] * det_abpo + o_p[4] * det_abpc;

                //if (a,b,c,p,o) regular, no need to flip

                if(fabs(det_abcpo) < eps)
                {
                        ValSort5(a, b, c, p, o, &ia, &ib, &ic, &id, &ie, &iswap);
                        sos->Minor5(vertexList, &ia, &ib, &ic, &id, &ie, &val);
                        det_abcpo = val * iswap;
                }

                if ((det_abcpo * itest_abcp) < 0)
                {
                        *regular = true;
                        return;
                }
                *regular = false;

                //If not regular, we test for convexity
                if(fabs(det_abpo) < eps)
                {
                        ValSort4(a, b, p, o, &ia, &ib, &ic, &id, &iswap);
                        sos->Minor4(vertexList, &ia, &ib, &ic, &id, &val);
                        det_abpo = val * iswap;
                }
                if(fabs(det_bcpo) < eps)
                {
                        ValSort4(b, c, p, o, &ia, &ib, &ic, &id, &iswap);
                        sos->Minor4(vertexList, &ia, &ib, &ic, &id, &val);
                        det_bcpo = val * iswap;
                }
                if(fabs(det_capo) < eps)
                {
                        ValSort4(c, a, p, o, &ia, &ib, &ic, &id, &iswap);
                        sos->Minor4(vertexList, &ia, &ib, &ic,&id, &val);
                        det_capo = val * iswap;
                }

                *test_abpo = det_abpo > 0;
                *test_bcpo = det_bcpo > 0;
                *test_capo = det_capo > 0;

                *convex = false;
                if((itest_abcp * det_abpo) > 0) return;
                if((itest_abcp * det_bcpo) > 0) return;
                if((itest_abcp * det_capo) > 0) return;
                *convex = true;

        }
        else if (ninf == 1)
        {       /*
                Define X as infinite point, and (i,j) the pair of finite points.
                If X = a, (i,j) = (b,c)
                If X = b, (i,j) = (c,a)
                If X = c, (i,j) = (a,b)
                If we define inf(a) = 1 if a infinite, 0 otherwise,
                then idx_X  = 2 - inf(a) + inf(c)
                */
                idx = 2 - infPoint[a] + infPoint[c];
                infp = list[idx] - npoints;
                i = list[order[idx][1]];
                j = list[order[idx][2]];

                //Get the coordinates
                for (m = 1; m <= 3; m++)
                {
                        i_p[m] = vertexList[i].Coordinates[m] - coordp[m];
                        j_p[m] = vertexList[j].Coordinates[m] - coordp[m];
                }

                if(infPoint[o] == 0)
                {
                        //First case:	o is finite
                        for (m = 1; m <= 3; m++)
                        {
                                o_p[m] = vertexList[o].Coordinates[m] - coordp[m];
                        }

                        icol1 = inf4_1[infp];
                        sign1 = sign4_1[infp];
                        /*
                        The three 4x4 determinants become:
                        -det(i,p,o) [X missing], det(j,p,o) [X missing],det(i,j,p,o)
                        And the 5x5 determinant becomes:
                        - det(i,j,p,o)
                        */
                        Mjo[1] = j_p[1] * o_p[3] - j_p[3] * o_p[1];
                        Mjo[2] = j_p[2] * o_p[3] - j_p[3] * o_p[2];
                        Mjo[3] = j_p[1] * o_p[2] - j_p[2] * o_p[1];
                        /*
                        The correspondence between a,b,c and i,j is not essential
                        We use corresponce for a infinite; in thetwo other cases
                        (b infinite or c infinite),computed determinants are the
                        same , but not in the same order
                        */
                        det_abpo = i_p[icol1] * o_p[3] - i_p[3] * o_p[icol1];
                        if(fabs(det_abpo) < eps)
                        {
                                int temp = 3;
                                ValSort3(i, p, o, &ia, &ib, &ic, &iswap);
                                sos->Minor3(vertexList, &ia, &ib, &ic,&icol1,&temp,&val);
                                det_abpo = -val * iswap;
                        }
                        det_abpo = det_abpo * sign1;
                        det_capo = - Mjo[icol1];
                        if(fabs(det_capo) < eps)
                        {
                                int temp = 3;
                                ValSort3(j, p, o, &ia, &ib, &ic, &iswap);
                                sos->Minor3(vertexList, &ia, &ib, &ic, &icol1, &temp, &val);
                                det_capo = val * iswap;
                        }
                        det_capo = det_capo * sign1;
                        det_bcpo =-i_p[1] * Mjo[2] + i_p[2] * Mjo[1] - i_p[3] * Mjo[3];
                        if(fabs(det_bcpo) < eps)
                        {
                                ValSort4(i, j, p, o, &ia, &ib, &ic, &id, &iswap);
                                sos->Minor4(vertexList, &ia, &ib, &ic, &id, &val);
                                det_bcpo = val * iswap;
                        }
                        det_abcpo = -det_bcpo;
                }
                else
                {
                        //Second case: o is infinite
                        info = o - npoints;
                        /*
                        The three 4x4 determinants become:
                        -det(i,p) [o,X missing]
                        det(j,p) [o,X missing]
                        det(i,j,p) [o missing]
                        And the 5x5 determinant becomes:
                        det(i,j,p) [o,X missing]
                        */
                        icol1 = inf4_2[infp][info];
                        sign1 = sign4_2[infp][info];

                        icol2 = inf4_1[info];
                        sign2 = sign4_1[info];

                        icol5 = inf5_2[infp][info];
                        sign5 = sign5_2[infp][info];

                        det_abpo = -i_p[icol1] * sign1;
                        if(fabs(det_abpo) < eps)
                        {
                                ValSort2(i, p, &ia, &ib, &iswap);
                                sos->Minor2(vertexList, &ia, &ib, icol1, &val);
                                det_abpo = -val * iswap * sign1;
                        }
                        det_capo = j_p[icol1] * sign1;
                        if(fabs(det_capo) < eps)
                        {
                                ValSort2(j, p, &ia, &ib, &iswap);
                                sos->Minor2(vertexList, &ia, &ib, icol1, &val);
                                det_capo = val * iswap * sign1;
                        }
                        det_bcpo = i_p[icol2] * j_p[3] - i_p[3] * j_p[icol2];
                        if(fabs(det_bcpo) < eps)
                        {
                                int temp = 3;
                                ValSort3(i, j, p, &ia, &ib, &ic, &iswap);
                                sos->Minor3(vertexList, &ia, &ib, &ic, &icol2, &temp,&val);
                                det_bcpo = val * iswap;
                        }
                        det_bcpo = det_bcpo * sign2;
                        det_abcpo = i_p[icol5] * j_p[3] - i_p[3] * j_p[icol5];
                        if(fabs(det_abcpo) < eps)
                        {
                                int temp = 3;
                                ValSort3(i, j, p, &ia, &ib, &ic, &iswap);
                                sos->Minor3(vertexList, &ia, &ib, &ic, &icol5, &temp, &val);
                                det_abcpo = val * iswap;
                        }
                        det_abcpo= det_abcpo * sign5;
                }

                //Test if (a,b,c,p,o) regular, in which case there is no need to flip
                if((det_abcpo * itest_abcp) < 0)
                {
                        *regular = true;
                        return;
                }
                *regular = false;

                //If not regular, we test for convexity

                testc[1] = det_abpo > 0;
                testc[2] = det_bcpo > 0;
                testc[3] = det_capo > 0;
                *test_abpo = testc[order_1[idx][1]];
                *test_bcpo = testc[order_1[idx][2]];
                *test_capo = testc[order_1[idx][3]];

                *convex = false;
                if ((itest_abcp * det_abpo) > 0) return;
                if ((itest_abcp * det_bcpo) > 0) return;
                if ((itest_abcp * det_capo) > 0) return;
                *convex = true;
        }
        else if (ninf == 2)
        {	/*
                define(k,l) as the two infinite points, and i be finite
                If i = a, (k,l) = (b,c)
                If i = b, (k,l) = (c,a)
                If i = c, (k,l) = (a,b)

                Again: idx = 2 + inf(a) - inf(c)
                */
                idx = 2 + infPoint[a] - infPoint[c];
                i = list[idx];
                k = list[order[idx][1]] - npoints;
                l = list[order[idx][2]] - npoints;

                //Get the coordinates

                for( m = 1; m <= 3; m++ )
                {
                        i_p[m] = vertexList[i].Coordinates[m] - coordp[m];
                }

                if (infPoint[o] == 0)
                {
                        //First case: o is finite
                        /*
                        The three 4x4 determinants become:
                        det(i,p,o) [k missing]
                        -det(i,p,o) [l missing]
                        S*det(p,o) [k,l missing, with S =1 if k<l, -1 otherwise]
                        The 5x5 determinants become:
                        S*det(i,p,o) [k,l missing, with S=1 if k<l, -1 otherwise]
                        */
                        for (m = 1; m <= 3; m++)
                        {
                                o_p[m] = vertexList[o].Coordinates[m] - coordp[m];
                        }
                        icol1 = inf4_1[k];
                        sign1 = sign4_1[k];
                        icol2 = inf4_1[l];
                        sign2 = sign4_1[l];
                        icol4 = inf4_2[l][k];
                        sign4 = sign4_2[l][k];
                        icol5 = inf5_2[l][k];
                        sign5 = sign5_2[l][k];

                        Mio[1] = i_p[1] * o_p[3] - i_p[3] * o_p[1];
                        Mio[2] = i_p[2] * o_p[3] - i_p[3] * o_p[2];
                        Mio[3] = i_p[1] * o_p[2] - i_p[2] * o_p[1];
                        /*
                        The correspondence between a,b,c and i,j,k is not essential
                        use the correspondence for a finite; in the two other cases
                        (b finite or c finite),have computed the same determinants,
                        but not in the same order
                        */
                        det_abpo = -Mio[icol1] * sign1;
                        if(fabs(det_abpo) < eps)
                        {
                                int temp = 3;
                                ValSort3(i, p, o, &ia, &ib, &ic, &iswap);
                                sos->Minor3(vertexList, &ia, &ib, &ic, &icol1, &temp, &val);
                                det_abpo = val * iswap * sign1;
                        }
                        det_capo = Mio[icol2] * sign2;
                        if(fabs(det_capo) < eps)
                        {
                                int temp = 3;
                                ValSort3(i, p, o, &ia, &ib, &ic, &iswap);
                                sos->Minor3(vertexList, &ia, &ib, &ic, &icol2, &temp, &val);
                                det_capo = -val * iswap * sign2;
                        }
                        det_bcpo = -o_p[icol4] * sign4;
                        if(fabs(det_bcpo) < eps)
                        {
                                ValSort2(p, o, &ia, &ib, &iswap);
                                sos->Minor2(vertexList, &ia, &ib, icol4, &val);
                                det_bcpo = val * sign4 * iswap;
                        }
                        det_abcpo = -Mio[icol5] * sign5;
                        if(fabs(det_abcpo) < eps)
                        {
                                int temp = 3;
                                ValSort3(i, p, o, &ia, &ib, &ic, &iswap);
                                sos->Minor3(vertexList, &ia, &ib, &ic, &icol5, &temp, &val);
                                det_abcpo = val * iswap * sign5;
                        }
                }
                else
                {
                        //Second case: o is infinite
                        info = o - npoints;
                        /*
                        The three 4x4 determinants become:
                        det(i,p) [o,k missing]
                        -det(i,p) [o,l missing]
                        Const [o,k,l missing]
                        The 5x5 determinants become:
                        Const*det(i,p) [o,k,l missing]
                        */
                        icol1 = inf4_2[k][info];
                        sign1 = sign4_2[k][info];
                        icol2 = inf4_2[l][info];
                        sign2 = sign4_2[l][info];

                        MissInfSign(info, k, l, &icol4, &iswap);

                        det_abpo = i_p[icol1] * sign1;
                        if(fabs(det_abpo) < eps)
                        {
                                ValSort2(i, p, &ia, &ib, &iswap2);
                                sos->Minor2(vertexList, &ia, &ib, icol1, &val);
                                det_abpo = val * iswap2 * sign1;
                        }
                        det_capo = -i_p[icol2] * sign2;
                        if(fabs(det_capo) < eps)
                        {
                                ValSort2(i, p, &ia, &ib, &iswap2);
                                sos->Minor2(vertexList, &ia, &ib,icol2,&val);
                                det_capo = -val * iswap2 * sign2;
                        }
                        det_bcpo = sign4_3[icol4] * iswap;
                        det_abcpo = sign5_3[icol4] * iswap * i_p[inf5_3[icol4]];
                        if(fabs(det_abcpo) < eps)
                        {
                                ValSort2(i, p, &ia, &ib, &iswap2);
                                sos->Minor2(vertexList,&ia,&ib,inf5_3[icol4],&val);
                                det_abcpo = val* iswap2 * iswap * sign5_3[icol4];
                        }
                }
                //if (a,b,c,p,o) regular, no need to flip

                if((det_abcpo * itest_abcp) < 0)
                {
                        *regular = true;
                        return;
                }
                *regular = false;

                //If not regular, we test for convexity

                testc[1] = det_abpo > 0;
                testc[2] = det_bcpo > 0;
                testc[3] = det_capo > 0;
                *test_abpo = testc[order_1[idx][1]];
                *test_bcpo = testc[order_1[idx][2]];
                *test_capo = testc[order_1[idx][3]];

                *convex = false;
                if((itest_abcp * det_abpo) > 0) return;
                if((itest_abcp * det_bcpo) > 0) return;
                if((itest_abcp * det_capo) > 0) return;
                *convex = true;
        }
        else if (ninf == 3)
        {
                assert(true);
                //this should not happen
        }
}

/*!
    \fn DeluanayComplex::DefineFacet(int itetra, int jtetra, int idx_o, int itouch[], int idx[], int jtouch[], int jdx[])
 */
void DeluanayComplex::DefineFacet(int itetra, int jtetra, int idx_o, int itouch[], int idx[], int jtouch[], int jdx[])
{
        int	i,k;
        int	Ia,Ib,Ie,If;

        int corresp[4];

        /*
        find the three vertices that define their common face.these vertices are stored in the array triangle and also find the vertices p and o
        To define the common face of the two tetrahedra itetra and jtetra,look at the neighbours of itetra : one of them is jtetra! this also provides p.
        The same procedure is repeated for jtetra,to get o
        */

        for( i = 1; i <= 3; i++ )
        {
                itouch[i] = DeluanayTet[itetra].Neighbours[i];
                idx[i] = DeluanayTet[itetra].Nindex[i];
        }

        Ia = DeluanayTet[itetra].Corners[1];

        for( i = 1; i <= 3; i++ )
        {
                k = other[idx_o][i];
                Ie = DeluanayTet[jtetra].Corners[k];
                if(Ia == Ie)
                {
                        corresp[1] = k;
                        break;
                }
        }
        Ib = DeluanayTet[itetra].Corners[2];
        Ie = other2[idx_o][corresp[1]][1];
        If = other2[idx_o][corresp[1]][2];
        if(Ib == DeluanayTet[jtetra].Corners[Ie])
        {
                corresp[2] = Ie;
                corresp[3] = If;
        }
        else
        {
                corresp[2] = If;
                corresp[3] = Ie;
        }

        for( i = 1; i <= 3; i++ )
        {
                k = corresp[i];
                jtouch[i] = DeluanayTet[jtetra].Neighbours[k];
                jdx[i] = DeluanayTet[jtetra].Nindex[k];
        }
}

/*!
    \fn DeluanayComplex::FindTetra(int itetra, int idx_c, int a, int b, int o, int *ifind, int *tetra_loc, int *idx_a, int *idx_b)
 */
void DeluanayComplex::FindTetra(int itetra, int idx_c, int a, int b, int o, int *ifind, int *tetra_loc, int *idx_a, int *idx_b)
{
        int i, ot, otx, otest;

        ot = DeluanayTet[itetra].Neighbours[idx_c];
        otx = DeluanayTet[itetra].Nindex[idx_c];

        otest = DeluanayTet[ot].Corners[otx];

        if (otest == o)
        {
                *ifind = 1;
                *tetra_loc = ot;

                //found the tetrahedron, define the position
                //of a and b in this tetrahedron

                for (i = 1; i <= 4; i++)
                {
                        if (DeluanayTet[*tetra_loc].Corners[i] == a)
                        {
                                *idx_a = i;
                        }
                        else if (DeluanayTet[*tetra_loc].Corners[i] == b)
                        {
                                *idx_b = i;
                        }
                }
        }
        else
        {
                *ifind = 0;
        }
}

/*!
    \fn DeluanayComplex::FlipJW2_3(int itetra, int jtetra, int p, int o, int itetra_touch[], int itetra_idx[], int jtetra_touch[], int jtetra_idx[], bool *test_abpo, bool *test_bcpo, bool *test_capo, int *ierr)
 */
void DeluanayComplex::FlipJW2_3(int itetra, int jtetra, int p, int o, int itetra_touch[], int itetra_idx[], int jtetra_touch[], int jtetra_idx[], bool *test_abpo, bool *test_bcpo, bool *test_capo, int *ierr)
{
        int i, j, k;
        int it, jt, idx, jdx;
        int newtetra;
        int tests[4];
        int position[4];

        ierr = 0;

        //if itetra or jtetra are inactive, cannot flip

        if((DeluanayTet[itetra].Status != 1) || (DeluanayTet[jtetra].Status != 1))
        {
                *ierr = 1;
                return;
        }
        /*
        The three new tetrahedra are going to be stored in : any free space in the tetrahedron list,
        and at the end of the list of known tetrahedra if needed
        */
        k = 0;
        for (i = nfree; i >= myMax(nfree - 2, 1); i--)
        {
                k++;
                position[k] = free[i];
                free.erase(free.begin()+i);
        }
        nfree = myMax(nfree - 3, 0);

        for (i = k+1; i <= 3; i++)
        {
                position[i] = DeluanayTet.size();    //add dummy tets at the end as placeholders and fill them later
                Tetrahedron tet;

                DeluanayTet.push_back(tet);
                ntetra++;
        }

        //mark itetra and jtetra to be 'killed' later and returned free list :P

        kill.push_back(itetra);
        kill.push_back(jtetra);
        nkill += 2;

        /*
        things of need:the vertices a,b,c are the first vertices of itetra, the vertices p and o, for each vertex in the triangle,
        define opposing faces in itetra and jtetra, & the tetrahedra that share that faces with itetra & jtetra, respectively.store
        this info in two arrays, itetra_touch and jtetra_touch.These information are provided by the calling program. for convenience
        always store p as the last vertex.

        define the three new tetrahedra: (bcop), (acop) and (abop) 	as well as their neighbours

        tetrahedron bcop : neighbours -> acop, abop, neighbour of (abcp) on face bcp, neighbour of (abco) on face bco
        tetrahedron acop : neighbours -> bcop, abop, neighbour of (abcp) on face acp, neighbour of (abco) on face aco
        tetrahedron abop : neighbours -> bcop, acop, neighbour of (abcp) on face abp, neighbour of (abco) on face abo
        */
        tests[1] = 1;
        if (*test_bcpo) { tests[1] = -1; }

        tests[2] = -1;
        if (*test_capo) { tests[2] = 1; }

        tests[3] = 1;
        if (*test_abpo) { tests[3] = -1; }

        //now fill the new tetrahedra already added

        for (i = 1; i <= 3; i++)
        {
                newtetra = position[i];
                k = 0;
                for (j = 1; j <= 3; j++)
                {
                        if (j == i) continue;
                        k++;
                        DeluanayTet[newtetra].Corners[k] = DeluanayTet[itetra].Corners[j];
                        DeluanayTet[newtetra].Neighbours[k] = position[j];
                        DeluanayTet[newtetra].Nindex[k] = idx_list[i][k];
                }

                DeluanayTet[newtetra].Corners[3] = o;
                it = itetra_touch[i];
                idx = itetra_idx[i];
                DeluanayTet[newtetra].Neighbours[3] = it;
                DeluanayTet[newtetra].Nindex[3] = idx;
                if ((idx != 0) && (it != 0))
                {
                        DeluanayTet[it].Neighbours[idx] = newtetra;
                        DeluanayTet[it].Nindex[idx] = 3;
                }

                DeluanayTet[newtetra].Corners[4] = p;
                jt = jtetra_touch[i];
                jdx = jtetra_idx[i];
                DeluanayTet[newtetra].Neighbours[4] = jt;
                DeluanayTet[newtetra].Nindex[4] = jdx;
                if ((jdx != 0) && (jt != 0))
                {
                        DeluanayTet[jt].Neighbours[jdx] = newtetra;
                        DeluanayTet[jt].Nindex[jdx] = 4;
                }

                DeluanayTet[newtetra].Status = 1;
                DeluanayTet[newtetra].Orient = tests[i];
        }
        /*
        add all three faces of jtetra containing o in the link_facet queue.
        Each link_facet (a triangle) is implicitly defined as the intersection of two tetrahedra

        link_facet:	bco     -> bcop & neighbour of (abco) on bco
        link_facet:	aco	-> acop & neighbour of (abco) on aco
        link_facet:	abo	-> abop & neighbour of (abco) on abo
        */
        for (i = 1; i <= 3; i++)
        {
                newtetra = position[i];
                LinkFacet lf;
                LinkIndex li;
                lf.f1 = newtetra;
                lf.f2 = DeluanayTet[newtetra].Neighbours[4];
                li.i1 = 4;
                li.i2 = DeluanayTet[newtetra].Nindex[4];
                linkFacet.push_back(lf);
                linkIndex.push_back(li);
        }
        //kill the old tetrahedra
        DeluanayTet[itetra].Kill();
        DeluanayTet[jtetra].Kill();
}

/*!
    \fn DeluanayComplex::FlipJW3_2(int itetra, int jtetra, int ktetra, int vertices[], int itetra_touch[], int itetra_idx[], int jtetra_touch[], int jtetra_idx[], int ktetra_touch[], int ktetra_idx[], bool *test_bcpo, bool *test_acpo,  int *ierr)
 */
void DeluanayComplex::FlipJW3_2(int itetra, int jtetra, int ktetra, int vertices[], int itetra_touch[], int itetra_idx[], int jtetra_touch[], int jtetra_idx[], int ktetra_touch[], int ktetra_idx[], bool *test_bcpo, bool *test_acpo,  int *ierr)
{
        int i, j, k, p, o, c;
        int it, jt, kt, idx, jdx, kdx, newtetra;
        int edge[3];
        int tests[3];
        int position[3];

        tests[1] = 1;
        if (*test_bcpo) { tests[1] = -1; }

        tests[2] = 1;
        if (*test_acpo) { tests[2] = -1; }

        *ierr = 0;

        //If itetra, jtetra or ktetra are inactive, cannot flip

        if ((DeluanayTet[itetra].Status != 1) || (DeluanayTet[jtetra].Status != 1) || (DeluanayTet[ktetra].Status != 1))
        {
                *ierr = 1;
                return;
        }

        edge[1] = vertices[1];
        edge[2] = vertices[2];
        c = vertices[3];
        p = vertices[4];
        o = vertices[5];

        k = 0;
        for (i = nfree; i >= myMax(nfree - 1, 1); i--)
        {
                k++;
                position[k] = free[i];
                free.erase(free.begin()+i);
        }
        nfree = myMax(nfree - 2, 0);

        for (i = k+1; i <= 2; i++)
        {
                position[i] = DeluanayTet.size(); //as usual add two dummy tets and fill them later
                Tetrahedron tet;

                DeluanayTet.push_back(tet);
                ntetra++;
        }

        //add itetra, jtetra and ktetra to the "kill" list so that it can be returned to free list:P
        kill.push_back(itetra);
        kill.push_back(jtetra);
        kill.push_back(ktetra);
        nkill += 3;

        /*
        things needed:two vertices that define common edge (ab) these are stored in the array edge,the vertices c, p and o that form the new triangle
        for each vertex in the edge (ab), define the opposing faces in the three tetrahedra itetra, jtetra & ktetra,the tetrahedron that share these faces
        with itetra, jtetra & ktetra, respectively.
        This information is stored in three arrays, itetra_touch, jtetra_touch and ktetra_touch.These information are given by the calling program

        for convenience set p to be the last vertex of the new tetrahedra

        define the two new tetrahedra: (bcop) and (acop) as well as their neighbours

        tetrahedron bcop : neighbours -> acop, neighbour of (abop) on face bpo, neighbour of (abcp) on face bcp & neighbour of (abco) on face (bco)
        tetrahedron acop : neighbours -> bcop, neighbour of (abop) on face apo, neighbour of (abcp) on face acp & neighbour of (abco) on face (aco)
        */

        for (i = 1; i <= 2; i++)
        {
                newtetra = position[i];
                k = 0;
                for (j = 1; j <= 2; j++)
                {
                        if (j == i) continue;
                        k++;
                        DeluanayTet[newtetra].Corners[k] = edge[j];
                        DeluanayTet[newtetra].Neighbours[k] = position[j];
                        DeluanayTet[newtetra].Nindex[k] = 1;
                }

                DeluanayTet[newtetra].Corners[2] = c;
                kt = ktetra_touch[i];
                kdx = ktetra_idx[i];
                DeluanayTet[newtetra].Neighbours[2] = kt;
                DeluanayTet[newtetra].Nindex[2] = kdx;
                if((kdx != 0) && (kt != 0))
                {
                        DeluanayTet[kt].Neighbours[kdx] = newtetra;
                        DeluanayTet[kt].Nindex[kdx] = 2;
                }

                DeluanayTet[newtetra].Corners[3] = o;
                it = itetra_touch[i];
                idx = itetra_idx[i];
                DeluanayTet[newtetra].Neighbours[3] = it;
                DeluanayTet[newtetra].Nindex[3] = idx;
                if ((idx != 0) && (it != 0))
                {
                        DeluanayTet[it].Neighbours[idx] = newtetra;
                        DeluanayTet[it].Nindex[idx] = 3;
                }

                DeluanayTet[newtetra].Corners[4] = p;
                jt = jtetra_touch[i];
                jdx = jtetra_idx[i];
                DeluanayTet[newtetra].Neighbours[4] = jt;
                DeluanayTet[newtetra].Nindex[4] = jdx;
                if ((jdx != 0) && (jt != 0))
                {
                        DeluanayTet[jt].Neighbours[jdx] = newtetra;
                        DeluanayTet[jt].Nindex[jdx] = 4;
                }

                DeluanayTet[newtetra].Status = 1;
                DeluanayTet[newtetra].Orient = tests[i];
        }

        /*
        Now add the two faces of ktetra containing (co) in the link_facet queue.
        Each link_facet is defined as the intersection of two tetrahedra

        link_facet:	bco	-> bcop and neighbour of (abco) on bco
        link_facet:	aco	-> acop and neighbour of (abco) on aco
        */

        for (i = 1; i <= 2; i++)
        {
                newtetra = position[i];
                LinkFacet lf;
                LinkIndex li;
                lf.f1 = newtetra;
                lf.f2 = DeluanayTet[newtetra].Neighbours[4];
                li.i1 = 4;
                li.i2 = DeluanayTet[newtetra].Nindex[4];
                linkFacet.push_back(lf);
                linkIndex.push_back(li);
        }
        DeluanayTet[itetra].Kill();
        DeluanayTet[jtetra].Kill();
        DeluanayTet[ktetra].Kill();
}

/*!
    \fn DeluanayComplex::FlipJW4_1(std::vector<Vertex>& vertexList,int itetra, int jtetra, int ktetra, int ltetra, int vertices[], int ishare, int idx, int jshare, int jdx, int kshare, int kdx, int lshare, int ldx, bool *test_acpo, int *ierr)
 */
void DeluanayComplex::FlipJW4_1(std::vector<Vertex>& vertexList,int itetra, int jtetra, int ktetra, int ltetra, int vertices[], int ishare, int idx, int jshare, int jdx, int kshare, int kdx, int lshare, int ldx, bool *test_acpo, int *ierr)
{
        int p, o, a, b, c;
        int test1, newtetra;

        *ierr = 0;

        test1 = 1;
        if (*test_acpo) { test1 = -1; }

        //if any one of 4 tets are inactive cant flip

        if ((DeluanayTet[itetra].Status != 1) || (DeluanayTet[jtetra].Status != 1) || (DeluanayTet[ktetra].Status != 1) || (DeluanayTet[ltetra].Status != 1))
        {
                *ierr = 1;
                return;
        }

        if (nfree != 0)
        {
                newtetra = free[nfree];
                free.erase(free.begin()+nfree);
                nfree--;
        }
        else
        {
                newtetra = DeluanayTet.size();
                Tetrahedron tet;

                DeluanayTet.push_back(tet);               //as usual add dummy tet first and update later
                ntetra++;
        }
        //set all tets to be killed later added to free list :P
        kill.push_back(itetra);
        kill.push_back(jtetra);
        kill.push_back(ktetra);
        kill.push_back(ltetra);
        nkill += 4;
        /*
        things needed : the vertex b shared by all 4 tetrahedra,the vertices a, c, p and o,for each tetrahedron, neighbour attached to the face opposite to b;
        this information is stored in ishare,jshare,kshare,lshare.These information are provided by the calling program
        */
        a = vertices[1];
        b = vertices[2];
        c = vertices[3];
        p = vertices[4];
        o = vertices[5];
        /*
        for convenience p is set to be last vertex of tetrahedron

        new tetrahedron: (acop) -> neighbor of (bcop) on face cpo, neighbor of (abop) on face apo,
        neighbor of (abcp) on face acp, neighbor of (abco) on face aco
        */
        //mark vertex b as redundant
        vertexList[b].Redinfo = 1;

        DeluanayTet[newtetra].Corners[1] = a;
        DeluanayTet[newtetra].Neighbours[1] = lshare;
        DeluanayTet[newtetra].Nindex[1] = ldx;
        if ((lshare != 0) && (ldx != 0))
        {
                DeluanayTet[lshare].Neighbours[ldx] = newtetra;
                DeluanayTet[lshare].Nindex[ldx] = 1;
        }

        DeluanayTet[newtetra].Corners[2] = c;
        DeluanayTet[newtetra].Neighbours[2] = kshare;
        DeluanayTet[newtetra].Nindex[2] = kdx;
        if ((kshare != 0) && (kdx != 0))
        {
                DeluanayTet[kshare].Neighbours[kdx] = newtetra;
                DeluanayTet[kshare].Nindex[kdx] = 2;
        }

        DeluanayTet[newtetra].Corners[3] = o;
        DeluanayTet[newtetra].Neighbours[3] = ishare;
        DeluanayTet[newtetra].Nindex[3] = idx;
        if ((ishare != 0) && (idx != 0))
        {
                DeluanayTet[ishare].Neighbours[idx] = newtetra;
                DeluanayTet[ishare].Nindex[idx] = 3;
        }

        DeluanayTet[newtetra].Corners[4] = p;
        DeluanayTet[newtetra].Neighbours[4] = jshare;
        DeluanayTet[newtetra].Nindex[4] = jdx;
        if ((jshare != 0) && (jdx != 0))
        {
                DeluanayTet[jshare].Neighbours[jdx] = newtetra;
                DeluanayTet[jshare].Nindex[jdx] = 4;
        }

        DeluanayTet[newtetra].Status = 1;
        DeluanayTet[newtetra].Orient = test1;
        /*
        Now add one link facet
        link_facet:	aco	-> acop and neighbour of (abco) on aco
        */
        LinkFacet lf;
        LinkIndex li;
        lf.f1 = newtetra;
        lf.f2 = jshare;
        li.i1 = 4;
        li.i2 = jdx;
        linkFacet.push_back(lf);
        linkIndex.push_back(li);

        DeluanayTet[itetra].Kill();
        DeluanayTet[jtetra].Kill();
        DeluanayTet[ktetra].Kill();
        DeluanayTet[ltetra].Kill();
}

/*!
    \fn DeluanayComplex::FlipJW(std::vector<Vertex> & vertexList)
 */
void DeluanayComplex::FlipJW(std::vector<Vertex> & vertexList)
{
        uint j;
        int i,p,o,a,b,c;
        int ierr = 0,ifind = 0;
        int itetra,jtetra;
        int tetra_ab = 0, tetra_ac = 0, tetra_bc = 0;
        int iorder;
        int ireflex,iflip;
        int idx_p,idx_o,itest_abcp;
        int idx_a = 0, idx_b = 0, idx_c = 0;
        int ishare,jshare,kshare,lshare;
        int idxi,idxj,idxk,idxl;
        int ia,ib,ic,ii,ij;
        int itetra_touch[4];
        int jtetra_touch[4];
        int itetra_idx[4];
        int jtetra_idx[4];
        int ktetra_touch[4][3];
        int ktetra_idx[4][3];
        int itouch[3];
        int jtouch[3];
        int ktouch[3];
        int i_idx[3];
        int j_idx[3];
        int k_idx[3];
        int tetra_flip[4];
        int list_flip[4];
        int vert_flip[6];

        bool test_or[4][3];
        bool test,regular = false,convex = false;
        bool test_abpo = false,test_abpc,test_capo = false,test_acpb;
        bool test_bcpo = false,test_bcpa,test_acpo;

        for (j = 1; j < linkFacet.size(); j++)
        {
                //get the tetras of which the given linkfacet is common

                itetra = linkFacet[j].f1;
                jtetra = linkFacet[j].f2;
                idx_p = linkIndex[j].i1;
                idx_o = linkIndex[j].i2;

                //if linkfacet is on convex hull discard them
                if ((itetra == 0) || (jtetra == 0)) continue;

                //if these tetrahedra are already discarded then discard the linkfacet
                if (DeluanayTet[itetra].Status == 0)
                {
                        if (DeluanayTet[jtetra].Status == 0)
                                continue;
                        else
                        {
                                itetra = DeluanayTet[jtetra].Neighbours[idx_o];
                                idx_p = DeluanayTet[jtetra].Nindex[idx_o];
                        }
                }
                if(DeluanayTet[jtetra].Status == 0)
                {
                        jtetra = DeluanayTet[itetra].Neighbours[idx_p];
                        idx_o = DeluanayTet[itetra].Nindex[idx_p];
                }
                /*Define the vertices of the two tetrahedra:
                itetra:		a,b,c,p
                jtetra:		a,b,c,o
                */
                a = DeluanayTet[itetra].Corners[1];
                b = DeluanayTet[itetra].Corners[2];
                c = DeluanayTet[itetra].Corners[3];
                p = DeluanayTet[itetra].Corners[4];

                o = DeluanayTet[jtetra].Corners[idx_o];

                itest_abcp = DeluanayTet[itetra].Orient;

                //Check for local regularity (and convexity, at very little extra cost)

                RegularConvex(vertexList,a, b, c, p, o, itest_abcp, &regular, &convex, &test_abpo, &test_bcpo, &test_capo);

                //if link facet is locally regular there is nothing to do but discard
                if (regular) continue;

                //define neighbours of the facet in itetra and jtetra
                DefineFacet(itetra, jtetra, idx_o, itetra_touch, itetra_idx, jtetra_touch, jtetra_idx);

                test_abpc = (itest_abcp != 1);

                /*
                After discarding the trivial case, we now test if the tetrahedra
                can be flipped.

                At this stage, we know that the link facet is not locally
                regular. but we still don't know if it is "flippable"

                if {itetra} U {jtetra} is convex (using the convexity test performed
                with the regularity test) perform a 2-3 flip.
                */
                if (convex)
                {
                        FlipJW2_3(itetra, jtetra, p, o, itetra_touch, itetra_idx, jtetra_touch, jtetra_idx, &test_abpo, &test_bcpo, &test_capo, &ierr);
                        continue;
                }
                /*
                The union of the two tetrahedra is not convex...
                check the edges of the triangle in the link facet, and
                check if they are "reflexes" (see definition in Edelsbrunner and
                Shah, Algorithmica (1996), 15:223-241)
                */
                ireflex = 0;
                iflip = 0;

                /*
                First check edge (ab):
                - (ab) is reflex iff o and c lies on opposite sides of
                the hyperplane defined by (abp). We therefore test the
                orientation of (abpo) and (abpc): if they differ (ab)
                is reflex
                - if (ab) is reflex, we test if it is of degree 3.
                (ab) is of degree 3 if it is shared by 3 tetrahedra,
                namely (abcp), (abco) and (abpo). The first two are itetra
                and jtetra, so we only need to check if (abpo) exists.
                since (abpo) contains p, (abp) should then be a link facet
                of p, so we test all tetrahedra that define link facets
                */

                if(test_abpo != test_abpc)
                {
                        ireflex = ireflex + 1;

                        FindTetra(itetra, 3, a, b, o, &ifind, &tetra_ab, &idx_a, &idx_b);

                        if(ifind == 1)
                        {
                                iflip = iflip + 1;
                                tetra_flip[iflip] = tetra_ab;
                                list_flip[iflip] = 1;
                                ktetra_touch[iflip][1] = DeluanayTet[tetra_ab].Neighbours[idx_a];
                                ktetra_touch[iflip][2] = DeluanayTet[tetra_ab].Neighbours[idx_b];
                                ktetra_idx[iflip][1] = DeluanayTet[tetra_ab].Nindex[idx_a];
                                ktetra_idx[iflip][2]= DeluanayTet[tetra_ab].Nindex[idx_b];
                                test_or[iflip][1] = test_bcpo;
                                test_or[iflip][2] = (!test_capo);
                        }
                }

                /*
                Now check edge (ac):
                - (ac) is reflex iff o and b lies on opposite sides of
                the hyperplane defined by (acp). We therefore test the
                orientation of (acpo) and (acpb): if they differ (ac)
                is reflex
                - if (ac) is reflex, we test if it is of degree 3.
                (ac) is of degree 3 if it is shared by 3 tetrahedra,
                namely (abcp), (abco) and (acpo). The first two are itetra
                and jtetra, so we only need to check if (acpo) exists.
                since (acpo) contains p, (acp) should then be a link facet
                of p, so we test all tetrahedra that define link facets
                */
                test_acpo = !test_capo;
                test_acpb = !test_abpc;

                if(test_acpo != test_acpb)
                {
                        ireflex = ireflex + 1;

                        FindTetra(itetra, 2, a, c, o, &ifind, &tetra_ac, &idx_a, &idx_c);

                        if(ifind == 1)
                        {
                                iflip = iflip + 1;
                                tetra_flip[iflip] = tetra_ac;
                                list_flip[iflip] = 2;
                                ktetra_touch[iflip][1] = DeluanayTet[tetra_ac].Neighbours[idx_a];
                                ktetra_touch[iflip][2] = DeluanayTet[tetra_ac].Neighbours[idx_c];
                                ktetra_idx[iflip][1] = DeluanayTet[tetra_ac].Nindex[idx_a];
                                ktetra_idx[iflip][2] = DeluanayTet[tetra_ac].Nindex[idx_c];
                                test_or[iflip][1] = !test_bcpo;
                                test_or[iflip][2] = test_abpo;
                        }
                }
                /*
                Now check edge (bc):
                - (bc) is reflex iff o and a lies on opposite sides of
                the hyperplane defined by (bcp). We therefore test the
                orientation of (bcpo) and (bcpa): if they differ (bc)
                is reflex
                - if (bc) is reflex, we test if it is of degree 3.
                (bc) is of degree 3 if it is shared by 3 tetrahedra,
                namely (abcp), (abco) and (bcpo). The first two are itetra
                and jtetra, so we only need to check if (bcpo) exists.
                since (bcpo) contains p, (bcp) should then be a link facet
                of p, so we test all tetrahedra that define link facets
                */

                test_bcpa = test_abpc;

                if(test_bcpo != test_bcpa)
                {
                        ireflex = ireflex + 1;

                        FindTetra(itetra, 1, b, c, o, &ifind, &tetra_bc, &idx_b, &idx_c);

                        if(ifind == 1)
                        {
                                iflip = iflip + 1;
                                tetra_flip[iflip] = tetra_bc;
                                list_flip[iflip] = 3;
                                ktetra_touch[iflip][1] = DeluanayTet[tetra_bc].Neighbours[idx_b];
                                ktetra_touch[iflip][2] = DeluanayTet[tetra_bc].Neighbours[idx_c];
                                ktetra_idx[iflip][1]  = DeluanayTet[tetra_bc].Nindex[idx_b];
                                ktetra_idx[iflip][2]  = DeluanayTet[tetra_bc].Nindex[idx_c];
                                test_or[iflip][1] = test_capo;
                                test_or[iflip][2] = !test_abpo;
                        }
                }
                if (ireflex != iflip) continue;

                if (iflip == 1)
                {
                        //Only one edge is "flippable": we do a 3-2 flip
                        iorder = list_flip[iflip];
                        ia = table32[iorder][1];
                        ib = table32[iorder][2];
                        ic = table32[iorder][3];
                        vert_flip[ia] = a;
                        vert_flip[ib] = b;
                        vert_flip[ic] = c;
                        vert_flip[4] = p;
                        vert_flip[5] = o;
                        ia = table32_2[iorder][1];
                        ib = table32_2[iorder][2];
                        itouch[1] = itetra_touch[ia];
                        itouch[2] = itetra_touch[ib];
                        i_idx[1]  = itetra_idx[ia];
                        i_idx[2]  = itetra_idx[ib];
                        jtouch[1] = jtetra_touch[ia];
                        jtouch[2] = jtetra_touch[ib];
                        j_idx[1]  = jtetra_idx[ia];
                        j_idx[2]  = jtetra_idx[ib];
                        ktouch[1] = ktetra_touch[iflip][1];
                        ktouch[2] = ktetra_touch[iflip][2];
                        k_idx[1]  = ktetra_idx[iflip][1];
                        k_idx[2]  = ktetra_idx[iflip][2];
                        FlipJW3_2(itetra, jtetra, tetra_flip[1], vert_flip, itouch, i_idx, jtouch, j_idx, ktouch, k_idx, &test_or[iflip][1], &test_or[iflip][2], &ierr);
                }
                else if (iflip == 2)
                {
                /*
                        In this case, one point is redundant: the point common to
                        the two edges that can be flipped.perform a 4-1 flip
                */
                        iorder = list_flip[1] + list_flip[2] - 2;
                        vert_flip[table41[iorder][1]] = a;
                        vert_flip[table41[iorder][2]] = b;
                        vert_flip[table41[iorder][3]] = c;
                        vert_flip[4] = p;
                        vert_flip[5] = o;
                        ii = table41_2[iorder][1];
                        ij = table41_2[iorder][2];
                        ishare = itetra_touch[iorder];
                        idxi   = itetra_idx[iorder];
                        jshare = jtetra_touch[iorder];
                        idxj   = jtetra_idx[iorder];
                        kshare = ktetra_touch[1][ii];
                        idxk   = ktetra_idx[1][ii];
                        lshare = ktetra_touch[2][ij];
                        idxl   = ktetra_idx[2][ij];
                        if(iorder == 1)
                        {
                                test = test_bcpo;
                        }
                        else if(iorder == 2)
                        {
                                test = !test_capo;
                        }
                        else
                        {
                                test = test_abpo;
                        }
                        FlipJW4_1(vertexList,itetra, jtetra, tetra_flip[1], tetra_flip[2], vert_flip, ishare, idxi, jshare, idxj, kshare, idxk, lshare, idxl, &test, &ierr);
                }
                else
                {
                        assert(true);
                }
        }
        //very important: empty the linkfacet list after everything is done
        linkFacet.clear();
        linkIndex.clear();
        AddDummyLF();

        //Add killed tetra_ab to free zone
        for (i = 1; i <= nkill; i++)
        {
                free.push_back(kill[i]);
        }
        nfree += nkill;
}

/*!
    \fn DeluanayComplex::FlipJW1_4(int ipoint, int itetra)
 */
void DeluanayComplex::FlipJW1_4(int ipoint, int itetra)
{
        int vertex[5];
        int neighbours[5];
        int nindex[5];
        int position[5];
        int fact,newtetindex,k,jtetra,idx;

        for (int i = 1; i <= 4; i++)
        {
                vertex[i] = DeluanayTet[itetra].Corners[i];
                neighbours[i] = DeluanayTet[itetra].Neighbours[i];
                nindex[i] = DeluanayTet[itetra].Nindex[i];
        }

        //
        k = 0;
        for (int i = nfree; i >= myMax((nfree - 3), 1); i--)
        {
                k++;
                position[k] = free[i];
                free.erase(free.begin()+i);
        }

        nfree = myMax((nfree - 4), 0);

        for (int i = k + 1; i <= 4; i++)
        {
                //new tetra added will be added towards end only if there is no slots free
                position[i] = DeluanayTet.size();
                Tetrahedron tet;

                DeluanayTet.push_back(tet);
                ntetra++;
        }

        fact = DeluanayTet[itetra].Orient;

        //The four new tetrahedra are going to be stored
        //in : any free space in the tetrahedron list,
        //and at the end of the list of known tetrahedra
        //no free list is maintained so no need of this processing

        //set itetra to be killed :P
        kill.clear();
        kill.push_back(0);
        kill.push_back(itetra);
        nkill = 1;

        /*
        The tetrahedron is defined as (ijkl);
        four new tetrahedra are created: jklp, iklp, ijlp, and ijkp,
        where p is the new point to be included

        For each new tetrahedron, define all four neighbours:
        For each neighbour, store the index of the vertex opposite to
        the common face in Nindex

        tetrahedron jklp : neighbours are iklp, ijlp, ijkp and neighbour of (ijkl) on face jkl
        tetrahedron iklp : neighbours are jklp, ijlp, ijkp and neighbour of (ijkl) on face ikl
        tetrahedron ijlp : neighbours are jklp, iklp, ijkp and neighbour of (ijkl) on face ijl
        tetrahedron ijkp : neighbours are jklp, iklp, ijlp and neighbour of (ijkl) on face ijk
        */

        for (int i = 1; i <= 4; i++)
        {
                newtetindex = position[i];

                k = 0;

                for (int j = 1; j <= 4; j++)
                {
                        if (j == i) continue;

                        k++;
                        DeluanayTet[newtetindex].Corners[k] = vertex[j];
                        DeluanayTet[newtetindex].Neighbours[k] = position[j];
                        DeluanayTet[newtetindex].Nindex[k] = idxList[i][k];
                }

                jtetra = neighbours[i];
                idx = nindex[i];
                DeluanayTet[newtetindex].Corners[4] = ipoint;
                DeluanayTet[newtetindex].Neighbours[4] = jtetra;
                DeluanayTet[newtetindex].Nindex[4] = idx;

                //update the neighbors of the neighbour of itetra!
                if(jtetra != 0 && idx != 0)
                {
                        DeluanayTet[jtetra].Neighbours[idx] = newtetindex;
                        DeluanayTet[jtetra].Nindex[idx] = 4;
                }
                DeluanayTet[newtetindex].Status = 1;
                fact = -fact;
                DeluanayTet[newtetindex].Orient = fact;
        }
        /*
        Now add all fours faces of itetra in the link_facet queue.
        Each link_facet (a triangle) is implicitly defined as the
        intersection of two tetrahedra

        link_facet:	jkl	tetrahedra:	jklp and neighbour of (ijkl) on jkl
        link_facet:	ikl	tetrahedra:	iklp and neighbour of (ijkl) on ikl
        link_facet:	ijl	tetrahedra:	ijlp and neighbour of (ijkl) on ijl
        link_facet:	ijk	tetrahedra:	ijkp and neighbour of (ijkl) on ijk
        */
        for (int i = 1; i <= 4; i++)
        {
                newtetindex = position[i];
                LinkFacet lf;
                LinkIndex li;
                lf.f1 = newtetindex;
                lf.f2 = DeluanayTet[newtetindex].Neighbours[4];
                li.i1 = 4;
                li.i2 = DeluanayTet[newtetindex].Nindex[4];
                linkFacet.push_back(lf);
                linkIndex.push_back(li);
        }

        //also kill the tet in question
        DeluanayTet[itetra].Kill();
}

/*!
    \fn DeluanayComplex::InitFreeKill()
 */
void DeluanayComplex::InitFreeKill()
{
        free.push_back(0);
        kill.push_back(0);
}

/*!
    \fn DeluanayComplex::Jump(std::vector<Vertex> & vertexList, int *iseed, int ival, int *itetra)
 */
void DeluanayComplex::Jump(std::vector<Vertex> & vertexList, int *iseed, int ival, int *itetra)
{
        int MaxEntry, n, j, nkeep = 0;
        int tempList[30];
        double xval[4];
        double xa[4];
        double dist, distmin;

        MaxEntry = 20;

        n = myMin(MaxEntry, ntetra);

        while (true)
        {
                for (int i = 0; i < n; i++)
                {
                        j = (int)(ntetra * ran2(iseed)) + 1;
                        j = myMin(j, ntetra);
                        if (DeluanayTet[j].Status == 0) continue;
                        nkeep++;
                        tempList[nkeep] = j;
                }
                if (nkeep != 0) break;
        }

        for (int i = 1; i <= 3; i++)
        {
                xval[i] = vertexList[ival].Coordinates[i];
                if (DeluanayTet[tempList[1]].Corners[1] >= vertexList.size())
                {
                        xa[i] = vertexList[0].Coordinates[i];
                }
                else
                {
                        xa[i] = vertexList[DeluanayTet[tempList[1]].Corners[1]].Coordinates[i];
                }
        }

        distmin = 0.0;
        for (int i = 1; i <= 3; i++)
        {
                distmin = distmin + pow((xval[i] - xa[i]), 2.0);
        }
        *itetra = tempList[1];

        for (int i = 2; i < nkeep; i++)
        {
                for (j = 1; j <= 3; j++)
                {
                        if (DeluanayTet[tempList[i]].Corners[1] >= vertexList.size())
                        {
                                xa[j] = vertexList[0].Coordinates[j];
                        }
                        else
                        {
                                xa[j] = vertexList[DeluanayTet[tempList[i]].Corners[1]].Coordinates[j];
                        }
                }
                dist = 0.0;
                for (j = 1; j <= 3; j++)
                {
                        dist = dist + pow((xval[j] - xa[j]), 2.0);
                }
                if (dist < distmin)
                {
                        distmin = dist;
                        *itetra = tempList[i];
                }
        }
}

/*!
    \fn DeluanayComplex::InsideTetraJW(std::vector<Vertex> & vertexList,int p, int a, int b, int c, int d, int iorient, int *is_in, int *redundant, int *ifail)
 */
void DeluanayComplex::InsideTetraJW(std::vector<Vertex> & vertexList,int p, int a, int b, int c, int d, int iorient, int *is_in, int *redundant, int *ifail)
{
        int	i,j,k,l;
        int	ia = 0,ib = 0,ic = 0,id = 0,ie = 0,idx;
        int	ic1,ic5,ic1_k,ic1_l,sign,sign5,sign_k,sign_l;
        int	npoints = vertexList.size() - 1;//is  - 1 needed;
        int	nswap = 0,iswap = 0,ninf;
        int	val = 0;

        int list[5];

        double	Sij_1,Sij_2,Sij_3,Skl_1,Skl_2,Skl_3;
        double	det_pijk,det_pjil,det_pkjl,det_pikl,det_pijkl;

        double detij[4];
        double i_p[5];
        double j_p[5];
        double k_p[5];
        double l_p[5];
        bool test_pijk,test_pjil,test_pkjl,test_pikl;
        /*
        If (i,j,k,l) is the tetrahedron in positive orientation, we need to test:
        (p,i,j,k)
        (p,j,i,l)
        (p,k,j,l)
        (p,i,k,l)
        If all four are positive, than p is inside the tetrahedron.
        All four tests relies on the sign of the corresponding 4x4
        determinant. Interestingly, these four determinants share
        some common lines, which can be used to speed up the computation.

        Let us consider or example:

        det(p,i,j,k) =  | p(1) p(2) p(3) 1|
                        | i(1) i(2) i(3) 1|
                        | j(1) j(2) j(3) 1|
                        | k(1) k(2) k(3) 1|

        p appears in each determinant. The corresponding line can therefore
        be substraced from all 3 other lines . Using the example above,
        we find:

        det(i,j,k,l) = -|ip(1) ip(2) ip(3)|
                        |jp(1) jp(2) jp(3)|
                        |kp(1) kp(2) kp(3)|

        where :xp(m) = x(m) - p(m) for x = i,j,k and m = 1,2,3

        Now we notice that the first two lines of det(p,i,j,k) and det(p,i,j,l) are the same.

        Let us define: Sij_3 = |ip(1) ip(2)| Sij_2 = |ip(1) ip(3)| and Sij_1 = |ip(2) ip(3)|
                               |jp(1) jp(2)|         |jp(1) jp(3)|             |jp(2) jp(3)|

        We find:
        det(p,i,j,k) = - kp(1)*Sij_1 + kp(2)*Sij_2 - kp(3)*Sij_3
        and:
        det(p,j,i,l) =   lp(1)*Sij_1 - lp(2)*Sij_2 + lp(3)*Sij_3

        Similarly, if we define:

        Skl_3 = |kp(1) kp(2)|	Skl_2 = |kp(1) kp(3)|	Skl_1 = |kp(2) kp(3)|
                |lp(1) lp(2)|		|lp(1) lp(3)|		|lp(2) lp(3)|

        We find:
        det(p,k,j,l) = jp(1)*Skl_1 - jp(2)*Skl_2 + jp(3)*Skl_3
        and:
        det(p,i,k,l) = -ip(1)*Skl_1 + ip(2)*Skl_2 - ip(3)*Skl_3

        Furthermore:
        det(p,i,j,k,l) = -ip(4)*det(p,k,j,l)-jp(4)*det(p,i,k,l)
                        -kp(4)*det(p,j,i,l)-lp(4)*det(p,i,j,k)

        The equations above hold for the general case; special care is
        required to take in account infinite points (see below)
        */

        list[1] = a;
        list[2] = b;
        list[3] = c;
        list[4] = d;

        ninf = infPoint[a] + infPoint[b] + infPoint[c] + infPoint[d];

        if (ninf == 0) //no infinite points
        {
                for( i = 1; i < 4; i++ )
                {
                        i_p[i] = vertexList[a].Coordinates[i] - vertexList[p].Coordinates[i];
                        j_p[i] = vertexList[b].Coordinates[i] - vertexList[p].Coordinates[i];
                        k_p[i] = vertexList[c].Coordinates[i] - vertexList[p].Coordinates[i];
                        l_p[i] = vertexList[d].Coordinates[i] - vertexList[p].Coordinates[i];
                }

                //Now compute 2x2 determinants Sij and Skl
                Sij_1 = i_p[2] * j_p[3] - i_p[3] * j_p[2];
                Sij_2 = i_p[1] * j_p[3] - i_p[3] * j_p[1];
                Sij_3 = i_p[1] * j_p[2] - i_p[2] * j_p[1];

                Skl_1 = k_p[2] * l_p[3] - k_p[3] * l_p[2];
                Skl_2 = k_p[1] * l_p[3] - k_p[3] * l_p[1];
                Skl_3 = k_p[1] * l_p[2] - k_p[2] * l_p[1];

                //Now perform tests
                //Start with is_in = .false. :
                *is_in = 0;//false

                //We check all other four determinants
                det_pijk = -k_p[1] * Sij_1 + k_p[2] * Sij_2 - k_p[3] * Sij_3;
                det_pijk = det_pijk * iorient;
                test_pijk = (fabs(det_pijk) > eps);
                if(test_pijk && det_pijk > 0)
                {
                        *ifail = 4;
                        return;
                }

                det_pjil = l_p[1] * Sij_1 - l_p[2] * Sij_2 + l_p[3] * Sij_3;
                det_pjil = det_pjil * iorient;
                test_pjil = (fabs(det_pjil) > eps);
                if(test_pjil && det_pjil > 0)
                {
                        *ifail = 3;
                        return;
                }

                det_pkjl = j_p[1] * Skl_1 - j_p[2] * Skl_2 + j_p[3] * Skl_3;
                det_pkjl = det_pkjl * iorient;
                test_pkjl = (fabs(det_pkjl) > eps);
                if(test_pkjl && det_pkjl > 0)
                {
                        *ifail = 1;
                        return;
                }

                det_pikl = -i_p[1] * Skl_1 + i_p[2] * Skl_2 - i_p[3] * Skl_3;
                det_pikl = det_pikl * iorient;
                test_pikl = (fabs(det_pikl) > eps);
                if(test_pikl && det_pikl > 0)
                {
                        *ifail = 2;
                        return;
                }
                /*
                At this stage, either all four determinants are positive,
                or one of the determinant is not precise enough, and
                we need to switch the MP
                In this case, since we may need SoS, we have to rank
                the indices
                */
                if(!test_pijk)
                {
                        ValSort4(p, a, b, c, &ia, &ib, &ic, &id, &nswap);
                        sos->Minor4(vertexList,&ia,&ib,&ic,&id,&val);
                        val = val * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = 4;
                                return;
                        }
                }

                if(!test_pjil)
                {
                        ValSort4(p, b, a, d, &ia, &ib, &ic, &id, &nswap);
                        sos->Minor4(vertexList, &ia, &ib, &ic, &id, &val);
                        val = val * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = 3;
                                return;
                        }
                }

                if(!test_pkjl)
                {
                        ValSort4(p, c, b, d, &ia, &ib, &ic, &id, &nswap);
                        sos->Minor4(vertexList, &ia, &ib, &ic, &id, &val);
                        val = val * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = 1;
                                return;
                        }
                }

                if(!test_pikl)
                {
                        ValSort4(p, a, c, d, &ia, &ib, &ic, &id, &nswap);
                        sos->Minor4(vertexList, &ia, &ib, &ic, &id, &val);
                        val = val * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = 2;
                                return;
                        }
                }

                //If we have gone that far, p is inside the tetrahedron
                *is_in = 1;//true

                //Now we check if p is redundant
                i_p[4] = vertexList[a].Weight - vertexList[p].Weight;
                j_p[4] = vertexList[b].Weight - vertexList[p].Weight;
                k_p[4] = vertexList[c].Weight - vertexList[p].Weight;
                l_p[4] = vertexList[d].Weight - vertexList[p].Weight;

                det_pijkl = -i_p[4] * det_pkjl - j_p[4] * det_pikl - k_p[4] * det_pjil - l_p[4] * det_pijk;

                //No need to multiply by iorient, since all minors contais iorient
                if(fabs(det_pijkl) < eps)
                {
                        ValSort5(p, a, b, c, d, &ia, &ib, &ic, &id, &ie, &nswap);
                        sos->Minor5(vertexList, &ia, &ib, &ic, &id, &ie, &val);
                        det_pijkl = val * nswap * iorient;
                }
                *redundant = (det_pijkl < 0)?1:0;
        }
        else if (ninf == 1)
        {       /*
                We know that one of the 4 vertices a,b,c or d is infinite
                To find which one it is, we use a map between
                (inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
                is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4
                if a,b,c or d are infinite, respectively.
                A good mapping function is:
                X = 3 - inf(a) - inf(a) -inf(b) + inf(d)
                */
                idx = (3 - infPoint[a] - infPoint[a] - infPoint[b] + infPoint[d]);
                l = list[idx] - npoints;

                i = list[order1[idx][1]];
                j = list[order1[idx][2]];
                k = list[order1[idx][3]];

                ic1 = inf4_1[l];
                sign = sign4_1[l];
                /*
                let us look at the four determinant we need to compute:
                det_pijk	: unchanged
                det_pjil	: 1 infinite point (l), becomes det3_pji where
                det3_pij =        |p(ic1) p(ic2) 1|
                                  |i(ic1) i(ic2) 1|
                                  |j(ic1) j(ic2) 1|
                and ic1 and ic2 depends on which infinite ( ic2 is always 3)
                point is considered
                det_pkjl	: 1 infinite point (l), becomes det3_pkj
                det_pikl	: 1 infinite point (l), becomes det3_pik
                */
                //get Coordinates
                for (int _i = 1; _i < 4; _i++)
                {
                        i_p[_i] = vertexList[i].Coordinates[_i] - vertexList[p].Coordinates[_i];
                        j_p[_i] = vertexList[j].Coordinates[_i] - vertexList[p].Coordinates[_i];
                        k_p[_i] = vertexList[k].Coordinates[_i] - vertexList[p].Coordinates[_i];
                }

                detij[1] = i_p[1] * j_p[3] - i_p[3] * j_p[1];
                detij[2] = i_p[2] * j_p[3] - i_p[3] * j_p[2];
                detij[3] = i_p[1] * j_p[2] - i_p[2] * j_p[1];

                *is_in = 0;//false

                det_pijk = -k_p[1] * detij[2] + k_p[2] * detij[1] - k_p[3] * detij[3];
                det_pijk = det_pijk * iorient;
                test_pijk = (fabs(det_pijk) > eps);
                if(test_pijk && det_pijk > 0)
                {
                        *ifail = idx;
                        return;
                }

                det_pjil = -detij[ic1] * sign * iorient;
                test_pjil = (fabs(det_pjil) > eps);
                if(test_pjil && det_pjil > 0)
                {
                        *ifail = order1[idx][3];
                        return;
                }

                det_pkjl = k_p[ic1] * j_p[3] - k_p[3] * j_p[ic1];
                det_pkjl = sign * det_pkjl * iorient;
                test_pkjl = (fabs(det_pkjl) > eps);
                if(test_pkjl && det_pkjl > 0)
                {
                        *ifail = order1[idx][1];
                        return;
                }

                det_pikl = i_p[ic1] * k_p[3] - i_p[3] * k_p[ic1];
                det_pikl = sign * det_pikl * iorient;
                test_pikl = (fabs(det_pikl) > eps);
                if(test_pikl && det_pikl > 0)
                {
                        *ifail = order1[idx][2];
                        return;
                }
                /*
                At this stage, either all four determinants are positive,
                or one of the determinant is not precise enough, and
                we need to switch the MP
                */
                if(!test_pijk)
                {
                        ValSort4(p, i, j, k, &ia, &ib, &ic, &id, &nswap);
                        sos->Minor4(vertexList, &ia, &ib, &ic, &id, &val);
                        val = val * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = idx;
                                return;
                        }
                }

                if(!test_pjil)
                {
                        ValSort3(p, j, i, &ia, &ib, &ic, &nswap);
                        int temp = 3;
                        sos->Minor3(vertexList, &ia, &ib, &ic, &ic1, &temp, &val);
                        val = val * sign * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = order1[idx][3];
                                return;
                        }
                }

                if(!test_pkjl)
                {
                        ValSort3(p, k, j, &ia, &ib, &ic, &nswap);
                        int temp = 3;
                        sos->Minor3(vertexList, &ia, &ib, &ic, &ic1, &temp, &val);
                        val = val * sign * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = order1[idx][1];
                                return;
                        }
                }

                if(!test_pikl)
                {
                        ValSort3(p, i, k, &ia, &ib, &ic, &nswap);
                        int temp = 3;
                        sos->Minor3(vertexList, &ia, &ib, &ic, &ic1, &temp, &val);
                        val = val * sign * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = order1[idx][2];
                                return;
                        }
                }

                //If we have gone so far, p is inside the tetrahedron
                *is_in = 1;//true
                /*
                Now we check if p is redundant
                since det_pijkl = det_pijk > 1
                p cannot be redundant !
                */
                *redundant = 0;//false
        }
        else if(ninf == 2)
        {	/*
                We know that two of the 4 vertices a,b,c or d are infinite
                To find which one it is, we use a map between
                (inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
                is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4,5,6
                if (a,b), (a,c), (a,d), (b,c), (b,d), or (c,d) are
                infinite, respectively
                A good mapping function is:
                X = 3 - inf(a) - inf(a) +inf(c) + inf(d) + inf(d)
                */
                idx = (3 -infPoint[a] - infPoint[a] + infPoint[c] + 2 * infPoint[d]);

                //The two infinite points :
                k = list[order3[idx][1]] - npoints;
                l = list[order3[idx][2]] - npoints;

                //The two finite points
                i = list[order2[idx][1]];
                j = list[order2[idx][2]];

                ic1_k = inf4_1[k];
                ic1_l = inf4_1[l];
                sign_k = sign4_1[k];
                sign_l = sign4_1[l];
                ic1 = inf4_2[l][k];
                sign = sign4_2[l][k];

                //Get coordinates
                for (int _i = 1; _i < 4; _i++)
                {
                        i_p[_i] = vertexList[i].Coordinates[_i] - vertexList[p].Coordinates[_i];
                        j_p[_i] = vertexList[j].Coordinates[_i] - vertexList[p].Coordinates[_i];
                }

                //Perform test; first set is_in .false.
                *is_in = 0;//false

                //det_pijk is now det3_pij with k as infinite point
                det_pijk = i_p[ic1_k] * j_p[3] - i_p[3] * j_p[ic1_k];
                det_pijk = det_pijk * sign_k * iorient;
                test_pijk = (fabs(det_pijk) > eps);
                if(test_pijk && det_pijk > 0)
                {
                        *ifail = order3[idx][2];
                        return;
                }

                //det_pjil is now det3_pji with l as infinite point
                det_pjil = i_p[3] * j_p[ic1_l] - i_p[ic1_l] * j_p[3];
                det_pjil = det_pjil * sign_l * iorient;
                test_pjil = (fabs(det_pjil) > eps);
                if(test_pjil && det_pjil > 0)
                {
                        *ifail = order3[idx][1];
                        return;
                }

                //det_pkjl is now -det2_pj (k,l infinite)
                det_pkjl = j_p[ic1] * sign * iorient;
                test_pkjl = (fabs(det_pkjl) > eps);
                if(test_pkjl && det_pkjl > 0)
                {
                        *ifail = order2[idx][1];
                        return;
                }

                //det_pikl is now det2_pi (k,l infinite)
                det_pikl = -i_p[ic1] * sign * iorient;
                test_pikl = (fabs(det_pikl) > eps);
                if(test_pikl && det_pikl > 0)
                {
                        *ifail = order2[idx][2];
                        return;
                }
                /*
                At this stage, either all four determinants are positive,
                or one of the determinant is not precise enough, and
                we need to switch the MP
                */
                if(!test_pijk)
                {
                        ValSort3(p, i, j, &ia, &ib, &ic, &nswap);
                        int temp = 3;
                        sos->Minor3(vertexList, &ia, &ib, &ic, &ic1_k, &temp, &val);
                        val = val * sign_k * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = order3[idx][2];
                                return;
                        }
                }

                if(!test_pjil)
                {
                        ValSort3(p, j, i, &ia, &ib, &ic, &nswap);
                        int temp = 3;
                        sos->Minor3(vertexList, &ia, &ib, &ic, &ic1_l, &temp, &val);
                        val = val * sign_l * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = order3[idx][1];
                                return;
                        }
                }

                if(!test_pkjl)
                {
                        ValSort2(p, j, &ia, &ib, &nswap);
                        sos->Minor2(vertexList,&ia,&ib,ic1,&val);
                        val = -val * sign * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = order2[idx][1];
                                return;
                        }
                }

                if(!test_pikl)
                {
                        ValSort2(p, i, &ia, &ib, &nswap);
                        sos->Minor2(vertexList,&ia,&ib,ic1,&val);
                        val = val * sign * nswap * iorient;
                        if(val == 1)
                        {
                                *ifail = order2[idx][2];
                                return;
                        }
                }

                //Again, if we have gone so far, p is inside the tetrahedron
                *is_in = 1;//true

                //Now we check if p is redundant
                //det_pijkl becomes det3_pij
                ic5 = inf5_2[l][k];
                sign5 = sign5_2[l][k];
                det_pijkl = i_p[ic5] * j_p[3] - i_p[3] * j_p[ic5];
                if(fabs(det_pijkl) < eps)
                {
                        ValSort3(p, i, j, &ia, &ib, &ic, &nswap);
                        int temp = 3;
                        sos->Minor3(vertexList,&ia,&ib,&ic,&ic5,&temp,&val);
                        det_pijkl = val*nswap;
                }
                det_pijkl = det_pijkl * sign5 * iorient;

                *redundant = (det_pijkl < 0)?1:0;
        }
        else if (ninf == 3)
        {       /*
                We know that three of the 4 vertices a,b,c or d are infinite
                To find which one is finite, we use a map between
                (inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
                is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4
                if a,b,c or d are finite, respectively.
                A good mapping function is:
                X = 1 + inf(a) + inf(a) +inf(b) - inf(d)
                */
                idx = (1 + 2 * infPoint[a] + infPoint[b] - infPoint[d]);
                i = list[idx];
                j = list[order1[idx][1]] - npoints;
                k = list[order1[idx][2]] - npoints;
                l = list[order1[idx][3]] - npoints;

                //Index of the "missing" infinite point (i.e. the fourth infinite point)
                MissInfSign(j,k,l,&ie,&iswap);

                //Get coordinates
                for (int _i = 1; _i < 4; _i++)
                {
                        i_p[_i] = vertexList[i].Coordinates[_i] - vertexList[p].Coordinates[_i];
                }

                //Perform test; first set is_in to .false.
                *is_in = 0;//false

                //det_pijk is now - det2_pi (missing j,k)
                det_pijk = i_p[inf4_2[k][j]] * iorient * sign4_2[k][j];
                test_pijk = (fabs(det_pijk) > eps);
                if(test_pijk && det_pijk > 0)
                {
                        *ifail = order1[idx][3];
                        return;
                }

                //det_pjil is now det2_pi (missing j,l)
                det_pjil = -i_p[inf4_2[l][j]] * iorient * sign4_2[l][j];
                test_pjil = (fabs(det_pjil) > eps);
                if(test_pjil && det_pjil > 0)
                {
                        *ifail = order1[idx][2];
                        return;
                }

                //det_pkjl is now det1_p
                det_pkjl = iorient * iswap * sign4_3[ie];
                if(det_pkjl > 0)
                {
                        *ifail = idx;
                        return;
                }

                //det_ikl is now - det2_pi (missing k,l)
                det_pikl = i_p[inf4_2[l][k]] * iorient * sign4_2[l][k];
                test_pikl = (fabs(det_pikl) > eps);
                if(test_pikl && det_pikl > 0)
                {
                        *ifail = order1[idx][1];
                        return;
                }
                /*
                At this stage, either all four determinants are positive,
                or one of the determinant is not precise enough, and
                we need to switch the MP
                */
                if(!test_pijk)
                {
                        ValSort2(p, i, &ia, &ib, &nswap);
                        sos->Minor2(vertexList,&ia,&ib,inf4_2[k][j],&val);
                        val = -val * sign4_2[k][j] * iorient * nswap;
                        if(val == 1)
                        {
                                *ifail = order1[idx][3];
                                return;
                        }
                }

                if(!test_pjil)
                {
                        ValSort2(p, i, &ia, &ib, &nswap);
                        sos->Minor2(vertexList,&ia,&ib,inf4_2[l][j],&val);
                        val = val * sign4_2[l][j] * iorient * nswap;
                        if(val == 1)
                        {
                                *ifail = order1[idx][2];
                                return;
                        }
                }

                if(!test_pikl)
                {
                        ValSort2(p, i, &ia, &ib, &nswap);
                        sos->Minor2(vertexList,&ia,&ib,inf4_2[l][k],&val);
                        val = -val * sign4_2[l][k] * iorient * nswap;
                        if(val == 1)
                        {
                                *ifail = order1[idx][1];
                                return;
                        }
                }

                *is_in = 1;//true

                //Now check for redundancy
                //det_pijkl becomes -det2_pi

                ic1 = inf5_3[ie];
                sign5 = sign5_3[ie];
                det_pijkl = -i_p[ic1];
                if(fabs(det_pijkl) < eps)
                {
                        ValSort2(p, i, &ia, &ib, &nswap);
                        sos->Minor2(vertexList,&ia,&ib,ic1,&val);
                        det_pijkl = val * nswap;
                }
                det_pijkl = - iorient * det_pijkl * sign5 * iswap;
                *redundant = (det_pijkl < 0)?1:0;

        }
        else
        {
                //In the case all four points ia,ib,ic, and id are infinite,
                //then is_in = .true. and redundant = .false.
                *is_in = 1;//true
                *redundant = 0;//false
        }
}

/*!
    \fn DeluanayComplex::LocateByJumpWalk(std::vector<Vertex> & vertexList,int *iseed, int ival,int *tetraLoc,int *iredundant)
 */
void DeluanayComplex::LocateByJumpWalk(std::vector<Vertex> & vertexList,int *iseed, int ival,int *tetraLoc,int *iredundant)
{
        int itetra = -1, iorient, intest = 0, redtest = 0, idx = -1;
        int a, b, c, d;

        *iredundant = 0;

        if (ntetra == 1)
        {
                *tetraLoc = 1;
                return;
        }

        Jump(vertexList,iseed, ival, &itetra);

        while(true)
        {

                a = DeluanayTet[itetra].Corners[1];
                b = DeluanayTet[itetra].Corners[2];
                c = DeluanayTet[itetra].Corners[3];
                d = DeluanayTet[itetra].Corners[4];

                iorient = DeluanayTet[itetra].Orient;

                InsideTetraJW(vertexList,ival, a, b, c, d, iorient, &intest, &redtest, &idx);

                if (intest != 0)
                {
                        break;
                }

                itetra = DeluanayTet[itetra].Neighbours[idx];
        }

        *tetraLoc = itetra;
        if (redtest != 0)
        {
                *iredundant = 1;
        }
}

/*!
    \fn DeluanayComplex::ConstructDT(std::vector<Vertex> & vertexList)
 */
void DeluanayComplex::ConstructDT(std::vector<Vertex> & vertexList)
{
        //FILE *fp;
        for (uint i = 0; i < vertexList.size() + 4; i++)
        {
                infPoint.push_back(0);
        }
        //add four extra points for dummy tetrahedra and mark them as infpoint
        infPoint[vertexList.size()] = infPoint[vertexList.size() + 1] = infPoint[vertexList.size() + 2] = infPoint[vertexList.size() + 3] = 1;

        //sort ranlist along x axis and eliminate duplicates if any
        std::vector<Vertex> sortedVertexList(vertexList);
        std::sort(sortedVertexList.begin(),sortedVertexList.end());

        EliminateDup(sortedVertexList);

        for(uint i = 1;i < vertexList.size();i++)
        {
                if(sortedVertexList[i].Redinfo == 1)
                {
                        vertexList[sortedVertexList[i].Index].Redinfo = 1;
                }
        }

        //initialise facet structures with a dummy entry
        AddDummyLF();
        // add four infinite points and initialize first tetrahedron
        nvertex = vertexList.size() + 4;

        ntetra = 1;

        Tetrahedron tet;
        DeluanayTet.push_back(tet);
        //dummy tet just for index shift
        tet.Corners[1] = vertexList.size() + 0;
        tet.Corners[2] = vertexList.size() + 1;
        tet.Corners[3] = vertexList.size() + 2;
        tet.Corners[4] = vertexList.size() + 3;

        for (int i = 1; i <= 4; i++)
        {
                tet.Neighbours[i] = 0;
        }

        tet.Status = 1;
        tet.Orient = -1;
        //infinite tetrahedron
        DeluanayTet.push_back(tet);
        InitFreeKill();

        //randomise the way in which points are added
        iseed = -1;

        for (uint i = 1; i < vertexList.size(); i++)
        {
                vertexList[i].ranValue = ran2(&iseed);
        }

        std::vector<Vertex> randomSortedList(vertexList);
        std::sort(randomSortedList.begin()+1,randomSortedList.end(),Vertex::ranSortFunc);

        //fp = fopen("tetnum.txt","w");
        //now add points one by one
        for (uint i = 1; i < vertexList.size(); i++)
        {
                ival = randomSortedList[i].Index;

                if (vertexList[ival].Redinfo == 1)
                {
                        //avoid redundant points
                        continue;
                }
                //first locate the point in the list of known tetrahedra
                LocateByJumpWalk(vertexList, &iseed, ival, &tetraLoc, &iredundant);
                //if point is redundant move to next point
                if (iredundant == 1)
                {
                        vertexList[ival].Redinfo = 1;
                        continue;
                }
                //otherwise add to the tetrahedron given by position tetraLoc and perform 1-4 flip
                FlipJW1_4(ival, tetraLoc);
                //scan link facet list and flip till the list is empty
                FlipJW(vertexList);
                //fprintf(fp,"%d\n",ival);
        }
        //fclose(fp);
        //after triangulation remove infinite points and define convex hull
        RemoveInf();
        //peel of flat tetrahedra
        Peel(vertexList);
}

/*!
    \fn DeluanayComplex::DefineTriangles()
 */
void DeluanayComplex::DefineTriangles()
{
        uint i,idx;
        int j, k, l;
        int jtetra, jindex, nswap = 0, ntrig = 0;
        int vertice[5];
        int neighbor[5];
        int index[5];
        int indx[5];
        int trig1, trig2, trig3, trig4;

        //as usual add dummy triangle first
        Triangle tri;
        DeluanayTrigs.push_back(tri);

        //loop over all tets to generate triangles
        for (i = 1; i < DeluanayTet.size(); i++)
        {
                //ignore dead tets
                if (DeluanayTet[i].Status == 0)
                {
                        redundantCount++;
                        continue;
                }

                //first sort corners,neighbours etc
                for (j = 1; j <= 4; j++)
                {
                        vertice[j] = DeluanayTet[i].Corners[j];
                }
                int temp = 4;
                ISortIdx(vertice, indx, &nswap,&temp);

                for (j = 1; j <= 4; j++)
                {
                        neighbor[j] = DeluanayTet[i].Neighbours[indx[j]];
                        index[j] = DeluanayTet[i].Nindex[indx[j]];
                }

                if (nswap % 2 == 1)
                {
                    DeluanayTet[i].Orient = -DeluanayTet[i].Orient;
                }

                for (j = 1; j <= 4; j++)
                {
                        DeluanayTet[i].Corners[j] = vertice[j];
                        DeluanayTet[i].Neighbours[j] = neighbor[j];
                        DeluanayTet[i].Nindex[j] = index[j];
                        if (neighbor[j] != 0)
                        {
                                DeluanayTet[neighbor[j]].Nindex[index[j]] = j;
                        }
                }

                //set hull to 0
                DeluanayTet[i].Hull = 0;

                //Generate all triangles
                for (j = 1; j <= 4; j++)
                {
                        jtetra = DeluanayTet[i].Neighbours[j];

                        if ((jtetra == 0) || (jtetra > i))
                        {
                                //first add dummy triangle fill later;
                                Triangle tri;
                                DeluanayTrigs.push_back(tri);

                                l = 0;
                                ntrig++;
                                for (k = 1; k <= 4; k++)
                                {
                                        if (k == j) continue;
                                        l++;
                                        DeluanayTrigs[ntrig].Corners[l] = DeluanayTet[i].Corners[k];
                                }
                                if (jtetra == 0)
                                {
                                        DeluanayTet[i].Hull = 1;
                                        DeluanayTrigs[ntrig].Hull = 1;
                                }
                                else
                                {
                                        DeluanayTrigs[ntrig].Hull = 0;
                                }
                                DeluanayTet[i].TetLink[j] = ntrig;
                        }
                        else
                        {
                                jindex = DeluanayTet[i].Nindex[j];
                                DeluanayTet[i].TetLink[j] = DeluanayTet[jtetra].TetLink[jindex];
                        }
                }
        }

        for(idx=1;idx<DeluanayTet.size();idx++)
        {
                if(DeluanayTet[idx].Status == 0)
                        continue;

                trig1 = DeluanayTet[idx].TetLink[1];
                trig2 = DeluanayTet[idx].TetLink[2];
                trig3 = DeluanayTet[idx].TetLink[3];
                trig4 = DeluanayTet[idx].TetLink[4];

                DeluanayTrigs[trig1].nLink++;
                if(DeluanayTrigs[trig1].nLink == 1)
                {
                        DeluanayTrigs[trig1].ReverseLink1 = idx;
                }
                else if(DeluanayTrigs[trig1].nLink == 2)
                {
                        DeluanayTrigs[trig1].ReverseLink2 = idx;
                }

                DeluanayTrigs[trig2].nLink++;
                if(DeluanayTrigs[trig2].nLink == 1)
                {
                        DeluanayTrigs[trig2].ReverseLink1 = idx;
                }
                else if(DeluanayTrigs[trig2].nLink == 2)
                {
                        DeluanayTrigs[trig2].ReverseLink2 = idx;
                }

                DeluanayTrigs[trig3].nLink++;
                if(DeluanayTrigs[trig3].nLink == 1)
                {
                        DeluanayTrigs[trig3].ReverseLink1 = idx;
                }
                else if(DeluanayTrigs[trig3].nLink == 2)
                {
                        DeluanayTrigs[trig3].ReverseLink2 = idx;
                }

                DeluanayTrigs[trig4].nLink++;
                if(DeluanayTrigs[trig4].nLink == 1)
                {
                        DeluanayTrigs[trig4].ReverseLink1 = idx;
                }
                else if(DeluanayTrigs[trig4].nLink == 2)
                {
                        DeluanayTrigs[trig4].ReverseLink2 = idx;
                }
        }
        /*FILE *fp = fopen("reverselink.txt","w");
        for(i=1;i<DeluanayTrigs.size();i++)
        {
                fprintf(fp,"nlink = %d ",DeluanayTrigs[i].nLink);
                if (DeluanayTrigs[i].nLink == 1)
                {
                        fprintf(fp,"reverselink1 = %d\n",DeluanayTrigs[i].ReverseLink1);
                }
                else if(DeluanayTrigs[i].nLink == 2)
                {
                        fprintf(fp,"reverselink1 = %d reverselink2 = %d\n",DeluanayTrigs[i].ReverseLink1,DeluanayTrigs[i].ReverseLink2);
                }
                else
                {
                        fprintf(fp,"\n");
                }
        }
        fclose(fp);*/
}

/*!
    \fn DeluanayComplex::DefineEdges(std::vector<Vertex> & vertexList)
 */
void DeluanayComplex::DefineEdges(std::vector<Vertex> & vertexList)
{
        uint i, idx;
        int ipos, jpos;
        int a, b, imask;
        int trig1, trig2, triga, trigb, ipair;
        int itrig1, itrig2;
        int ktetra, trig_in, trig_out, npass, itetra;
        int nedge;
        int flag = 0;

        std::vector<int> tetra_mask;

        for (i = 0; i < DeluanayTet.size(); i++)
        {
                tetra_mask.push_back(0);
        }

        nedge = 0;
        Edge edge;
        DeluanayEdges.push_back(edge);

        for (idx = 1; idx < DeluanayTet.size(); idx++)
        {
                //ignore dead tets as usual
                if (DeluanayTet[idx].Status == 0) continue;

                imask = tetra_mask[idx];

                for (i = 1; i <= 6; i++)
                {
                        if ((imask & mask[i]) != 0) continue;

                        a = DeluanayTet[idx].Corners[pair[i][1]];
                        b = DeluanayTet[idx].Corners[pair[i][2]];

                        nedge++;
                        DeluanayEdges.push_back(edge);
                        DeluanayEdges[nedge].Corners[1] = a;
                        DeluanayEdges[nedge].Corners[2] = b;

                        trig1 = trig_info[i][1];
                        trig2 = trig_info[i][2];

                        itrig1 = DeluanayTet[idx].TetLink[trig1];
                        itrig2 = DeluanayTet[idx].TetLink[trig2];

                        ipos = trig_pos[i][1];
                        jpos = trig_pos[i][2];

                        DeluanayTrigs[itrig1].TrigLink[ipos] = nedge;
                        DeluanayTrigs[itrig2].TrigLink[jpos] = nedge;

                        ktetra = idx;
                        npass = 1;
                        trig_out = trig1;
                        /*
                        Loop over all the tetrahedra in the link of the edge: start in one direction, and if we
                        hit the convex hull before cycling back to the initial tetrahedron, restart in the other
                        direction
                        */
                        while(true)
                        {
                                itetra = DeluanayTet[ktetra].Neighbours[trig_out];

                                //Leave this side of the link if we hit	the convex hull
                                if (itetra == 0)
                                {
                                        if (npass != 2)
                                        {
                                                npass++;
                                                ktetra = idx;
                                                trig_out = trig2;
                                                continue;
                                        }
                                        else
                                        {
                                                flag = 1;
                                                break;
                                        }
                                }

                                //Leave the loop completely if we have described the full cycle
                                if (itetra == idx) break;

                                //Identify the position of edge (ab) in	tetrahedron itetra
                                if (a == DeluanayTet[itetra].Corners[1])
                                {
                                        if (b == DeluanayTet[itetra].Corners[2])
                                        {
                                                ipair = 6;
                                        }
                                        else if (b == DeluanayTet[itetra].Corners[3])
                                        {
                                                ipair = 5;
                                        }
                                        else
                                        {
                                                ipair = 4;
                                        }
                                }
                                else if (a == DeluanayTet[itetra].Corners[2])
                                {
                                        if (b == DeluanayTet[itetra].Corners[3])
                                        {
                                                ipair = 3;
                                        }
                                        else
                                        {
                                                ipair = 2;
                                        }
                                }
                                else
                                {
                                        ipair = 1;
                                }

                                tetra_mask[itetra] += mask[ipair];
                                trig_in = DeluanayTet[ktetra].Nindex[trig_out];

                                triga = trig_info[ipair][1];
                                trigb = trig_info[ipair][2];

                                trig_out = triga;

                                ipos = trig_pos[ipair][1];
                                if (trig_in == triga)
                                {
                                        trig_out = trigb;
                                        ipos = trig_pos[ipair][2];
                                }

                                DeluanayTrigs[DeluanayTet[itetra].TetLink[trig_out]].TrigLink[ipos] = nedge;

                                ktetra = itetra;
                        }
                        if (flag == 1)
                        {
                                flag = 0;
                                continue;
                        }
                }
        }

        for (i = 1; i < DeluanayTrigs.size(); i++)
        {
                if (DeluanayTrigs[i].Hull == 1)
                {
                        DeluanayEdges[DeluanayTrigs[i].TrigLink[1]].Hull = 1;
                        DeluanayEdges[DeluanayTrigs[i].TrigLink[2]].Hull = 1;
                        DeluanayEdges[DeluanayTrigs[i].TrigLink[3]].Hull = 1;
                        vertexList[DeluanayTrigs[i].Corners[1]].Hull = 1;
                        vertexList[DeluanayTrigs[i].Corners[2]].Hull = 1;
                        vertexList[DeluanayTrigs[i].Corners[3]].Hull = 1;
                }
        }
}

void DeluanayComplex::CalculateNormals(std::vector<Vertex> &vertexList)
{
        uint i;
        int a,b,c,ia = 0, ib = 0, ic = 0;
        mpz_t a1,a2,a3,b1,b2,b3,nx,ny,nz,s,temp1,temp2;
        mpz_init(a1);mpz_init(a2);mpz_init(a3);
        mpz_init(b1);mpz_init(b2);mpz_init(b3);
        mpz_init(nx);mpz_init(ny);mpz_init(nz);
        mpz_init(s);mpz_init(temp1);mpz_init(temp2);

        for (i = 1; i < DeluanayTrigs.size(); i++)
        {
                a = DeluanayTrigs[i].Corners[1];
                b = DeluanayTrigs[i].Corners[2];
                c = DeluanayTrigs[i].Corners[3];

                //find a,b,c are ccw
                int nswaps = 0;
                int i1 = 1,i2 = 2,r = 0;
                ValSort3(a,b,c,&ia,&ib,&ic,&nswaps);
                sos->Minor3(vertexList,&ia,&ib,&ic,&i1,&i2,&r);

                r = r*nswaps;

                //r > 0 => ccw
                mpz_sub(a1,vertexList[b].V[1],vertexList[a].V[1]);
                mpz_sub(a2,vertexList[b].V[2],vertexList[a].V[2]);
                mpz_sub(a3,vertexList[b].V[3],vertexList[a].V[3]);

                mpz_sub(b1,vertexList[c].V[1],vertexList[a].V[1]);
                mpz_sub(b2,vertexList[c].V[2],vertexList[a].V[2]);
                mpz_sub(b3,vertexList[c].V[3],vertexList[a].V[3]);

                mpz_mul(temp1,a2,b3);
                mpz_mul(temp2,a3,b2);
                mpz_sub(nx,temp1,temp2);

                mpz_mul(temp1,a3,b1);
                mpz_mul(temp2,a1,b3);
                mpz_sub(ny,temp1,temp2);

                mpz_mul(temp1,a1,b2);
                mpz_mul(temp2,a2,b1);
                mpz_sub(nz,temp1,temp2);

                mpz_mul(s,nx,nx);
                mpz_mul(temp1,ny,ny);
                mpz_mul(temp2,nz,nz);

                mpz_add(s,s,temp1);
                mpz_add(s,s,temp2);

                double length = mpz_get_d(s);
                length = sqrt(length);

                DeluanayTrigs[i].Normal->X = (r > 0) ? (mpz_get_d(nx)/length) : (-mpz_get_d(nx)/length);
                DeluanayTrigs[i].Normal->Y = (r > 0) ? (mpz_get_d(ny)/length) : (-mpz_get_d(ny)/length);
                DeluanayTrigs[i].Normal->Z = (r > 0) ? (mpz_get_d(nz)/length) : (-mpz_get_d(nz)/length);

                //DeluanayTrigs[i].Normal->X = mpz_get_d(nx)/length;
                //DeluanayTrigs[i].Normal->Y = mpz_get_d(ny)/length;
                //DeluanayTrigs[i].Normal->Z = mpz_get_d(nz)/length;
        }
        mpz_clear(a1);mpz_clear(a2);mpz_clear(a3);
        mpz_clear(b1);mpz_clear(b2);mpz_clear(b3);
        mpz_clear(nx);mpz_clear(ny);mpz_clear(nz);
        mpz_clear(s);mpz_clear(temp1);mpz_clear(temp2);
}

/*!
    \fn DeluanayComplex::CNInit()
 */
void DeluanayComplex::CNInit()
{
        int i;
        Nvert = nvertex - 4;
        NTri = DeluanayTrigs.size();
        for(i=0;i<=NTri;i++)
        {
                mark_buf.push_back(0);
                tri_buf.push_back(0);
        }
}

/*!
    \fn DeluanayComplex::CNClean()
 */
void DeluanayComplex::CNClean()
{
        mark_buf.clear();
        tri_buf.clear();
}

/*!
    \fn DeluanayComplex::MarkCH(std::vector<Vertex> &vertexList)
 */
int DeluanayComplex::MarkCH(std::vector<Vertex> &vertexList)
{
        Vector3 *v1 = new Vector3();
        Vector3 *v2 = new Vector3();

        uint t;
        int i;
        tcount = 0;
        for(t = 1; t < DeluanayTrigs.size(); t++)
        {
                if(DeluanayTrigs[t].Hull == 1)
                {
                        double res = 0.0;
                        i = DeluanayTrigs[t].Corners[1];

                        v1->X = vertexList[i].Coordinates[1];
                        v1->Y = vertexList[i].Coordinates[2];
                        v1->Z = vertexList[i].Coordinates[3];

                        v2->X = DeluanayTrigs[t].Normal->X;
                        v2->Y = DeluanayTrigs[t].Normal->Y;
                        v2->Z = DeluanayTrigs[t].Normal->Z;

                        Vector3::DotProduct(v1,v2,&res);

                        if(res<0)
                        {
//                                v2->X = -v2->X;
//                                v2->Y = -v2->Y;
//                                v2->Z = -v2->Z;
                                DeluanayTrigs[t].Normal->X = -v2->X;
                                DeluanayTrigs[t].Normal->Y = -v2->Y;
                                DeluanayTrigs[t].Normal->Z = -v2->Z;
                        }
                        mark_buf[t] = 1;
                        tri_buf[tcount] = t;
                        tcount++;
                }
        }
        return 1;
}

/*!
    \fn DeluanayComplex::MarkNeighbours(int t,int n_iter)
 */
void DeluanayComplex::MarkNeighbours(int t,int n_iter)
{
        double res;

        std::vector <int> trig;

        if(DeluanayTrigs[t].nLink == 2)
        {
                int tet1 = DeluanayTrigs[t].ReverseLink1;
                int tet2 = DeluanayTrigs[t].ReverseLink2;

                for(int i = 1; i < 5; i++)
                {
                        if(DeluanayTet[tet1].TetLink[i] != t)
                        {
                                trig.push_back(DeluanayTet[tet1].TetLink[i]);
                        }
                        if(DeluanayTet[tet2].TetLink[i] != t)
                        {
                                trig.push_back(DeluanayTet[tet2].TetLink[i]);
                        }
                }
        }
        else if(DeluanayTrigs[t].nLink == 1)
        {
                int tet1 = DeluanayTrigs[t].ReverseLink1;
                for(int i = 1; i < 5; i++)
                {
                        if(DeluanayTet[tet1].TetLink[i] != t)
                        {
                                trig.push_back(DeluanayTet[tet1].TetLink[i]);
                        }
                }
        }

        for(uint i = 0;i< trig.size();i++)
        {
                int t1 = trig[i];

                if(!mark_buf[t1])
                {
                        res = 0.0;
                        Vector3::DotProduct(DeluanayTrigs[t].Normal,DeluanayTrigs[t1].Normal,&res);
                        if(res < -0.9)
                        {
                                DeluanayTrigs[t1].Normal->X *= -1;
                                DeluanayTrigs[t1].Normal->Y *= -1;
                                DeluanayTrigs[t1].Normal->Z *= -1;
                                res*=-1;
                        }
                        if(res < 0.9)
                        {
                                mark_buf[t1] = n_iter;
                                tri_buf[tcount] = t1;
                                tcount++;
                        }
                }
        }
        trig.clear();
}

/*!
    \fn DeluanayComplex::MarkNeighboursInSet(int n_iter,int i_beg,int i_end)
 */
void DeluanayComplex::MarkNeighboursInSet(int n_iter,int i_beg,int i_end)
{
        int i;
        for(i=i_beg;i<i_end;i++)
        {
                MarkNeighbours(tri_buf[i], n_iter);
        }
}

/*!
    \fn DeluanayComplex::MarkIterate(std::vector<Vertex> &vertexList)
 */
void DeluanayComplex::MarkIterate(std::vector<Vertex> &vertexList)
{
        int n_iter,i_beg,i_end;
        n_iter = MarkCH(vertexList);
        i_beg = 0;i_end = tcount;
        while((tcount<NTri) && (i_beg < i_end))
        {
                n_iter++;
                MarkNeighboursInSet(n_iter,i_beg,i_end);
                i_beg = i_end;i_end = tcount;
        }
}

/*!
    \fn DeluanayComplex::MarkNotMarked()
 */
void DeluanayComplex::MarkNotMarked()
{
        uint t;

        for(t = 1;t<DeluanayTrigs.size();t++)
        {
                if(!mark_buf[t])
                {
                        if(DeluanayTrigs[t].Normal->Z != 0)
                        {
                                if(DeluanayTrigs[t].Normal->Z < 0)
                                {
                                        DeluanayTrigs[t].Normal->X = -DeluanayTrigs[t].Normal->X;
                                        DeluanayTrigs[t].Normal->Y = -DeluanayTrigs[t].Normal->Y;
                                        DeluanayTrigs[t].Normal->Z = -DeluanayTrigs[t].Normal->Z;
                                }
                        }
                        else if(DeluanayTrigs[t].Normal->Y != 0)
                        {
                                if(DeluanayTrigs[t].Normal->Y < 0)
                                {
                                        DeluanayTrigs[t].Normal->X = -DeluanayTrigs[t].Normal->X;
                                        DeluanayTrigs[t].Normal->Y = -DeluanayTrigs[t].Normal->Y;
                                        DeluanayTrigs[t].Normal->Z = -DeluanayTrigs[t].Normal->Z;
                                }
                        }
                        else if(DeluanayTrigs[t].Normal->X != 0)
                        {
                                if(DeluanayTrigs[t].Normal->X < 0)
                                {
                                        DeluanayTrigs[t].Normal->X = -DeluanayTrigs[t].Normal->X;
                                        DeluanayTrigs[t].Normal->Y = -DeluanayTrigs[t].Normal->Y;
                                        DeluanayTrigs[t].Normal->Z = -DeluanayTrigs[t].Normal->Z;
                                }
                        }
                        mark_buf[t] = -1;
                }
        }
}

/*!
    \fn DeluanayComplex::CorrectNormals(std::vector<Vertex> &vertexList)
 */
void DeluanayComplex::CorrectNormals(std::vector<Vertex> &vertexList)
{
        CNInit();
        MarkIterate(vertexList);
        MarkNotMarked();
        CNClean();
}
