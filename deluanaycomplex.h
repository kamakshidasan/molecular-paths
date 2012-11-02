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

#ifndef DELUANAYCOMPLEX_H
#define DELUANAYCOMPLEX_H

#include <vector>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <tetrahedron.h>
#include <triangle.h>
#include <edge.h>
#include <vertex.h>
#include <sos.h>
#include <gmp.h>
#include <vector3.h>

class LinkFacet
{
        public:
        int f1;
        int f2;
};
class LinkIndex
{
        public:
        int i1;
        int i2;
};

class DeluanayComplex
{
    private:
        static int inf4_1[5];
        static int sign4_1[5];
        static int inf4_2[5][5];
        static int sign4_2[5][5];
        static int sign4_3[5];
        static int inf5_2[5][5];
        static int sign5_2[5][5];
        static int inf5_3[5];
        static int sign5_3[5];
        static int order1[5][4];
        static int order2[7][3];
        static int order3[7][3];
        static int idxList[5][4];
        static int table32[4][4];
        static int table32_2[4][3];
        static int table41[4][4];
        static int table41_2[4][3];
        static int order[4][3];
        static int order_1[4][4];
        static int other[5][4];
        static int other2[5][5][3];
        static int idx_list[4][3];
        static int mask[7];
        static int pair[7][3];
        static int trig_info[7][3];
        static int trig_pos[7][3];
        static double eps;
        static int ir[98];
        static int iy;

        int nvertex;
        int ntetra;
        int iseed;
        int nkill, nfree;
        std::vector <int> free;
        std::vector <int> kill;
        std::vector <LinkFacet> linkFacet;
        std::vector <LinkIndex> linkIndex;
        std::vector <int> infCount;
        int ival,iredundant,tetraLoc,iff;
        std::vector <int> infPoint;
        Sos *sos;

        int compare(Vertex *a,Vertex *b);

        void ValSort2(int a,int b,int *ia,int *ib,int *iswap);
        void ValSort3(int a,int b,int c,int *ia,int *ib,int *ic,int *iswap);
        void ValSort4(int a,int b,int c,int d,int *ia,int *ib,int *ic,int *id,int *iswap);
        void ValSort5(int a,int b,int c,int d,int e,int *ia,int *ib,int *ic,int *id,int *ie,int *iswap);

        void TetraVol(std::vector<Vertex>& vertexList, int ia, int ib, int ic, int id, double *vol);

        void ISortIdx(int vertice[],int idx[],int *nswap,int *n);

        void EliminateDup(std::vector <Vertex> &sortedVertexList);
        void AddDummyLF();
        void InitFreeKill();

        void Jump(std::vector<Vertex> & vertexList, int *iseed, int ival, int *itetra);
        void InsideTetraJW(std::vector<Vertex> & vertexList,int p, int a, int b, int c, int d, int iorient, int *is_in, int *redundant, int *ifail);
        void LocateByJumpWalk(std::vector<Vertex> & vertexList,int *iseed, int ival,int *tetraLoc,int *iredundant);
        void FindTetra(int itetra, int idx_c, int a, int b, int o, int *ifind, int *tetra_loc, int *idx_a, int *idx_b);
        void DefineFacet(int itetra, int jtetra, int idx_o, int itouch[], int idx[], int jtouch[], int jdx[]);
        void RegularConvex(std::vector<Vertex> & vertexList,int a, int b, int c, int p, int o, int itest_abcp, bool *regular, bool *convex, bool *test_abpo, bool *test_bcpo, bool *test_capo);

        void FlipJW1_4(int ipoint, int itetra);
        void FlipJW(std::vector<Vertex> & vertexList);
        void FlipJW2_3(int itetra, int jtetra, int p, int o, int itetra_touch[], int itetra_idx[], int jtetra_touch[], int jtetra_idx[], bool *test_abpo, bool *test_bcpo, bool *test_capo, int *ierr);
        void FlipJW3_2(int itetra, int jtetra, int ktetra, int vertices[], int itetra_touch[], int itetra_idx[], int jtetra_touch[], int jtetra_idx[], int ktetra_touch[], int ktetra_idx[], bool *test_bcpo, bool *test_acpo,  int *ierr);
        void FlipJW4_1(std::vector<Vertex>& vertexList,int itetra, int jtetra, int ktetra, int ltetra, int vertices[], int ishare, int idx, int jshare, int jdx, int kshare, int kdx, int lshare, int ldx, bool *test_acpo, int *ierr);

        void MissInfSign(int i,int j,int k,int *l,int *sign);
        void Peel(std::vector<Vertex> & vertexList);
        void MarkZero(int itetra, int ivertex);
        void RemoveInf();

        void CNInit();
        void CNClean();
        void MarkIterate(std::vector<Vertex> &vertexList);
        int MarkCH(std::vector<Vertex> &vertexList);
        void MarkNeighbours(int t,int n_iter);
        void MarkNeighboursInSet(int n_iter,int i_beg,int i_end);
        void MarkNotMarked();

    public:
        int redundantCount;
        std::vector <Tetrahedron> DeluanayTet;
        std::vector <Triangle> DeluanayTrigs;
        std::vector <Edge> DeluanayEdges;

        DeluanayComplex();
        ~DeluanayComplex();

        double ran2(int *idum);
        void ConstructDT(std::vector<Vertex> & vertexList);
        void DefineTriangles();
        void DefineEdges(std::vector<Vertex> & vertexList);
        void CalculateNormals(std::vector<Vertex> &vertexList);
        void CorrectNormals(std::vector<Vertex> &vertexList);
};

#endif // DELUANAYCOMPLEX_H
