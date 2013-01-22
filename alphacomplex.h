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

#ifndef ALPHACOMPLEX_H
#define ALPHACOMPLEX_H

#define ALF_MAXIMUM  ((float) 1E+37)
#define ALF_MINIMUM  ((float) 1E-37)
#define ALF_ZERO     ((float) 0.0)
#define ALF_INFINITY ((float) 1E+37)

/* f_type/r_type enumerators; f_type == dimension + 1 */
#define ALF_TETRA    4
#define ALF_TRIANGLE 3
#define ALF_EDGE     2
#define ALF_VERTEX   1    /* ^                                       */
#define ALF_BLANK    0    /* |- 3 bits f_type  and  2 bits |- r_type */
#define ALF_RHO      1    /*                               V         */
#define ALF_MU1      2
#define ALF_MU2      3

#include <vector>
#include <algorithm>
#include <cmath>
#include <alfmasternode.h>
#include <mheapnode.h>
#include <tlistnode.h>
#include <deluanaycomplex.h>
#include <size.h>
#include <disjointset.h>
#include <depth.h>
#include <deluanayflow.h>
#include <lightmaterial.h>
#include <GL/glew.h>
#include <skinsurface.h>
#include <rankmap.h>
#include <simplexmasterlistmap.h>
#include <QLineEdit>
#include <volume.h>
#include <powerdiagram.h>

class AlphaComplex
{
    private:
        AlfMasterNode *mlNode;
        std::vector<int> mlmarks;
        std::vector<int> mlrepeats;
        std::vector<int> betti0;
        Depth *depth;
        DeluanayFlow *dflow;
        SkinSurface *alphaSkin;
        DisJointSet *unionFind;
        std::vector<std::vector<int> > AllComponents;
        Volume *volume;
        double Epsilon;
        int LowRank;
        bool isFiltrationModified;
        Processor* processor;


        int MaxRank(int s, int r);
        int MinRank(int s, int r);

        int SoftInfinity(int last);
        int Add(std::vector <TlistNode> & sortedTListNode,int i, int last);
        void Put(mpz_t *p, mpz_t *q, int i, unsigned int ftype);
        void Push(int ix, uint fType, uint rType, int r);


        int MlEofSublist(int r);
        unsigned int MlRType();
        int MlPrev();
        int MlIsAttached();

        void SpectrumTetra(int i);
        void SpectrumTriangle(int i);
        void SpectrumEdge(int i);
        void SpectrumVertex(int i);

        void VertexMus();
        void EdgeMus();
        void TriangleMus();

        void CollectMaster();

        void MarkSimplices();
        void UpdateBoundaries();
        void PairSimplices();
        void UpdatePersistenceInfo();

        void CalculatePersistence();
        void WritePersistence();

        void fillCandidateTrigs(int lowRank,int highRank,std::vector <RankMap> &candidateTrigs);
        void fillCandidateTets(int lowRank,int highRank,std::vector <RankMap> &candidateTets);

    public:
        int MaximumRank;
        int MaximumPersistence;
        double eps;
        int sp_j;
        int currentRank;
        int entries;
        int hmax;
        int hn;
        int max;
        int ranks;
        int vals;
        unsigned int ftype;
        unsigned int rtype;

        std::vector <double> Spectrum;
        std::vector <int> MasterList;
        std::vector <int> auxi;
        std::vector <AlfMasterNode> MasterNode;
        std::vector <TlistNode> tListNode;
        std::vector <MHeapNode> mHeapNode;
        std::vector <Vertex> vertexList;
        std::vector <int> sortedTet;
        std::vector <int> sortedTrigs;

        DeluanayComplex *delcx;
        Size *size;

        AlphaComplex(Processor* process, std::vector <Vertex> &vertlist);
        AlphaComplex(std::vector <Vertex> &vertlist,double center[],double * scale);
        ~AlphaComplex();

        double AlphaSqrt(double v);
        double AlphaThresholdSqr(int r);
        int AlphaRank(double v);
        int AlphaRankSqr(double v);
        double AlphaThreshold(int r);
        double RealValue(mpz_t *a, mpz_t *b);

        void SpectrumOpen();
        void Tetra();
        void Trigs();
        void Edges();
        void Vertices();
        void SpectrumClose();

        void BuildSpectrum();
        int getMaxRank();
        int getMaxPersistence();

        AlfMasterNode * GetMLNode();

        int MlSublist (int r);
        unsigned int MlFType();
        int MlNext();
        int MlIsFirst();

        void BuildComplex(int Rank,std::vector <Vertex> &vertlist);
        void Adjust(double real_alpha,std::vector <Vertex> &vertlist);

        bool InComplex(std::vector<Vertex> &vertexList,unsigned int ftype, int rank,int i);

        void Render(int rank,bool al,bool wf,bool skin,bool ss,bool swf);
        void RenderUnModified(int rank,bool al,bool wf,bool skin,bool ss,bool swf);
        void RenderModified(int rank,bool al,bool wf,bool skin,bool ss,bool swf);
        void RenderComplement (bool wf);
        void PrintML();

        void FindProperties(std::vector<Vertex> &vertexList,QLineEdit *totVol,QLineEdit *totSurf,int Rank);

        void MarkRelevantSimplices(double epsilon,int lowRank,int highRank,std::vector<RankMap> &candidateTrigs,std::vector<RankMap> &candidateTets);
        void ModifyFiltration(std::vector <SimplexMasterListMap> &refinedCandidateTrigs,std::vector <SimplexMasterListMap> &refinedCandidateTets);
        void UndoModifyFiltration(std::vector <SimplexMasterListMap> &refinedCandidateTrigs,std::vector <SimplexMasterListMap> &refinedCandidateTets);
};

#endif // ALPHACOMPLEX_H
