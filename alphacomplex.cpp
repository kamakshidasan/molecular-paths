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

#include "alphacomplex.h"

#define Switch(i, j) { int tmp; tmp = i; i = j; j = tmp; }

typedef struct uf_rec_type
{
        int p;
} Uf_rec;

static Uf_rec* ufds = NULL;

FILE *fp1;

static void uf_create (int n)
{
    if (ufds == NULL)
    {
        ufds = (Uf_rec*)malloc(n*sizeof(Uf_rec));
    }
    bzero((char*)ufds, sizeof(Uf_rec)* n);
}

static int uf_find (int i)
{
    int j,t;
    j = i;
    while (ufds[j].p >= 0)
    {                                           /* find the root of this set                */
        j = ufds[j].p;
    }                                           /* j is the index of the root of set of i   */
    while (ufds[i].p >= 0)
    {                                           /* PATH COMPRESSION :  make nodes on path   */
        t = i;                                  /*   from i to root have root as parent     */
        i = ufds[i].p;
        ufds[t].p = j;
    }
    //fprintf(fp1,"find(%d) = %d\n",i,j);
    return (j);                                 /* j is the index of the root of set of i   */
}

static void uf_add (int i)
{
    ufds[i].p = -1;				/* element i defines a set by itself and has size 1 */
    //fprintf(fp1,"adding element %d\n",i);
}

static int uf_union (int i, int j)
{
    i = uf_find(i);				/* find the root of set containing i      */
    j = uf_find(j);				/* find the root of set containing j      */
    if (i == j)
    {
        return (0);                             /* no union effected                      */
    }
    else                                        /* if different sets                      */
    {
                                                /* WEIGHT BALANCING                       */
        if (ufds[i].p                           /* ufds[i].p = -1 *                       */
                                                /*            number of descendants of i  */
                      > ufds[j].p
           )
                Switch(i, j);                   /* i must be the larger set               */
        ufds[i].p = ufds[i].p + ufds[j].p;
                                                /* no. of descendants of i increases      */
        ufds[j].p = i;                          /* j now points to i                      */
        //fprintf(fp1,"union -> set1 = %d set2 = %d\n",i,j);
        return (1);                             /* union effected                         */
    }
}

static int myMax(int a, int b)
{
    if (b > a) return b;
    return a;
}

static int myMin(int a, int b)
{
    if (b < a) return b;
    return a;
}

/*!
    \fn AlphaComplex::AlphaComplex(std::vector <Vertex> &vertlist)
 */
AlphaComplex::AlphaComplex(std::vector <Vertex> &vertlist)
{
    double center[3] = {0.0,0.0,0.0};
    this->vertexList = vertlist;
    delcx = new DeluanayComplex();
    size = new Size();
    dflow = new DeluanayFlow();
    depth = new Depth();
    delcx->ConstructDT(vertexList);
    delcx->DefineTriangles();
    delcx->DefineEdges(vertexList);
    delcx->CalculateNormals(vertexList);
//    delcx->CorrectNormals(vertexList);
    MaximumRank = currentRank = vals = MaximumPersistence = 0;
    dflow->CalculateDF(delcx,vertexList);
    sortedTet.reserve(delcx->DeluanayTet.size()-delcx->redundantCount);
    alphaSkin = new SkinSurface(center,0.0,0);
    unionFind = new DisJointSet();
    volume = new Volume(vertexList,vertexList.size(),delcx->DeluanayTrigs.size(),delcx->DeluanayEdges.size());
    isFiltrationModified = false;
}
/*!
    \fn AlphaComplex::AlphaComplex(std::vector <Vertex> &vertlist)
 */
AlphaComplex::AlphaComplex(std::vector <Vertex> &vertlist,double center[],double *scale)
{
    //double center[3] = {0.0,0.0,0.0};
    this->vertexList = vertlist;
    delcx = new DeluanayComplex();
    size = new Size();
    dflow = new DeluanayFlow();
    depth = new Depth();
    delcx->ConstructDT(vertexList);
    delcx->DefineTriangles();
    delcx->DefineEdges(vertexList);
    delcx->CalculateNormals(vertexList);
    //delcx->CorrectNormals(vertexList);
    MaximumRank = currentRank = vals = MaximumPersistence = 0;
    dflow->CalculateDF(delcx,vertexList);
    sortedTet.reserve(delcx->DeluanayTet.size()-delcx->redundantCount);
    alphaSkin = new SkinSurface(center,*scale,0);
    unionFind = new DisJointSet();
    volume = new Volume(vertexList,vertexList.size(),delcx->DeluanayTrigs.size(),delcx->DeluanayEdges.size());
    isFiltrationModified = false;
}

/*!
    \fn AlphaComplex::~AlphaComplex()
 */
AlphaComplex::~AlphaComplex()
{
    delete delcx;
    delete size;
}

/*!
    \fn AlphaComplex::MaxRank(int s, int r)
 */
int AlphaComplex::MaxRank(int s, int r)
{
    return ((r > 0) ? myMax(s, r) : s);
}

/*!
    \fn AlphaComplex::MinRank(int s, int r)
 */
int AlphaComplex::MinRank(int s, int r)
{
    return ((r > 0) ? myMin(s, r) : s);
}

/*!
    \fn AlphaComplex::AlphaSqrt(double v)
 */
double AlphaComplex::AlphaSqrt(double v)
{
    bool negative = (v < 0);

    double abs_v = (negative) ? -v : v;

    if (abs_v >= ALF_INFINITY)
        return (v);
    else if (negative)
        return (-sqrt(abs_v));
    else
        return (sqrt(abs_v));
}

/*!
    \fn AlphaComplex::AlphaThreshold(int r)
 */
double AlphaComplex::AlphaThreshold(int r)
{
    return(AlphaSqrt(AlphaThresholdSqr(r)));
}

/*!
    \fn AlphaComplex::AlphaThresholdSqr(int r)
 */
double AlphaComplex::AlphaThresholdSqr(int r)
{
    r = myMax(0, r);
    r = myMin(r, MaximumRank);
    return (Spectrum[r]);
}

int AlphaComplex::AlphaRank(double v)
{
    double sign = 0.0;
    if (v < 0.0) sign = -1.0;
    else if( v > 0.0) sign = 1.0;
    return(AlphaRankSqr(sign * v * v));
}

int AlphaComplex::AlphaRankSqr(double v)
{
    int left = 1,right = MaximumRank,mid;
    if(v <= Spectrum[left])
        return left;
    else if(v >= Spectrum[right])
        return right;
    else
    {
        while(1)
        {
            if(right - left <= 1) break;
            else
            {
                mid = (left + right) / 2;
                if(v < Spectrum[mid])
                    right = mid;
                else
                    left = mid;
            }
        }
        return left;
    }
}

int AlphaComplex::getMaxRank()
{
    return MaximumRank;
}

int AlphaComplex::getMaxPersistence()
{
    return MaximumPersistence;
}

bool AlphaComplex::InComplex(std::vector<Vertex> &vertexList,unsigned int ftype, int rank,int i)
{
    //Assert(alpha_is_valid()) not required because of the way it is used
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

/*!
    \fn AlphaComplex::RealValue(mpz_t *a, mpz_t *b)
 */
double AlphaComplex::RealValue(mpz_t *a, mpz_t *b)
{
    float value;

    assert (mpz_sgn (*b) > 0);					//check no (hard) infinity

    if (mpz_sgn (*a) == 0)						//check for zero
    {
        value = ALF_ZERO;
    }
    else
    {
        double ratio = mpz_get_d (*a) / mpz_get_d (*b);		//normal ratio rho^2 = a/b

        assert (ratio > -ALF_MAXIMUM);				//no -infinity


        if (((ratio >=0) && (ratio < ALF_MINIMUM)) ||
                            ((ratio < 0) &&
                            (-ratio < ALF_MINIMUM)))		//(soft) zero
        {
            value = ALF_ZERO;
        }
        else if (ratio > ALF_MAXIMUM)				//(soft) infinity
        {
            value = ALF_INFINITY;
        }
        else
        {
            value =  ratio;
        }
    }
    return (value);
}

/*!
    \fn AlphaComplex::SoftInfinity(int last)
 */
int AlphaComplex::SoftInfinity(int last)
{
    int i = last;
    while (Spectrum[i] == ALF_INFINITY)
            i--;
    return (last - i);
}

/*!
    \fn AlphaComplex::Add(int i, int last)
 */
int AlphaComplex::Add(std::vector <TlistNode> & sortedTListNode,int i, int last)
{
    double newvalue;

    assert ((last ?	(ranks > 0):((ranks == 0) || (ranks == 1))));

    newvalue = RealValue(&sortedTListNode[i].a, &sortedTListNode[i].b);
    newvalue = newvalue / 1E16;

    if ((last == 0) || (size->AlfRatioCompare(&sortedTListNode[last].a,&sortedTListNode[last].b,&sortedTListNode[i].a,&sortedTListNode[i].b) < 0))
    {
        /*
        different rho-values (possibly same FP-values) ==> different rank
        */
        if((ranks > 0) && (Spectrum[ranks] > newvalue))
                newvalue = Spectrum[ranks];
        /*
        adjust new value if it's out of order due to floating-point error
        */
        assert (ranks < vals + 2);
        ranks++;
        Spectrum[ranks] = newvalue;
        last = i;
    }
    else
    {
        /*
        same rho-values ==> same rank :P of course
        */
        assert ((ranks > 0) && (Spectrum[ranks] == newvalue));
    }
    /*
    assign rank
    */
    sortedTListNode[i].r = ranks;
    switch (sortedTListNode[i].ftype)
    {
        case ALF_TETRA: delcx->DeluanayTet[sortedTListNode[i].ix].Rho = ranks;
        break;
        case ALF_TRIANGLE: delcx->DeluanayTrigs[sortedTListNode[i].ix].Rho = ranks;
        break;
        case ALF_EDGE: delcx->DeluanayEdges[sortedTListNode[i].ix].Rho = ranks;
        break;
        case ALF_VERTEX: vertexList[sortedTListNode[i].ix].Rho = ranks;
        break;
        default:
                break;
    }
    Push(sortedTListNode[i].ix, sortedTListNode[i].ftype, sortedTListNode[i].rtype, ranks);
    return (last);
}

/*!
    \fn AlphaComplex::Put(mpz_t *p, mpz_t *q, int i, unsigned int ftype)
 */
void AlphaComplex::Put(mpz_t *p, mpz_t *q, int i, unsigned int ftype)
{
    assert(max>vals);

    vals++;
    uint rtype = ALF_RHO;

    if (mpz_cmp_d(*p,0.0) == 0)
    {
        mpz_set_d(tListNode[vals].a,0.0);
        mpz_set_d(tListNode[vals].b,1.0);
    }
    else
    {
        mpz_set(tListNode[vals].a,*p);
        mpz_set(tListNode[vals].b,*q);
    }
    tListNode[vals].r = -1;
    tListNode[vals].ix = i;
    tListNode[vals].si = vals;
    tListNode[vals].rtype = rtype;
    tListNode[vals].ftype = ftype;
}

/*!
    \fn AlphaComplex::Push(int ix, uint fType, uint rType, int r)
 */
void AlphaComplex::Push(int ix, uint fType, uint rType, int r)
{
    if (r > 0)
    {
        assert(hn < hmax);
        hn++;
        MHeapNode mnode(ix, auxi[r], fType, rType);
        mHeapNode.push_back(mnode);
        auxi[r] = hn;
    }
}

/*!
    \fn AlphaComplex::MlSublist (int r)
 */
int AlphaComplex::MlSublist (int r)
{
    assert ((0 <= r) && (r <= MaximumRank));
    sp_j = MasterList[r];
    if ((sp_j != 0))
    {
        assert ((0 <= sp_j) && (sp_j <= entries));
        mlNode = &MasterNode[sp_j];
        ftype = mlNode->ftype;
        return (sp_j);
    }
    else
    {
        sp_j = (r >= ranks)? entries:0;
        mlNode = &MasterNode[sp_j];
        ftype = ALF_BLANK;
        return (0);
    }
}

/*!
    \fn AlphaComplex::MlEofSublist(int r)
 */
int AlphaComplex::MlEofSublist(int r)
{
    assert ((0 <= r) and (r <= MaximumRank));
    sp_j = MasterList[r];
    if (sp_j)
    {
        assert ((0 <= sp_j) && (sp_j <= entries));                  /* compute "next" j */
        {
            int j_prime;
            r++;
            while ((r <= MaximumRank) && (!MasterList[r]))
            {
                r++;
            }
            j_prime = ((r > MaximumRank)? entries: MasterList[r] - 1);
            assert (j_prime >= sp_j);
            sp_j = j_prime;
        }
        mlNode = &MasterNode[sp_j];
        ftype = mlNode->ftype;
        assert((0 <= sp_j) && (sp_j <= entries));
        return (sp_j);
    }
    else
    {
        sp_j = ((r > MaximumRank)? entries: 0);
        mlNode = &MasterNode[sp_j];
        ftype = ALF_BLANK;
        return (0);
    }
}

/*!
    \fn AlphaComplex::MlIsFirst()
 */
int AlphaComplex::MlIsFirst()
{
    return((MlRType() == ALF_RHO) || ((MlRType() == ALF_MU1) && (MlIsAttached())));
}

/*!
    \fn AlphaComplex::MlIsAttached()
 */
int AlphaComplex::MlIsAttached()
{
    switch (ftype)
    {
        case ALF_TRIANGLE:
                return (delcx->DeluanayTrigs[mlNode->ix].Rho == 0);
        case ALF_EDGE:
                return (delcx->DeluanayEdges[mlNode->ix].Rho == 0);
        case ALF_VERTEX:
                return (vertexList[mlNode->ix].Rho == 0);
        default:
                return (0);
    }
}

/*!
    \fn AlphaComplex::MlFType()
 */
unsigned int AlphaComplex::MlFType()
{
    return(ftype);
}

/*!
    \fn AlphaComplex::MlRType()
 */
unsigned int AlphaComplex::MlRType()
{
    if(ftype == ALF_BLANK)
    {
        return(ALF_BLANK);
    }
    else
    {
        return(mlNode->rtype);
    }
}

/*!
    \fn AlphaComplex::MlNext ()
 */
int AlphaComplex::MlNext ()
{
    if ((ftype == ALF_BLANK) || mlNode->last)
        return (0);
    else
    {
        sp_j++;
        assert ((0 <= sp_j) && (sp_j <= entries));
        mlNode = &MasterNode[sp_j];
        ftype = mlNode->ftype;
        return (sp_j);
    }
}

/*!
    \fn AlphaComplex::MlPrev()
 */
int AlphaComplex::MlPrev()
{
    AlfMasterNode* Prev = mlNode - 1;
    assert(mlNode);
    if((ftype == ALF_BLANK) || Prev->last)
        return(0);
    else
    {
        sp_j--;
        assert ((0 <= sp_j) && (sp_j <= entries));
        mlNode = &MasterNode[sp_j];
        ftype = mlNode->ftype;
        return(sp_j);
    }
}

/*!
    \fn AlphaComplex::SpectrumOpen()
 */
void AlphaComplex::SpectrumOpen()
{
    max = delcx->DeluanayTet.size() + delcx->DeluanayTrigs.size() + delcx->DeluanayEdges.size() + vertexList.size();
    hmax = 3 * max - 2 * delcx->DeluanayTet.size();
    mpz_t p,q;
    mpz_init(p);mpz_init(q);
    mpz_set_d(p,0.0);mpz_set_d(q,1.0);

    tListNode.reserve(max+1);

    for(int i=0;i<max+1;i++)
    {
        mpz_init(tListNode[i].a);
        mpz_init(tListNode[i].b);
        mpz_set_d(tListNode[i].a,0.0);
        mpz_set_d(tListNode[i].b,1.0);

        tListNode[i].r = -1;
        tListNode[i].ix = 0;
        tListNode[i].si = 0;
        tListNode[i].rtype = rtype;
        tListNode[i].ftype = ftype;
    }
}

/*!
    \fn AlphaComplex::Tetra()
 */
void AlphaComplex::Tetra()
{
    for (uint i = 0; i < delcx->DeluanayTet.size(); i++)
    {
        if (delcx->DeluanayTet[i].Status == 1)
        {
            SpectrumTetra(i);
        }
    }
}

/*!
    \fn AlphaComplex::SpectrumTetra(int i)
 */
void AlphaComplex::SpectrumTetra(int i)
{
    int a, b, c, d;
    mpz_t p;mpz_init(p);
    mpz_t q;mpz_init(q);
    uint ftype = ALF_TETRA;

    a = delcx->DeluanayTet[i].Corners[1];
    b = delcx->DeluanayTet[i].Corners[2];
    c = delcx->DeluanayTet[i].Corners[3];
    d = delcx->DeluanayTet[i].Corners[4];

    size->TetraSize(delcx,vertexList, a, b, c, d, &p, &q, i);
    Put(&p, &q, i, ftype);
    mpz_clear(p);mpz_clear(q);
}

/*!
    \fn AlphaComplex::Trigs()
 */
void AlphaComplex::Trigs()
{
    for (uint i = 1; i < delcx->DeluanayTrigs.size(); i++)
    {
        SpectrumTriangle(i);
    }
}

/*!
    \fn AlphaComplex::SpectrumTriangle(int i)
 */
void AlphaComplex::SpectrumTriangle(int i)
{
    int a, b, c;
    mpz_t p;
    mpz_t q;
    uint ftype = ALF_TRIANGLE;
    if (!delcx->DeluanayTrigs[i].IsAttached)
    {
        mpz_init(p);mpz_init(q);

        a = delcx->DeluanayTrigs[i].Corners[1];
        b = delcx->DeluanayTrigs[i].Corners[2];
        c = delcx->DeluanayTrigs[i].Corners[3];

        size->TrigSize(vertexList, a, b, c, &p, &q);
        Put(&p, &q, i, ftype);
        mpz_clear(p);mpz_clear(q);
    }
}

/*!
    \fn AlphaComplex::Edges()
 */
void AlphaComplex::Edges()
{
    for (uint i = 1; i < delcx->DeluanayEdges.size(); i++)
    {
        //SpectrumEdge(i,writer);
        SpectrumEdge(i);
    }
}

/*!
    \fn AlphaComplex::SpectrumEdge(int i)
 */
void AlphaComplex::SpectrumEdge(int i)
{
    int a, b;
    mpz_t p;
    mpz_t q;
    uint ftype = ALF_EDGE;
    if (!delcx->DeluanayEdges[i].IsAttached)
    {
        mpz_init(p);mpz_init(q);

        a = delcx->DeluanayEdges[i].Corners[1];
        b = delcx->DeluanayEdges[i].Corners[2];

        size->EdgeSize(vertexList, a, b, &p, &q);
        Put(&p, &q, i, ftype);
        mpz_clear(p);mpz_clear(q);
    }
}

/*!
    \fn AlphaComplex::Vertices
 */
void AlphaComplex::Vertices()
{
    Vertex *rvi;
    Vertex *rvj;
    int p = 0, q = 0;

    for (uint i = 1; i < delcx->DeluanayEdges.size(); i++)
    {
        rvi = &vertexList[delcx->DeluanayEdges[i].Corners[1]];
        p = delcx->DeluanayEdges[i].Corners[1];

        rvj = &vertexList[delcx->DeluanayEdges[i].Corners[2]];
        q = delcx->DeluanayEdges[i].Corners[2];

        if (rvi->Rho >= 0)
        {
            if (size->CheckVertex(vertexList,p, q) > 0)
            {
                rvi->Rho = -1;
            }
            else
                rvi->Rho = 1;
        }
        if (rvj->Rho >= 0)
        {
            if (size->CheckVertex(vertexList,q, p) > 0)
            {
                rvj->Rho = -1;
            }
            else
                rvj->Rho = 1;
        }
    }
    for (uint i = 1; i < vertexList.size(); i++)
    {
        rvi = &vertexList[i];
        switch (rvi->Rho)
        {
            case 0:
                    rvi->Rho = -1;
                    break;
            case -1:
                    rvi->Rho = 0;
                    break;
            case 1:
                    SpectrumVertex(i);
                    break;
            default:
                    assert(true);
                //no need for break;
        }
    }
}

/*!
    \fn AlphaComplex::SpectrumVertex(int i)
 */
void AlphaComplex::SpectrumVertex(int i)
{
    mpz_t p;
    mpz_t q;
    uint ftype = ALF_VERTEX;

    mpz_init(p);mpz_init(q);

    size->VertSize(vertexList, i, &p, &q);
    Put(&p, &q, i, ftype);

    mpz_clear(p);mpz_clear(q);
}

/*!
    \fn AlphaComplex::CollectMaster()
 */
void AlphaComplex::CollectMaster()
{
    int hi, low, i, j, r, swaps = 0,mlindex = 0;
    AlfMasterNode keyNode;
    /*
    the [0] element
    */
    entries = 0;
    MasterNode[0].ix = 0;
    MasterNode[0].ftype = ALF_BLANK;
    MasterNode[0].rtype = ALF_BLANK;
    MasterNode[0].last = true;
    MasterNode[0].isValid = true;

    for (r = 1; r <= ranks; r++)
    {
        hi = auxi[r];
        MasterList[r] = ((hi != 0) ? (entries + 1) : 0);
        while (hi  != 0)
        {
            assert((hn > 0) && (entries < hmax));
            hn--;
            entries++;
            mlindex++;
            MasterNode[mlindex].ftype = mHeapNode[hi].ftype;
            MasterNode[mlindex].rtype = mHeapNode[hi].rtype;
            MasterNode[mlindex].ix = mHeapNode[hi].ix;
            MasterNode[mlindex].isValid = true;
            hi = mHeapNode[hi].nexti;
            MasterNode[mlindex].last = (hi == 0);
        }
    }
    assert(hn == 0);
    /*
    sort sublists w/ insertion sort, because they're already nearly sorted
    */
    hi = entries;
    for (r = ranks; r > 1; r--)
    {
        low = MasterList[r];
        if (low != 0)
        {
            MasterNode[hi].last = false;			//temporarily!
            for (j = low + 1; j <= hi; j++)			//insertion sort on master_####[low..high]
            {
                keyNode = MasterNode[j];
                i = j - 1;
                while ((i >= low) && (MasterNode[i].ftype > keyNode.ftype))
                {
                    MasterNode[i + 1] = MasterNode[i];
                    swaps++;
                    i--;
                }
                if (i != j - 1)
                {
                    MasterNode[i + 1] = keyNode;
                }
            }
            MasterNode[hi].last = true;
            hi = low - 1;
        }
    }
}

/*!
    \fn AlphaComplex::VertexMus()
 */
void AlphaComplex::VertexMus()
{
    int r, redundant = 0;
    uint i;
    Edge *re;
    Vertex *rv;
    uint ftype,rtype;

    for( i = 1; i < delcx->DeluanayEdges.size(); i++ )  		//(un)attached edge  ==> 2 bounding vertices
    {
        re = &delcx->DeluanayEdges[i];
        r = re->Rho;
        if (r != 0)                             			//unattached edge
        {
            rv = &vertexList[delcx->DeluanayEdges[i].Corners[1]];
            rv->Mu1 = MinRank (r, rv->Mu1);
            rv->Mu2 = MaxRank (re->Mu2, rv->Mu2);

            rv = &vertexList[delcx->DeluanayEdges[i].Corners[2]];
            rv->Mu1 = MinRank(r, rv->Mu1);
            rv->Mu2 = MaxRank(re->Mu2, rv->Mu2);
        }
        else                                    			//attached edge
        {
            rv = &vertexList[delcx->DeluanayEdges[i].Corners[1]];
            rv->Mu1 = MinRank (re->Mu1, rv->Mu1);
            rv->Mu2 = MaxRank (re->Mu2, rv->Mu2);

            rv = &vertexList[delcx->DeluanayEdges[i].Corners[2]];
            rv->Mu1 = MinRank (re->Mu1, rv->Mu1);
            rv->Mu2 = MaxRank (re->Mu2, rv->Mu2);
        }
    }

    for( i = 1; i < delcx->DeluanayTrigs.size(); i++ )
    {
        if (delcx->DeluanayTrigs[i].Hull == 1)				//Remove mu2 values for vertices on CH!
        {
            vertexList[delcx->DeluanayTrigs[i].Corners[1]].Mu2 = 0;
            vertexList[delcx->DeluanayTrigs[i].Corners[2]].Mu2 = 0;
            vertexList[delcx->DeluanayTrigs[i].Corners[3]].Mu2 = 0;

            assert    ((vertexList[delcx->DeluanayTrigs[i].Corners[1]].Mu1 != 0)
                    && (vertexList[delcx->DeluanayTrigs[i].Corners[2]].Mu1 != 0)
                    && (vertexList[delcx->DeluanayTrigs[i].Corners[3]].Mu1 != 0));

        }
    }

    for (i = 1; i < vertexList.size(); i++)
    {
        rv = &vertexList[i];
        assert(
                (-1 <= rv->Rho) &&
                (rv->Rho <= rv->Mu1) &&
                (
                    (rv->Mu1 <= rv->Mu2) ||
                    (rv->Mu2 == 0)
                ) &&
                (rv->Mu2 <= ranks)
              );

        if ((rv->Mu1 == 0) && (rv->Mu2 == 0))
        {
            assert(rv->Rho == -1);
            redundant++;
        }
        ftype = ALF_VERTEX;
        rtype = ALF_MU1;
        Push(i, ftype, rtype, rv->Mu1);
        rtype = ALF_MU2;
        Push(i, ftype, rtype, rv->Mu2);
    }
}

/*!
    \fn AlphaComplex::EdgeMus()
 */
void AlphaComplex::EdgeMus()
{
    int r;
    Triangle *rf;
    Edge *re;
    uint rtype,ftype,i;

    for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
    {
        rf = &delcx->DeluanayTrigs[i];
        r = rf->Rho;

        re = &delcx->DeluanayEdges[delcx->DeluanayTrigs[i].TrigLink[1]];
        if (r != 0)								//unattached triangle
        {
            re->Mu1 = MinRank(r, re->Mu1);
            re->Mu2 = MaxRank(rf->Mu2, re->Mu2);
        }
        else 									//attached triangle
        {
            re->Mu1 = MinRank(rf->Mu1, re->Mu1);
            re->Mu2 = MaxRank(rf->Mu2, re->Mu2);
        }

        re = &delcx->DeluanayEdges[delcx->DeluanayTrigs[i].TrigLink[2]];
        if (r != 0)								//unattached triangle
        {
            re->Mu1 = MinRank(r, re->Mu1);
            re->Mu2 = MaxRank(rf->Mu2, re->Mu2);
        }
        else 									//attached triangle
        {
            re->Mu1 = MinRank(rf->Mu1, re->Mu1);
            re->Mu2 = MaxRank(rf->Mu2, re->Mu2);
        }

        re = &delcx->DeluanayEdges[delcx->DeluanayTrigs[i].TrigLink[3]];
        if (r != 0)								//unattached triangle
        {
            re->Mu1 = MinRank(r, re->Mu1);
            re->Mu2 = MaxRank(rf->Mu2, re->Mu2);
        }
        else 									//attached triangle
        {
            re->Mu1 = MinRank(rf->Mu1, re->Mu1);
            re->Mu2 = MaxRank(rf->Mu2, re->Mu2);
        }
    }

    for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
    {
        if (delcx->DeluanayTrigs[i].Hull == 1)					//Remove mu2 values for edges on CH!
        {
            delcx->DeluanayEdges[delcx->DeluanayTrigs[i].TrigLink[1]].Mu2 = 0;
            delcx->DeluanayEdges[delcx->DeluanayTrigs[i].TrigLink[2]].Mu2 = 0;
            delcx->DeluanayEdges[delcx->DeluanayTrigs[i].TrigLink[3]].Mu2 = 0;
        }
    }

    for (i = 1; i < delcx->DeluanayEdges.size(); i++)
    {
        re = &delcx->DeluanayEdges[i];
        assert(
                (0 <= re->Rho) &&
                (re->Rho <= re->Mu1) &&
                (
                    (re->Mu1 <= re->Mu2) ||
                    (re->Mu2 == 0)
                ) &&
                (re->Mu2 <= ranks)
              );
        ftype = ALF_EDGE;
        rtype = ALF_MU1;
        Push(i, ftype, rtype, re->Mu1);
        rtype = ALF_MU2;
        Push(i, ftype, rtype, re->Mu2);
    }
}

/*!
    \fn AlphaComplex::TriangleMus()
 */
void AlphaComplex::TriangleMus()
{
    uint i,rtype,ftype;
    int r;
    Triangle *rf;

    for (i = 1; i < delcx->DeluanayTet.size(); i++)		//tetrahedron: ==> 4 bounding triangles
    {
        if (delcx->DeluanayTet[i].Status == 0)
        {
            assert (delcx->DeluanayTet[i].Rho == 0);
            continue;
        }
        else
        {
            assert (delcx->DeluanayTet[i].Rho > 0);
            r = delcx->DeluanayTet[i].Rho;

            rf = &delcx->DeluanayTrigs[delcx->DeluanayTet[i].TetLink[1]];
            rf->Mu1 = MinRank(r, rf->Mu1);
            rf->Mu2 = MaxRank(r, rf->Mu2);

            rf = &delcx->DeluanayTrigs[delcx->DeluanayTet[i].TetLink[2]];
            rf->Mu1 = MinRank(r, rf->Mu1);
            rf->Mu2 = MaxRank(r, rf->Mu2);

            rf = &delcx->DeluanayTrigs[delcx->DeluanayTet[i].TetLink[3]];
            rf->Mu1 = MinRank(r, rf->Mu1);
            rf->Mu2 = MaxRank(r, rf->Mu2);

            rf = &delcx->DeluanayTrigs[delcx->DeluanayTet[i].TetLink[4]];
            rf->Mu1 = MinRank(r, rf->Mu1);
            rf->Mu2 = MaxRank(r, rf->Mu2);

        }
    }
    for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
    {
        rf = &delcx->DeluanayTrigs[i];

        if (delcx->DeluanayTrigs[i].Hull == 1)      //remove mu2 value for CH triangles!
        {
            rf->Mu2 = 0;
        }
        assert (
                (
                    (0 <= rf->Rho) &&
                    (rf->Rho <= rf->Mu1) &&
                    (
                        (rf->Mu1 <= rf->Mu2) ||
                        (rf->Mu2 == 0)
                    ) &&
                    (rf->Mu2 <= ranks)
                )
               );
        ftype = ALF_TRIANGLE;
        rtype = ALF_MU1;
        Push(i, ftype, rtype, rf->Mu1);
        rtype = ALF_MU2;
        Push(i, ftype, rtype, rf->Mu2);
    }
}

/*!
    \fn AlphaComplex::SpectrumClose()
 */
void AlphaComplex::SpectrumClose()
{
    int i, j, k;
    Spectrum.reserve(vals + 2 + 1);
    MasterList.reserve(vals + 2 + 1);
    auxi.reserve(vals + 2 + 1);

    for(i = 0;i< vals+2+1;i++)
    {
        auxi[i] = 0;
    }

    FILE *fp = fopen("unsorted.txt","w");
    for (i = 1; i < vals + 1; i++)
    {
        //fprintf(fp,"%d:%d %d %lf\n",i,tListNode[i].ix,tListNode[i].ftype,RealValue(&tListNode[i].a,&tListNode[i].b)/1E16);
    }
    fclose(fp);

    std::vector <TlistNode> sortedTListNode(tListNode.begin(),tListNode.begin()+vals+1);
    tListNode.clear();

    std::sort(sortedTListNode.begin()+1,sortedTListNode.end());

    k = 0;
    fp = fopen("sortedtet.txt","w");
    for (i = 1; i < vals + 1; i++)
    {
        if (sortedTListNode[i].ftype == ALF_TETRA)
        {
            k++;
            sortedTet[k] = sortedTListNode[i].ix;
            fprintf(fp,"%d\n",sortedTet[k]);
        }
    }
    fclose(fp);

    fp = fopen("sorted.txt","w");
    for (i = 1; i < vals + 1; i++)
    {
        //fprintf(fp,"%d:%d %d %lf\n",sortedTListNode[i].si,sortedTListNode[i].ix,sortedTListNode[i].ftype,RealValue(&sortedTListNode[i].a,&sortedTListNode[i].b)/1E16);
    }
    fclose(fp);

    hn = entries = ranks = 0;

    MHeapNode mnode(-1, -1, ALF_BLANK, ALF_BLANK);
    mHeapNode.push_back(mnode);

    Spectrum[0] = -ALF_INFINITY;
    MasterList[0] = 0;

    i = 1;
    j = Add(sortedTListNode,i, 0);
    while (i < vals)
    {
        i++;
        j = Add(sortedTListNode,i, j);
    }

    i = SoftInfinity(ranks);

    ranks++;
    Spectrum[ranks] = ALF_INFINITY;
    MasterList[ranks] = 0;
    assert (ranks < (vals + 2));

    TriangleMus();
    EdgeMus();
    VertexMus();

    MasterNode.reserve(hn + 1);

    CollectMaster();

    MaximumRank = ranks - 1;

    auxi.clear();
    mHeapNode.clear();

    //realloc is it correct?
    /*Spectrum.TrimExcess();
    MasterList.TrimExcess();
    MasterNode.TrimExcess();*/

    sortedTListNode.clear();
    vals = hmax = max = 0;
}

/*!
    \fn AlphaComplex::BuildSpectrum()
 */
void AlphaComplex::BuildSpectrum()
{
    SpectrumOpen();
    Tetra();
    Trigs();
    Edges();
    Vertices();
    SpectrumClose();
    PrintML();
    depth->CalculateDepth(delcx,sortedTet);
    CalculatePersistence();
}

/*!
    \fn AlphaComplex::BuildComplex(int Rank,std::vector <Vertex> &vertlist)
 */
void AlphaComplex::BuildComplex(int Rank,std::vector <Vertex> &vertlist)
{
    int j;
    uint i;
    for (i = 1; i < vertexList.size(); i++)
    {
        vertexList[i].AlphaStatus = -1;
        vertlist[i].AlphaStatus = -1;
    }
    for (i = 1; i < delcx->DeluanayEdges.size(); i++)
    {
        delcx->DeluanayEdges[i].AlphaStatus = 0;
    }
    for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
    {
        delcx->DeluanayTrigs[i].AlphaStatus = 0;
        delcx->DeluanayTrigs[i].isValid = false;
        delcx->DeluanayTrigs[i].trigCoef = 0;
    }
    for (i = 1; i < delcx->DeluanayTet.size(); i++)
    {
        delcx->DeluanayTet[i].AlphaStatus = 0;
        delcx->DeluanayTet[i].isValid = false;
    }

    assert (Rank < MaximumRank);

    for (j = 1; j <= Rank; j++)
    {
        if (MlSublist(j) != 0)
        {
            do
            {
                if(MlIsFirst())
                {
                    if((Rank == MaximumRank-1) || (mlNode->isValid))
                    {
                        /*
                        first check if the rank is max rank -1 then we need not worry about modifications and
                        hence can add all simplices. otherwise add it to the complex only if it is valid i.e. even
                        after modification it remains in the spectrum for the current rank
                        */
                        switch (MlFType())
                        {
                            case ALF_VERTEX:    {
                                                    vertexList[mlNode->ix].AlphaStatus = 0;
                                                    vertlist[mlNode->ix].AlphaStatus = 0;
                                                }
                            break;
                            case ALF_EDGE:      {
                                                    delcx->DeluanayEdges[mlNode->ix].AlphaStatus = 1;
                                                    delcx->DeluanayEdges[mlNode->ix].RenderFlag = 1;
                                                }

                            break;
                            case ALF_TRIANGLE:  {
                                                    delcx->DeluanayTrigs[mlNode->ix].AlphaStatus = 1;
                                                    delcx->DeluanayTrigs[mlNode->ix].isValid = true;
                                                    int e = delcx->DeluanayTrigs[mlNode->ix].TrigLink[1];
                                                    delcx->DeluanayEdges[e].RenderFlag = 2;
                                                    e = delcx->DeluanayTrigs[mlNode->ix].TrigLink[2];
                                                    delcx->DeluanayEdges[e].RenderFlag = 2;
                                                    e = delcx->DeluanayTrigs[mlNode->ix].TrigLink[3];
                                                    delcx->DeluanayEdges[e].RenderFlag = 2;
                                                }
                            break;
                        case ALF_TETRA:         {
                                                    delcx->DeluanayTet[mlNode->ix].AlphaStatus = 1;
                                                    delcx->DeluanayTet[mlNode->ix].isValid = true;
                                                }
                            break;
                            default:            break;
                        }
                    }
                    /*else
                    {
                        switch (MlFType())
                        {
                           case ALF_TRIANGLE:  {
                                                    for(uint _i = 1; _i <= 3; _i++)
                                                    {
                                                        int vert = delcx->DeluanayTrigs[mlNode->ix].Corners[_i];
                                                        vertexList[vert].valid = false;
                                                    }
                                               }
                            break;
                            case ALF_TETRA:     {
                                                    for(uint _i = 1; _i <= 4; _i++)
                                                    {
                                                        int vert = delcx->DeluanayTet[mlNode->ix].Corners[_i];
                                                        vertexList[vert].valid = false;
                                                    }
                                                }
                            break;
                            default:            break;
                        }
                    }*/
                    /*AlfMasterNode * tempNode = mlNode->next;
                    while(tempNode)
                    {
                        switch (tempNode->ftype)
                        {
                            case ALF_VERTEX:    {
                                                    vertexList[tempNode->ix].AlphaStatus = 0;
                                                    vertlist[tempNode->ix].AlphaStatus = 0;
                                                }
                            break;
                            case ALF_EDGE:      {
                                                    delcx->DeluanayEdges[tempNode->ix].AlphaStatus = 1;
                                                    delcx->DeluanayEdges[tempNode->ix].RenderFlag = 1;
                                                }

                            break;
                            case ALF_TRIANGLE:  {
                                                    delcx->DeluanayTrigs[tempNode->ix].AlphaStatus = 1;
                                                    int e = delcx->DeluanayTrigs[tempNode->ix].TrigLink[1];
                                                    delcx->DeluanayEdges[e].RenderFlag = 2;
                                                    e = delcx->DeluanayTrigs[tempNode->ix].TrigLink[2];
                                                    delcx->DeluanayEdges[e].RenderFlag = 2;
                                                    e = delcx->DeluanayTrigs[tempNode->ix].TrigLink[3];
                                                    delcx->DeluanayEdges[e].RenderFlag = 2;

                                                    for(uint _i = 1; _i <= 3; _i++)
                                                    {
                                                        int vert = delcx->DeluanayTrigs[mlNode->ix].Corners[_i];
                                                        vertexList[vert].valid = true;
                                                    }
                                                }
                            break;
                            case ALF_TETRA:     {
                                                    delcx->DeluanayTet[tempNode->ix].AlphaStatus = 1;

                                                    for(uint _i = 1; _i <= 4; _i++)
                                                    {
                                                        int vert = delcx->DeluanayTet[mlNode->ix].Corners[_i];
                                                        vertexList[vert].valid = true;
                                                    }
                                                }
                            break;
                            default:            break;
                        }
                    }*/
                }
            }
            while ((MlNext() != 0));
        }
    }

    for (i = 1; i < delcx->DeluanayTrigs.size(); i++)
    {
        if (delcx->DeluanayTrigs[i].AlphaStatus == 1)
        {
            delcx->DeluanayTrigs[i].trigCoef = 2;
        }
        else
        {
            delcx->DeluanayTrigs[i].trigCoef = 0;
        }
    }

    for (i = 1; i < delcx->DeluanayTet.size(); i++)
    {
        if (delcx->DeluanayTet[i].AlphaStatus == 1)
        {
            for (j = 1; j <= 4; j++)
            {
                delcx->DeluanayTrigs[delcx->DeluanayTet[i].TetLink[j]].trigCoef = delcx->DeluanayTrigs[delcx->DeluanayTet[i].TetLink[j]].trigCoef - 1;
            }
        }
        else
        {
            if (delcx->DeluanayTet[i].Status != 0)
            {
                delcx->DeluanayTet[i].AlphaStatus = -1;
                delcx->DeluanayTet[i].isValid = false;
            }
        }
    }
}

/*!
    \fn AlphaComplex::MarkSimplices()
 */
void AlphaComplex::MarkSimplices()
{
    int mlix = 0;
    int r;
    mlmarks.reserve(entries+1);
    mlrepeats.reserve(entries+1);
    betti0.reserve(MaximumRank+1);
    betti0[0] = betti0[1] = 0;

    int repeats = 0;

    uf_create(vertexList.size()+1);
    for(r = 1;r<=MaximumRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                mlix++;
                if(MlFType() == ALF_VERTEX)
                {
                    vertexList[mlNode->ix].Repeats = repeats;
                    if(MlIsFirst())
                    {
                        mlmarks[mlix] = 1;
                        uf_add(mlNode->ix);
                        betti0[r]++;
                    }
                    else
                    {
                        mlmarks[mlix] = -1;
                    }
                }
                else if(MlFType() == ALF_EDGE)
                {
                    delcx->DeluanayEdges[mlNode->ix].Repeats = repeats;
                    int i = delcx->DeluanayEdges[mlNode->ix].Corners[1];
                    int j = delcx->DeluanayEdges[mlNode->ix].Corners[2];
                    if(MlIsFirst())
                    {
                        int did_unite = uf_union(i,j);
                        if(did_unite)//union affected
                        {
                            mlmarks[mlix] = 0;
                            betti0[r]--;
                        }
                        else
                        {
                            mlmarks[mlix] = 1;
                        }
                    }
                    else
                    {
                        mlmarks[mlix] = -1;
                    }
                }
            }while(MlNext());
        }
        if(r<MaximumRank)
        {
            betti0[r+1] = betti0[r];
        }
    }
    free(ufds);
    ufds = NULL;

    uf_create(delcx->DeluanayTet.size()+1);
    mlix = entries+1;
    uf_add(0);
    for(r=MaximumRank;r>0;r--)
    {
        if(MlEofSublist(r))
        {
            do
            {
                mlix--;
                if(MlFType() == ALF_TETRA)
                {
                    mlmarks[mlix] = 0;
                    uf_add(mlNode->ix);
                }
                else if(MlFType() == ALF_TRIANGLE)
                {
                    if(MlIsFirst())
                    {
                        int i,j;
                        if(delcx->DeluanayTrigs[mlNode->ix].nLink == 2)
                        {
                            i = delcx->DeluanayTrigs[mlNode->ix].ReverseLink1;
                            j = delcx->DeluanayTrigs[mlNode->ix].ReverseLink2;
                        }
                        else if(delcx->DeluanayTrigs[mlNode->ix].nLink == 1)
                        {
                            i = delcx->DeluanayTrigs[mlNode->ix].ReverseLink1;
                            j = 0;
                        }
                        else
                        {
                            assert(true);
                        }
                        int did_unite = uf_union(i,j);
                        if(did_unite)//union affected
                        {
                            mlmarks[mlix] = 1;
                        }
                        else
                        {
                            mlmarks[mlix] = 0;
                        }
                    }
                    else
                    {
                        mlmarks[mlix] = -1;
                    }
                }
            }while(MlPrev());
        }
    }

    mlix = 0;
    for(r = 1;r<=MaximumRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                mlix++;

                if(MlFType() == ALF_VERTEX)
                {
                    mlrepeats[mlix] = repeats;
                    if(mlmarks[mlix] == -1)
                    {
                        repeats++;
                    }
                }
                else if(MlFType() == ALF_EDGE)
                {
                    mlrepeats[mlix] = repeats;
                    if(mlmarks[mlix] == -1)
                    {
                        repeats++;
                    }
                }
                else if(MlFType() == ALF_TRIANGLE)
                {
                    mlrepeats[mlix] = repeats;
                    if(mlmarks[mlix] == -1)
                    {
                        repeats++;
                    }
                }
                else if(MlFType() == ALF_TETRA)
                {
                    mlrepeats[mlix] = repeats;
                    if(mlmarks[mlix] == -1)
                    {
                        repeats++;
                    }
                }

            }while(MlNext());
        }
    }
}

/*!
    \fn AlphaComplex::UpdateBoundaries()
 */
void AlphaComplex::UpdateBoundaries()
{
    int mlix = 0;
    for(int r=1;r<MaximumRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                mlix++;
                if(MlFType() == ALF_TRIANGLE)
                {
                    if(MlIsFirst())
                    {
                        int tet;
                        if(delcx->DeluanayTrigs[mlNode->ix].nLink == 1)
                        {
                            tet = delcx->DeluanayTrigs[mlNode->ix].ReverseLink1;
                            delcx->DeluanayTet[tet].BoundarySimplices.push_back(mlix);
                        }
                        else if(delcx->DeluanayTrigs[mlNode->ix].nLink == 2)
                        {
                            tet = delcx->DeluanayTrigs[mlNode->ix].ReverseLink1;
                            delcx->DeluanayTet[tet].BoundarySimplices.push_back(mlix);
                            tet = delcx->DeluanayTrigs[mlNode->ix].ReverseLink2;
                            delcx->DeluanayTet[tet].BoundarySimplices.push_back(mlix);
                        }
                    }
                }
            }while(MlNext());
        }
    }
}

/*!
    \fn AlphaComplex::PairSimplices()
 */
void AlphaComplex::PairSimplices()
{
    int mlix = 0;
    std::vector<int> ptemp1;
    std::vector<int> ptemp2;
    int i,j,k,kk,jj,maxIndex;
    AlfMasterNode * tempNode;

    for(int r=1;r<MaximumRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                mlix++;
                if((mlmarks[mlix]!=0) || (MlFType()!=ALF_TETRA))
                {
                    continue;
                }
                else
                {
                    jj = mlix;
                    int size1 = delcx->DeluanayTet[mlNode->ix].BoundarySimplices.size();
                    for(i=0;i<size1;i++)
                    {
                        ptemp1.push_back(delcx->DeluanayTet[mlNode->ix].BoundarySimplices[i]);
                    }
                    maxIndex = size1 - 1;

                    jj = ptemp1[maxIndex];

                    while(jj>=0)
                    {
                        tempNode = &MasterNode[jj];

                        if(delcx->DeluanayTrigs[tempNode->ix].Entry == -1)
                        {
                            delcx->DeluanayTrigs[tempNode->ix].Entry = mlix;
                            delcx->DeluanayTet[mlNode->ix].posIndex = jj;
                            size1 = ptemp1.size();
                            for(i=0;i<size1;i++)
                            {
                                    delcx->DeluanayTrigs[tempNode->ix].BoundarySimplices.push_back(ptemp1[i]);
                            }
                            break;
                        }
                        else
                        {
                            size1 = ptemp1.size();
                            int size2 = delcx->DeluanayTrigs[tempNode->ix].BoundarySimplices.size();
                            i = j = kk = 0;
                            for(k=0;k<size1+size2;k++)
                            {
                                if(i == size1)
                                {
                                    for(;j<size2;j++)
                                    {
                                        ptemp2.push_back(delcx->DeluanayTrigs[tempNode->ix].BoundarySimplices[j]);
                                        kk++;
                                    }
                                    break;
                                }
                                if(j == size2)
                                {
                                    for(;i<size1;i++)
                                    {
                                        ptemp2.push_back(ptemp1[i]);
                                        kk++;
                                    }
                                    break;
                                }
                                if (ptemp1[i] < delcx->DeluanayTrigs[tempNode->ix].BoundarySimplices[j])
                                {
                                    ptemp2.push_back(ptemp1[i]);
                                    kk++;
                                    i++;
                                }
                                else if (ptemp1[i] > delcx->DeluanayTrigs[tempNode->ix].BoundarySimplices[j])
                                {
                                    ptemp2.push_back(delcx->DeluanayTrigs[tempNode->ix].BoundarySimplices[j]);
                                    j++;
                                    kk++;
                                }
                                else
                                {
                                    i++;
                                    j++;
                                }
                            }
                            maxIndex = ptemp2.size() - 1;
                            ptemp1.clear();
                            for(uint q =0;q< ptemp2.size();q++) //ptemp1 = ptemp2;
                            {
                                ptemp1.push_back(ptemp2[q]);
                            }
                            ptemp2.clear();
                            jj = ptemp1[maxIndex];
                        }
                    }
                }
                ptemp1.clear();
                ptemp2.clear();

            }while(MlNext());
        }
    }
}

/*!
    \fn AlphaComplex::UpdatePersistenceInfo()
 */
void AlphaComplex::UpdatePersistenceInfo()
{
    int r,j;int mlix = 0;
    AlfMasterNode * tempNode;

    mlix = 0;
    for(r=1;r<MaximumRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                mlix++;
                if((mlmarks[mlix] == 1) && (MlFType() == ALF_TRIANGLE))
                {
                    if(delcx->DeluanayTrigs[mlNode->ix].Entry != -1)
                    {
                        j = delcx->DeluanayTrigs[mlNode->ix].Entry;
                        tempNode = &MasterNode[j];
                        int mainRepeats = mlrepeats[j] - mlrepeats[mlix];
                        delcx->DeluanayTrigs[mlNode->ix].MainRepeats = delcx->DeluanayTet[tempNode->ix].MainRepeats = mainRepeats;
                        delcx->DeluanayTrigs[mlNode->ix].Persistence = delcx->DeluanayTet[tempNode->ix].Persistence = j - mlix - mainRepeats;

                        if(MaximumPersistence < delcx->DeluanayTet[tempNode->ix].Persistence)
                        {
                            MaximumPersistence = delcx->DeluanayTet[tempNode->ix].Persistence;
                        }
                    }
                }
            }while(MlNext());
        }
    }
}

/*!
    \fn AlphaComplex::WritePersistence()
 */
void AlphaComplex::WritePersistence()
{
    int r;int mlix = 0;
    FILE *fp = fopen("persist.txt","w");
    for(r=1;r<MaximumRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                mlix++;
                if(MlFType() == ALF_TRIANGLE)
                {
                    if(mlmarks[mlix] == 1)
                    {
                        if(delcx->DeluanayTrigs[mlNode->ix].Persistence == -1)
                        {
                            fprintf(fp,"Rank = %d Dimension = %d, Positive = %d, Negative = Infinity, Repeats = Invalid, Persistence = Infinity\n",r,2,mlNode->ix);
                        }
                        else
                        {
                            fprintf(fp,"Rank = %d Dimension = %d, Positive = %d, Negative = %d, Repeats = %d, Persistence = %d\n",r,2,mlix,delcx->DeluanayTrigs[mlNode->ix].Entry,delcx->DeluanayTrigs[mlNode->ix].MainRepeats,delcx->DeluanayTrigs[mlNode->ix].Persistence);
                        }
                    }
                }
            }while(MlNext());
        }
    }
    mlmarks.clear();
    mlrepeats.clear();
    fclose(fp);
}

/*!
    \fn AlphaComplex::CalculatePersistence()
 */
void AlphaComplex::CalculatePersistence()
{
    MarkSimplices();
    UpdateBoundaries();
    PairSimplices();
    UpdatePersistenceInfo();
    WritePersistence();
}

/*!
    \fn AlphaComplex::Adjust(double real_alpha,std::vector <Vertex> &vertlist)
 */
void AlphaComplex::Adjust(double real_alpha,std::vector <Vertex> &vertlist)
{
    double temp,abs_v;
    double temp_alpha = -Epsilon * Epsilon;
    int negative;
    for(uint i=0;i<vertexList.size();i++)
    {
        if(vertexList[i].valid)
        {
            temp = vertexList[i].BackRadius * vertexList[i].BackRadius + real_alpha;
            vertexList[i].Weight = vertlist[i].Weight = -temp;
            negative = (temp<0);

            abs_v = ((negative)?(-temp):(temp));
            if(abs_v >= ALF_INFINITY)
            {
                vertlist[i].Radius = vertexList[i].Radius = temp;
            }
            else if(negative)
            {
                vertlist[i].Radius = vertexList[i].Radius = -sqrt(abs_v);
            }
            else
            {
                vertlist[i].Radius = vertexList[i].Radius = sqrt(abs_v);
            }

            for(int j=1;j<4;j++)
            {
                temp = vertexList[i].Coordinates[j];
                vertlist[i].Weight = vertexList[i].Weight = vertexList[i].Weight + temp * temp;
            }
            vertlist[i].Coef = vertexList[i].Coef = 1;
        }
        else
        {
            temp = vertexList[i].BackRadius * vertexList[i].BackRadius + temp_alpha;
            vertexList[i].Weight = vertlist[i].Weight = -temp;
            negative = (temp<0);

            abs_v = ((negative)?(-temp):(temp));
            if(abs_v >= ALF_INFINITY)
            {
                vertlist[i].Radius = vertexList[i].Radius = temp;
            }
            else if(negative)
            {
                vertlist[i].Radius = vertexList[i].Radius = -sqrt(abs_v);
            }
            else
            {
                vertlist[i].Radius = vertexList[i].Radius = sqrt(abs_v);
            }

            for(int j=1;j<4;j++)
            {
                temp = vertexList[i].Coordinates[j];
                vertlist[i].Weight = vertexList[i].Weight = vertexList[i].Weight + temp * temp;
            }
            vertlist[i].Coef = vertexList[i].Coef = 1;
        }
    }
}
/*!
    \fn AlphaComplex::RenderUnModified(int rank,bool al,bool wf,bool skin,bool ss,bool swf)
 */
void AlphaComplex::RenderUnModified (int rank, bool al, bool wf, bool skin, bool ss, bool swf)
{
    uint i;
    int j;
    int mlix = 0;
    int r;
    double a[3],b[3],c[3],v[3];
    std::vector<Vertex> skinList;
    FILE *fp;

    if(!skin)
    {
        if(al)
        {
            if(!wf)
            {
                glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
            }
            else
            {
                glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
            }

            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[3]);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);

            glBegin(GL_TRIANGLES);
            for(i=1;i<delcx->DeluanayTrigs.size();i++)
            {
                if(delcx->DeluanayTrigs[i].AlphaStatus == 1)
                {
                    for(j=0;j<3;j++)
                    {
                        a[j] = vertexList[delcx->DeluanayTrigs[i].Corners[1]].NormCoordinates[j+1];
                        b[j] = vertexList[delcx->DeluanayTrigs[i].Corners[2]].NormCoordinates[j+1];
                        c[j] = vertexList[delcx->DeluanayTrigs[i].Corners[3]].NormCoordinates[j+1];
                    }
                    v[0] = delcx->DeluanayTrigs[i].Normal->X;
                    v[1] = delcx->DeluanayTrigs[i].Normal->Y;
                    v[2] = delcx->DeluanayTrigs[i].Normal->Z;

                    glNormal3dv(v);
                    glVertex3dv(a);
                    glVertex3dv(b);
                    glVertex3dv(c);

                }
            }
            glEnd();

            glBegin(GL_LINES);
            for(i=1;i<delcx->DeluanayEdges.size();i++)
            {
                if((delcx->DeluanayEdges[i].AlphaStatus == 1) && (delcx->DeluanayEdges[i].RenderFlag == 1))
                {
                    for(j=0;j<3;j++)
                    {
                        a[j] = vertexList[delcx->DeluanayEdges[i].Corners[1]].NormCoordinates[j+1];
                        b[j] = vertexList[delcx->DeluanayEdges[i].Corners[2]].NormCoordinates[j+1];
                    }

                    glVertex3dv(a);
                    glVertex3dv(b);

                }
            }
            glEnd();

            glBegin(GL_POINTS);
            for(i=1;i<vertexList.size();i++)
            {
                if(vertexList[i].AlphaStatus == 0)
                {
                    for(j=0;j<3;j++)
                    {
                        a[j] = vertexList[i].NormCoordinates[j+1];
                    }

                    glVertex3dv(a);

                }
            }
            glEnd();
        }
    }
    else
    {
        if(rank!=currentRank)
        {
            system("rm -f alpha_level0.off");
            system("rm -f alpha");
            skinList.clear();
            printf("Calculating SkinSurface for the Molecule..........................\n");

            if(betti0[rank] == 1)
            {
                for(i = 1;i<vertexList.size();i++)
                {
                    if(vertexList[i].AlphaStatus == 0)
                    {
                        skinList.push_back(vertexList[i]);
                    }
                }
                fp = fopen("alpha","w");
                fprintf(fp,"%d\n",skinList.size());
                fprintf(fp,"#junk\n");
                for(uint ll = 0 ; ll < skinList.size() ; ll++)
                {
                        fprintf(fp,"%d %f %f %f %f\n",ll+1,skinList[ll].Coordinates[1],skinList[ll].Coordinates[2],skinList[ll].Coordinates[3],skinList[ll].Radius);
                }
                fclose(fp);
                skinList.clear();
                system("./smesh alpha -s alpha.off -t alpha.tet");
                alphaSkin->Read("alpha_lev0.off",0);
                alphaSkin->Process();
            }
            else if(betti0[rank]>1)
            {
                unionFind->Clear();
                for(r = 1;r<=rank;r++)
                {
                    if(MlSublist(r))
                    {
                        do
                        {
                            mlix++;
                            if(MlFType() == ALF_VERTEX)
                            {
                                if(MlIsFirst())
                                {
                                    vertexList[mlNode->ix].ufKey = unionFind->ElementCount();
                                    unionFind->Add(mlNode->ix);
                                }
                            }
                            else if(MlFType() == ALF_EDGE)
                            {
                                int i = delcx->DeluanayEdges[mlNode->ix].Corners[1];
                                int j = delcx->DeluanayEdges[mlNode->ix].Corners[2];
                                int s_i = vertexList[i].ufKey;
                                int s_j = vertexList[j].ufKey;
                                if(MlIsFirst())
                                {
                                    unionFind->Union(s_i,s_j);
                                }

                            }
                        }while(MlNext());
                    }
                }
                AllComponents = unionFind->Consolidate();
                for(uint _i = 0;_i<AllComponents.size();_i++)
                {
                    for(uint _j = 0;_j<AllComponents[_i].size();_j++)
                    {
                        skinList.push_back(vertexList[AllComponents[_i][_j]]);
                    }
                    fp = fopen("alpha","w");
                    fprintf(fp,"%d\n",skinList.size());
                    fprintf(fp,"#junk\n");
                    for(uint ll = 0 ; ll < skinList.size() ; ll++)
                    {
                            fprintf(fp,"%d %f %f %f %f\n",ll+1,skinList[ll].Coordinates[1],skinList[ll].Coordinates[2],skinList[ll].Coordinates[3],skinList[ll].Radius);
                    }
                    fclose(fp);
                    skinList.clear();
                    system("./smesh alpha -s alpha.off -t alpha.tet");
                    alphaSkin->Read("alpha_lev0.off",_i);
                }
                alphaSkin->Process();
            }
            currentRank = rank;
        }
        alphaSkin->Draw(ss,swf);
    }
}

/*!
    \fn AlphaComplex::Render(int rank,bool al,bool wf,bool skin,bool ss,bool swf)
 */
void AlphaComplex::RenderModified (int rank, bool al, bool wf, bool skin, bool ss, bool swf)
{
    uint i;
    int j;
    int mlix = 0;
    int r;
    double a[3],b[3],c[3],v[3];
    std::vector<Vertex> skinList;
    FILE *fp;

    if(!skin)
    {
        if(al)
        {
            if(!wf)
            {
                glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
            }
            else
            {
                glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
            }

            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[3]);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);

            for(i=1;i<delcx->DeluanayTrigs.size();i++)
            {
                if((delcx->DeluanayTrigs[i].AlphaStatus == 1) && (delcx->DeluanayTrigs[i].isValid == true))
                {
                    for(j=0;j<3;j++)
                    {
                        a[j] = vertexList[delcx->DeluanayTrigs[i].Corners[1]].NormCoordinates[j+1];
                        b[j] = vertexList[delcx->DeluanayTrigs[i].Corners[2]].NormCoordinates[j+1];
                        c[j] = vertexList[delcx->DeluanayTrigs[i].Corners[3]].NormCoordinates[j+1];
                    }
                    v[0] = delcx->DeluanayTrigs[i].Normal->X;
                    v[1] = delcx->DeluanayTrigs[i].Normal->Y;
                    v[2] = delcx->DeluanayTrigs[i].Normal->Z;
                    glBegin(GL_TRIANGLES);
                    glNormal3dv(v);
                    glVertex3dv(a);
                    glVertex3dv(b);
                    glVertex3dv(c);
                    glEnd();
                }
                else
                {
                    if(delcx->DeluanayTrigs[i].isValid == false)
                    {
                        printf("reached false\n");
                    }
                }
            }

            for(i=1;i<delcx->DeluanayEdges.size();i++)
            {
                if((delcx->DeluanayEdges[i].AlphaStatus == 1) && (delcx->DeluanayEdges[i].RenderFlag == 1))
                {
                    for(j=0;j<3;j++)
                    {
                        a[j] = vertexList[delcx->DeluanayEdges[i].Corners[1]].NormCoordinates[j+1];
                        b[j] = vertexList[delcx->DeluanayEdges[i].Corners[2]].NormCoordinates[j+1];
                    }
                    glBegin(GL_LINES);
                    glVertex3dv(a);
                    glVertex3dv(b);
                    glEnd();
                }
            }

            for(i=1;i<vertexList.size();i++)
            {
                if(vertexList[i].AlphaStatus == 0)
                {
                    for(j=0;j<3;j++)
                    {
                        a[j] = vertexList[i].NormCoordinates[j+1];
                    }
                    glBegin(GL_POINTS);

                    glVertex3dv(a);
                    glEnd();
                }
            }
        }
    }
    else
    {
        if(rank!=currentRank)
        {
            system("rm -f alpha_level0.off");
            system("rm -f alpha");
            skinList.clear();
            printf("Calculating SkinSurface for the Molecule..........................\n");

            if(betti0[rank] == 1)
            {
                for(i = 1;i<vertexList.size();i++)
                {
                    if(vertexList[i].AlphaStatus == 0)
                    {
                        skinList.push_back(vertexList[i]);

                    }
                }
                fp = fopen("alpha","w");
                fprintf(fp,"%d\n",skinList.size());
                fprintf(fp,"#junk\n");
                for(uint ll = 0 ; ll < skinList.size() ; ll++)
                {
                        fprintf(fp,"%d %f %f %f %f\n",ll+1,skinList[ll].Coordinates[1],skinList[ll].Coordinates[2],skinList[ll].Coordinates[3],skinList[ll].Radius);
                }
                fclose(fp);
                skinList.clear();
                system("./smesh alpha -s alpha.off -t alpha.tet");
                alphaSkin->Read("alpha_lev0.off",0);
                alphaSkin->Process();
            }
            else if(betti0[rank]>1)
            {
                unionFind->Clear();
                for(r = 1;r<=rank;r++)
                {
                    if(MlSublist(r))
                    {
                        do
                        {
                            mlix++;
                            if(MlFType() == ALF_VERTEX)
                            {
                                if(MlIsFirst())
                                {
                                    vertexList[mlNode->ix].ufKey = unionFind->ElementCount();
                                    unionFind->Add(mlNode->ix);
                                }
                            }
                            else if(MlFType() == ALF_EDGE)
                            {
                                int i = delcx->DeluanayEdges[mlNode->ix].Corners[1];
                                int j = delcx->DeluanayEdges[mlNode->ix].Corners[2];
                                int s_i = vertexList[i].ufKey;
                                int s_j = vertexList[j].ufKey;
                                if(MlIsFirst())
                                {
                                    unionFind->Union(s_i,s_j);
                                }

                            }
                        }while(MlNext());
                    }
                }
                AllComponents = unionFind->Consolidate();
                for(uint _i = 0;_i<AllComponents.size();_i++)
                {
                    for(uint _j = 0;_j<AllComponents[_i].size();_j++)
                    {
                        skinList.push_back(vertexList[AllComponents[_i][_j]]);
                    }
                    fp = fopen("alpha","w");
                    fprintf(fp,"%d\n",skinList.size());
                    fprintf(fp,"#junk\n");
                    for(uint ll = 0 ; ll < skinList.size() ; ll++)
                    {
                            fprintf(fp,"%d %f %f %f %f\n",ll+1,skinList[ll].Coordinates[1],skinList[ll].Coordinates[2],skinList[ll].Coordinates[3],skinList[ll].Radius);
                    }
                    fclose(fp);
                    skinList.clear();
                    system("./smesh alpha -s alpha.off -t alpha.tet");
                    alphaSkin->Read("alpha_lev0.off",_i);
                }
                alphaSkin->Process();
            }
            currentRank = rank;
        }
        alphaSkin->Draw(ss,swf);
    }
}

/*!
    \fn AlphaComplex::Render(int rank,bool al,bool wf,bool skin,bool ss,bool swf)
 */
void AlphaComplex::Render(int rank,bool al,bool wf,bool skin,bool ss,bool swf)
{
    if(isFiltrationModified)
    {
        RenderModified (rank,al,wf,skin,ss,swf);
    }
    else
    {
        RenderUnModified (rank,al,wf,skin,ss,swf);
    }
}

/*!
    \fn AlphaComplex::PrintML()
 */
void AlphaComplex::PrintML()
{
    int r,mlix;
    mlix = 0;
    FILE *fp = fopen("masterlist.txt","w");
    for(r = 1;r<=MaximumRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                mlix++;
                if(MlFType() == ALF_VERTEX)
                {
                    if(MlIsFirst())
                    {
                        fprintf(fp,"masterlist index = %d ",mlix);
                        fprintf(fp,"vertexlist index = %d ",mlNode->ix);
                        fprintf(fp,"Rank = %d ",r);
                        fprintf(fp,"Dimension = %d ",MlFType() -1);
                        fprintf(fp,"Corners = %d\n",mlNode->ix);
                    }
                }
                else if(MlFType() == ALF_EDGE)
                {
                    if(MlIsFirst())
                    {
                        fprintf(fp,"masterlist index = %d ",mlix);
                        fprintf(fp,"edgelist index = %d ",mlNode->ix);
                        fprintf(fp,"Rank = %d ",r);
                        fprintf(fp,"Dimension = %d ",MlFType() -1);
                        int i = delcx->DeluanayEdges[mlNode->ix].Corners[1];
                        int j = delcx->DeluanayEdges[mlNode->ix].Corners[2];
                        fprintf(fp,"Corners = %d,%d\n",i,j);
                    }
                }
                else if(MlFType() == ALF_TRIANGLE)
                {
                    if(MlIsFirst())
                    {
                        fprintf(fp,"masterlist index = %d ",mlix);
                        fprintf(fp,"triglist index = %d ",mlNode->ix);
                        fprintf(fp,"Rank = %d ",r);
                        fprintf(fp,"Dimension = %d ",MlFType() -1);
                        fprintf(fp,"Corners = %d,%d,%d\n ",delcx->DeluanayTrigs[mlNode->ix].Corners[1],delcx->DeluanayTrigs[mlNode->ix].Corners[2],delcx->DeluanayTrigs[mlNode->ix].Corners[3]);
                    }
                }
                else if(MlFType() == ALF_TETRA)
                {
                    if(MlIsFirst())
                    {
                        fprintf(fp,"masterlist index = %d ",mlix);
                        fprintf(fp,"tetlist index = %d ",mlNode->ix);
                        fprintf(fp,"Rank = %d ",r);
                        fprintf(fp,"Dimension = %d ",MlFType() -1);
                        fprintf(fp,"Corners = %d,%d,%d,%d\n ",delcx->DeluanayTet[mlNode->ix].Corners[1],delcx->DeluanayTet[mlNode->ix].Corners[2],delcx->DeluanayTet[mlNode->ix].Corners[3],delcx->DeluanayTet[mlNode->ix].Corners[4]);
                    }
                }
            }while(MlNext());
        }
    }
    fclose(fp);
}

void AlphaComplex::fillCandidateTrigs(int lowRank, int highRank,std::vector <RankMap> &candidateTrigs)
{
    int r;
    RankMap rankmap(-1,-1,-1);
    for(r=lowRank;r<=highRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                if(MlFType() == ALF_TRIANGLE)
                {
                    if(MlIsFirst())
                    {
                        rankmap.rank = r;
                        rankmap.simplex = mlNode->ix;
                        rankmap.mlIndex = MasterList[r] - 1;
                        candidateTrigs.push_back(rankmap);
                    }
                }
            }while(MlNext());
        }
    }
}

void AlphaComplex::fillCandidateTets(int lowRank, int highRank,std::vector <RankMap> &candidateTets)
{
    int r;
    RankMap rankmap(-1,-1,-1);
    for(r=lowRank;r<=highRank;r++)
    {
        if(MlSublist(r))
        {
            do
            {
                if(MlFType() == ALF_TETRA)
                {
                    if(MlIsFirst())
                    {
                        rankmap.rank = r;
                        rankmap.simplex = mlNode->ix;
                        rankmap.mlIndex = MasterList[r] - 1;
                        candidateTets.push_back(rankmap);
                    }
                }
            }while(MlNext());
        }
    }
}

void AlphaComplex::MarkRelevantSimplices(double epsilon,int lowRank, int highRank,std::vector<RankMap> &candidateTrigs,std::vector<RankMap> &candidateTets)
{
    Epsilon = epsilon;
    LowRank = lowRank;
    fillCandidateTrigs(lowRank,highRank,candidateTrigs);
    fillCandidateTets(lowRank,highRank,candidateTets);
}

void AlphaComplex::ModifyFiltration(std::vector<SimplexMasterListMap> &refinedCandidateTrigs, std::vector<SimplexMasterListMap> &refinedCandidateTets)
{
    isFiltrationModified = true;
    FILE *fp = fopen("modified.txt","w");
    for(uint i =0;i<refinedCandidateTrigs.size();i++)
    {
        //find the triangle in the filtration
        int mlIndex = refinedCandidateTrigs[i].mlIndex;
        int simplex = refinedCandidateTrigs[i].simplex;

        //mark it as invalid
        MasterNode[mlIndex].isValid = false;

        //also mark it invalid for rendering
        //delcx->DeluanayTrigs[simplex].isValid = false;
        delcx->DeluanayTrigs[simplex].AlphaStatus = 0;
        fprintf(fp,"simplex = %d validity = %d\n",simplex,delcx->DeluanayTrigs[simplex].isValid);

        //important: mark the corners of this triangle so as to update its radius and weight later
        for(uint _i = 1; _i <= 3; _i++)
        {
            int vert = delcx->DeluanayTrigs[simplex].Corners[_i];
            vertexList[vert].valid = false;
        }

        //now find the tets of which the above triangle is a common face
        int tet1 = delcx->DeluanayTrigs[simplex].ReverseLink1;
        int tet2 = delcx->DeluanayTrigs[simplex].ReverseLink2;

        int tetrank1 = delcx->DeluanayTet[tet1].Rho;
        int tetrank2 = delcx->DeluanayTet[tet2].Rho;

        //mark them also as invalid so as to not violate the filtration
        int nextIndex = MasterList[tetrank1];
        MasterNode[nextIndex].isValid = false;
        //delcx->DeluanayTet[tet1].isValid = false;
        delcx->DeluanayTet[tet1].AlphaStatus = -1;
        fprintf(fp,"tet1 = %d validity = %d\n",tet1,delcx->DeluanayTet[tet1].isValid);

        nextIndex = MasterList[tetrank2];
        MasterNode[nextIndex].isValid = false;
        //delcx->DeluanayTet[tet2].isValid = false;
        delcx->DeluanayTet[tet2].AlphaStatus = -1;
        fprintf(fp,"tet2 = %d validity = %d\n",tet2,delcx->DeluanayTet[tet2].isValid);

        //important:mark the corners of this tet so as to update its radius and weight after this loop
        for(uint _i = 1; _i <= 4; _i++)
        {
            int vert = delcx->DeluanayTet[tet1].Corners[_i];
            vertexList[vert].valid = false;
            vert = delcx->DeluanayTet[tet2].Corners[_i];
            vertexList[vert].valid = false;
        }
    }

    /*for(uint i=0;i<refinedCandidateTets.size();i++)
    {
        int mlIndex = refinedCandidateTets[i].mlIndex;
        int simplex = refinedCandidateTets[1].simplex;
        //mark the tets as invalid
        MasterNode[mlIndex].isValid = false;

        //important:mark the corners of this tet so as to update its radius and weight after this loop
        for(uint _i = 1; _i <= 4; _i++)
        {
            int vert = delcx->DeluanayTet[simplex].Corners[_i];
            vertexList[vert].valid = false;
        }
    }*/
    fclose(fp);
}

void AlphaComplex::UndoModifyFiltration(std::vector<SimplexMasterListMap> &refinedCandidateTrigs, std::vector<SimplexMasterListMap> &refinedCandidateTets)
{
    //reverse the process done in AlphaComplex::ModifyFiltration ();
    for(uint i =0;i<refinedCandidateTrigs.size();i++)
    {
        int mlIndex = refinedCandidateTrigs[i].mlIndex;
        int simplex = refinedCandidateTrigs[i].simplex;

        MasterNode[mlIndex].isValid = true;

        //delcx->DeluanayTrigs[simplex].isValid = true;
        delcx->DeluanayTrigs[simplex].AlphaStatus = 1;
        for(uint _i = 1; _i <= 3; _i++)
        {
            int vert = delcx->DeluanayTrigs[simplex].Corners[_i];
            vertexList[vert].valid = true;
        }

        int tet1 = delcx->DeluanayTrigs[simplex].ReverseLink1;
        int tet2 = delcx->DeluanayTrigs[simplex].ReverseLink2;

        int tetrank1 = delcx->DeluanayTet[tet1].Rho;
        int tetrank2 = delcx->DeluanayTet[tet2].Rho;

        int nextIndex = MasterList[tetrank1];
        MasterNode[nextIndex].isValid = true;
        //delcx->DeluanayTet[tet1].isValid = true;
        delcx->DeluanayTet[tet1].AlphaStatus = 1;

        nextIndex = MasterList[tetrank2];
        MasterNode[nextIndex].isValid = true;
        //delcx->DeluanayTet[tet2].isValid = true;
        delcx->DeluanayTet[tet2].AlphaStatus = 1;

        //important:mark the corners of this tet so as to update its radius and weight after this loop
        for(uint _i = 1; _i <= 4; _i++)
        {
            int vert = delcx->DeluanayTet[tet1].Corners[_i];
            vertexList[vert].valid = true;
            vert = delcx->DeluanayTet[tet2].Corners[_i];
            vertexList[vert].valid = true;
        }
    }

    /*for(uint i=0;i<refinedCandidateTets.size();i++)
    {
        int mlIndex = refinedCandidateTets[i].mlIndex;
        int simplex = refinedCandidateTets[1].simplex;

        MasterNode[mlIndex].isValid = true;

        for(uint _i = 1; _i <= 4; _i++)
        {
            int vert = delcx->DeluanayTet[simplex].Corners[_i];
            vertexList[vert].valid = true;
        }
    }*/
    isFiltrationModified = false;
}

/*AlfMasterNode * AlphaComplex::GetMLNode()
{
    return mlNode;
}*/

void AlphaComplex::FindProperties (std::vector<Vertex> &vertexList,QLineEdit *totVol,QLineEdit *totSurf,int Rank)
{
    int idx,r,mlix = 0;

    double TotSurfaceArea = 0.0,TotVolume = 0.0;

     FILE *fp1 = fopen("totvolume.txt","w");

    for(r = 1;r<=Rank;r++)
    {
        if(MlSublist (r))
        {
            do
            {
                mlix++;
                if(MlFType () == ALF_VERTEX)
                {
                    if(MlIsFirst ())
                    {
                        idx = mlNode->ix;
                        volume->VertexProperties (vertexList,idx,fp1);
                    }
                }
                else if(MlFType () == ALF_EDGE)
                {
                    if(MlIsFirst ())
                    {
                        idx = mlNode->ix;

                        fprintf(fp1,"edge = %d\n",idx);
                        int i = delcx->DeluanayEdges[idx].Corners[1];
                        int j = delcx->DeluanayEdges[idx].Corners[2];
                        fprintf(fp1,"i = %d | j = %d\n",i,j);
                        volume->EdgeProperties (vertexList,i,j,fp1);
                    }
                }
                else if(MlFType () == ALF_TRIANGLE)
                {
                    if(MlIsFirst ())
                    {
                        idx = mlNode->ix;
                        fprintf(fp1,"trig = %d\n",idx);
                        int i = delcx->DeluanayTrigs[idx].Corners[1];
                        int j = delcx->DeluanayTrigs[idx].Corners[2];
                        int k = delcx->DeluanayTrigs[idx].Corners[3];
                        fprintf(fp1,"i = %d | j = %d | k =%d\n",i,j,k);
                        volume->TriangleProperties (vertexList,i,j,k,fp1);
                    }
                }
                else if(MlFType () == ALF_TETRA)
                {
                    if(MlIsFirst ())
                    {
                        idx = mlNode->ix;
                        fprintf(fp1,"tet = %d\n",idx);
                        int i = delcx->DeluanayTet[idx].Corners[1];
                        int j = delcx->DeluanayTet[idx].Corners[2];
                        int k = delcx->DeluanayTet[idx].Corners[3];
                        int l = delcx->DeluanayTet[idx].Corners[4];
                        fprintf(fp1,"i = %d | j = %d | k =%d | l = %d\n",i,j,k,l);
                        volume->TetProperties (vertexList,i,j,k,l,fp1);
                    }
                }
            }while(MlNext ());
        }
    }

    volume->FindTotal (vertexList,&TotVolume,&TotSurfaceArea);

    fclose(fp1);

    QString insertionString;
    QString &ran = insertionString.setNum(TotVolume,'f',5);
    totVol->setText(ran);

    ran = insertionString.setNum(TotSurfaceArea,'f',5);
    totSurf->setText(ran);
}

