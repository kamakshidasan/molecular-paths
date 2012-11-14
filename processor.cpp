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
#include "processor.h"

Processor::Processor()
{
        fr = new FileReader();
        alcx = 0;
        mouth = new Mouth();
        isRenderable = false;
        powerDiagram = NULL;
        proteinRenderer = NULL;
}

Processor::~Processor()
{
}

void Processor::read(const char *filename,double centre[],double *size, bool constantRadius, bool incrementRadius)
{
        double min[3], max[3];
        fr->ReadVertices(filename,vertexList,centre,&scale, min, max, constantRadius, incrementRadius);
        *size = scale;
        //alcx = new AlphaComplex(vertexList,centre,&scale);
        alcx = new AlphaComplex(vertexList);
        alcx->BuildSpectrum();
        powerDiagram = new PowerDiagram(alcx->delcx, vertexList, min, max);
        proteinRenderer = new ProteinRenderer(vertexList);
//        powerDiagram->printGraph();
        pocket = new Pocket(vertexList,centre,scale,alcx->delcx->DeluanayTrigs.size (),alcx->delcx->DeluanayEdges.size ());
        volume = new Volume(vertexList,vertexList.size(),alcx->delcx->DeluanayTrigs.size(),alcx->delcx->DeluanayEdges.size());
        isRenderable = true;
}

/*!
    \fn Processor::getMaxRank()
 */
int Processor::getMaxRank()
{
    return alcx->getMaxRank();
}

/*!
    \fn Processor::getRankForAlpha(double alphavalue)
 */
int Processor::getRankForAlpha(double alphavalue)
{
    return alcx->AlphaRank(alphavalue);
}

/*!
    \fn Processor::getMaxPersistence()
 */
int  Processor::getMaxPersistence()
{
    return alcx->getMaxPersistence();
}

/*!
    \fn Processor::getAlphaValue(int rank)
 */
double Processor::getAlphaValue(int rank)
{
    return alcx->AlphaThreshold(rank);
}

/*!
    \fn Processor::CalculateEverythingFor(int rank)
 */
void Processor::CalculateEverythingFor(int rank, int persistence)
{
        double real_alpha;
        if(alcx->MaximumRank!=0)
        {
                alcx->BuildComplex(rank,vertexList);
                real_alpha = alcx->AlphaThresholdSqr(rank);
                if(rank!=0)
                {
                        alcx->Adjust(real_alpha,vertexList);
                }
                if(rank != 0 && rank != alcx->MaximumRank)
                {
                        pocket->FindPockets(alcx->delcx,alcx->sortedTet);
                        pocket->MeasurePockets(alcx->delcx,vertexList);
                        //pocket->PocketProperties (alcx->delcx,vertexList,rank);
                        std::vector<std::vector<int> > pockets = pocket->GetPockets();
                        mouth->FindMouths(pockets,alcx->delcx);
                        mouth->MeasureMouths(alcx->delcx,vertexList);
                        volume->MeasureVolume(alcx->delcx,vertexList,totVol,totSurf,0);
                        //volume->MeasureVolume(alcx,vertexList,totVol,totSurf,0);
                        MapPocketsToMouths();
                        FillTables();
                        pspin->setRange(1,pockets.size());
                }
        }
}

/*!
    \fn Processor::CalculateRelevant(int rank)
 */
void Processor::CalculateRelevant(int rank, int persistence)
{
        double real_alpha;
        if(alcx->MaximumRank!=0)
        {
                alcx->BuildComplex(rank,vertexList);
                real_alpha = alcx->AlphaThresholdSqr(rank);
                if(rank!=0)
                {
                        alcx->Adjust(real_alpha,vertexList);
                }
                if(rank != 0 && rank != alcx->MaximumRank)
                {
                        pocket->FindPockets(alcx->delcx,alcx->sortedTet);
                        std::vector<std::vector<int> > pockets = pocket->GetPockets();
                        mouth->FindMouths(pockets,alcx->delcx);
                        MapPocketsToMouths();
                        pspin->setRange(1,pockets.size());
                }
        }
}

/*!
    \fn Processor::CalculateVolumes(int rank)
 */
void Processor::CalculateVolumes(int rank, int persistence)
{
        //double real_alpha;
        if(alcx->MaximumRank!=0)
        {
                if(rank != 0 && rank != alcx->MaximumRank)
                {
                        pocket->MeasurePockets(alcx->delcx,vertexList);
                        //pocket->PocketProperties (alcx->delcx,vertexList,rank);
                        mouth->MeasureMouths(alcx->delcx,vertexList);
                        volume->MeasureVolume(alcx->delcx,vertexList,totVol,totSurf,0);
                        //volume->MeasureVolume (alcx,vertexList,totVol,totSurf,0,rank);
                        //volume->FindVolume (alcx,vertexList,totVol,totSurf,0,rank);
                        //alcx->FindProperties (vertexList,totVol,totSurf,rank);
                        FillTables();
                }
        }

}

/*!
    \fn Processor::MapPocketsToMouths()
 */
void Processor::MapPocketsToMouths()
{
        PocketNMouths.clear();
        PocketNMouths.reserve(pocket->npockets);

        for(int i=0;i<pocket->npockets;i++)
        {
                PocketNMouths[i] = 0;
        }
        for(int i =0;i<mouth->nmouths;i++)
        {
                int j = mouth->mouthPocket(i,alcx->delcx);
                PocketNMouths[j]++;
        }
        mouth->CleanUp(alcx->delcx);
}

/*!
    \fn Processor::FillTables()
 */
void Processor::FillTables()
{
    pocket->FillTable(table1,PocketNMouths,pocVol,pocSurf);
    mouth->FillTable(table2,mouSurf);
}

bool Processor::IsRenderable()
{
    return isRenderable;
}

GLuint cHullId =-1;
GLuint cHullNormId = -1;

/*!
    \fn Processor::Render(int persistence,bool al,bool pc,bool vo,bool sk,bool mo,bool allindflag,int pnum,int rank,bool wf,bool fs)
 */
void Processor::Render(int persistence,bool alphaShape,bool allPockets,bool onlyPockets,bool onlyVoids,bool skinSurface,bool mouths,bool allindflag,int pnum,int rank,bool wireFrame,
                       bool alphaSkinSurface,bool alphaSkinWireFrame,bool smoothShading,bool skinWireFrame,bool pocketWireFrame, bool powerDiag,
                       bool cHull, bool cHullWF, bool cHullNorm, bool showPath, bool showSpaceFill)
{
        glEnable(GL_NORMALIZE);
        if(alcx)
        {
                alcx->Render(rank,alphaShape,wireFrame,alphaSkinSurface,smoothShading,alphaSkinWireFrame);

                if (onlyVoids || onlyPockets || allPockets)
                {
                        if(allindflag)
                        {
                                pocket->Render(alcx->delcx,vertexList,PocketNMouths,allPockets,onlyPockets,onlyVoids,skinSurface,persistence,rank,smoothShading,skinWireFrame,mouths,pocketWireFrame);
                        }
                        else
                        {
                                pocket->RenderSingle(alcx->delcx,vertexList,PocketNMouths,allPockets,onlyPockets,onlyVoids,skinSurface,persistence,rank,pnum,smoothShading,skinWireFrame,mouths,pocketWireFrame);
                        }
                }
                if (mouths)
                {
                        if (!onlyVoids)
                        {
                                if(allindflag)
                                {
                                        mouth->Render(alcx->delcx,vertexList,persistence);
                                }
                                else
                                {
                                        mouth->RenderMouthOfPocket(alcx->delcx,vertexList,persistence,pnum);
                                }
                        }
                }

            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[3]);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);

            if(proteinRenderer->init() && showSpaceFill){
                proteinRenderer->render();
            }
            powerDiagram->render(powerDiag, showPath);

            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[3]);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[3]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);

            glEnable(GL_COLOR_MATERIAL);
            if(cHull){
                glLineWidth(1);
                glColor3d(0.2, 0.6, 0.2);
                if(cHullWF){
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                }else{
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                }

                if (glIsList(cHullId) == GL_FALSE){
                    cHullId = glGenLists(1);
                    glNewList(cHullId, GL_COMPILE_AND_EXECUTE);

                    glBegin(GL_TRIANGLES);
                    for(int i=1; i< alcx->delcx->DeluanayTrigs.size(); i++){
                        Triangle tri = alcx->delcx->DeluanayTrigs[i];
                        if(tri.Hull){
                            glNormal3d(tri.Normal->X, tri.Normal->Y, tri.Normal->Z);
                            for(int j=1;j<4;j++){
                                Vector3 v = alcx->vertexList[tri.Corners[j]].getCoordVector();
                                glVertex3d(v.X, v.Y, v.Z);
                            }
                        }
                    }
                    glEnd();

                    glEndList();
                } else {
                    glCallList(cHullId);
                }
                if(cHullNorm){
                    if(glIsList(cHullNormId) == GL_FALSE){
                        cHullNormId = glGenLists(1);
                        glNewList(cHullNormId, GL_COMPILE_AND_EXECUTE);

                        //======= Drawing Normals of Convex Hull Tris =================

                        glColor3d(1,0,0);
                        glLineWidth(2);
                        glBegin(GL_LINES);
                        for(int i=1; i< alcx->delcx->DeluanayTrigs.size(); i++){
                            Triangle tri = alcx->delcx->DeluanayTrigs[i];
                            if(tri.Hull){
                                glNormal3d(tri.Normal->X, tri.Normal->Y, tri.Normal->Z);
                                Vector3 center;
                                for(int j=1; j<4; j++){
                                    Vertex v = alcx->vertexList[tri.Corners[j]];
                                    Vector3 coord = v.getCoordVector();
                                    Vector3::Sum(&center, &coord, &center);
                                }
                                Vector3::Scale(&center, &center, 1/3.0);
                                glVertex3d(center.X, center.Y, center.Z);
                                Vector3::Sum(&center, tri.Normal, &center);
                                Vector3::Sum(&center, tri.Normal, &center);
                                glVertex3d(center.X, center.Y, center.Z);
                            }
                        }
                        glEnd();

                        glEndList();
                    } else {
                        glCallList(cHullNormId);
                    }
                }
            }
            glDisable(GL_COLOR_MATERIAL);
            glDisable(GL_NORMALIZE);
        }
}

void Processor::ProcessEpsilonInterval(double epsilon,int lowRank, int highRank)
{
    double real_alpha;

    alcx->MarkRelevantSimplices(epsilon,lowRank,highRank,candidateTrigs,candidateTets);

    for(int i=lowRank;i<=highRank;i++)
    {
        //first find the pockets for the rank
        std::vector<std::vector<int> > pockets;
        if(alcx->MaximumRank != 0)
        {
            alcx->BuildComplex(i,vertexList);
            real_alpha = alcx->AlphaThresholdSqr(i);
            if(i!=0)
            {
                alcx->Adjust(real_alpha,vertexList);
            }
            if(i!=0 && i!=alcx->MaximumRank)
            {
                pocket->FindPockets(alcx->delcx,alcx->sortedTet);
                pockets = pocket->GetPockets();
            }
        }


        //for each rank find if a triangle belongs to different pockets
        //if yes add it refined list
        for(uint j = 0;j<candidateTrigs.size();j++)
        {
            int candidatePoc1 = -1;
            int candidatePoc2 = -1;
            if(alcx->delcx->DeluanayTrigs[candidateTrigs[j].simplex].Hull == 1)//on hull
            {
                continue;
            }

            if(i == candidateTrigs[j].rank)//ranks match
            {
                //FILE *fp2 = fopen("tetcheck.txt","w");
                int tet1 = alcx->delcx->DeluanayTrigs[candidateTrigs[j].simplex].ReverseLink1;
                int tet2 = alcx->delcx->DeluanayTrigs[candidateTrigs[j].simplex].ReverseLink2;

                for(uint _i = 0;_i<pockets.size();_i++)
                {
                    for(uint _j = 0;_j<pockets[_i].size();_j++)
                    {
                        if(tet1 == pockets[_i][_j])
                        {
                            //fprintf(fp2,"trig = %d tet1 = %d candidate poc1 = %d",candidateTrigs[j].simplex,tet1,_i);
                            candidatePoc1 = _i;
                        }
                        if(tet2 == pockets[_i][_j])
                        {
                            //fprintf(fp2,"trig = %d tet2 = %d candidate poc2 = %d",candidateTrigs[j].simplex,tet2,_i);
                            candidatePoc2 = _i;
                        }
                    }
                }
                //fclose (fp2);
            }

            if((candidatePoc1 != candidatePoc2)  && (candidatePoc1 != -1) && (candidatePoc2 != -1))//two tets are indeed in different pockets
            {
                bool flag = false;
                for(uint _i=0;_i<refinedCandidateTrigs.size();_i++)//check for repetition
                {
                    if(candidateTrigs[j].simplex == refinedCandidateTrigs[_i].simplex)
                    {
                        flag = true;
                        break;
                    }
                }
                if(!flag)
                {
                    if(candidateTrigs[j].simplex == 2707){
                        //FILE *fp1 = fopen("trial","w");
                        //fprintf(fp1,"%d",candidateTrigs[j].rank);
                        //fclose(fp1);
                    }
                    SimplexMasterListMap smlm(candidateTrigs[j].simplex,candidateTrigs[j].mlIndex);
                    refinedCandidateTrigs.push_back(smlm);
                }
            }
        }
        //no need to check anything any tet is a destroyer of a void and hence can be added unless its already there in the list
        for(uint j = 0;j<candidateTets.size();j++)
        {
            bool flag = false;
            for(unsigned int _i=0;_i<refinedCandidateTets.size ();_i++)//check for repetitions
            {
                if(candidateTets[j].simplex == refinedCandidateTets[_i].simplex)
                {
                    flag = true;
                    break;
                }
            }
            if(!flag)
            {
                SimplexMasterListMap smlm(candidateTets[j].simplex,candidateTets[j].mlIndex);
                refinedCandidateTets.push_back(smlm);
            }
        }
    }

    //FILE* fp = fopen("candidates","w");
    /*for(uint i = 0;i<refinedCandidateTrigs.size();i++)
    {
        int tet1 = alcx->delcx->DeluanayTrigs[refinedCandidateTrigs[i].simplex].ReverseLink1;
        int tet2 = alcx->delcx->DeluanayTrigs[refinedCandidateTrigs[i].simplex].ReverseLink2;
        //fprintf(fp,"trig = %d tet1 = %d tet2 = %d\n",refinedCandidateTrigs[i].simplex,tet1,tet2);
        //fprintf(fp,"tet1 = %d n1 = %d n2 = %d n3 = %d n4 = %d\n",tet1,alcx->delcx->DeluanayTet[tet1].Neighbours[1],alcx->delcx->DeluanayTet[tet1].Neighbours[2],alcx->delcx->DeluanayTet[tet1].Neighbours[3],alcx->delcx->DeluanayTet[tet1].Neighbours[4]);
        //fprintf(fp,"tet2 = %d n1 = %d n2 = %d n3 = %d n4 = %d\n\n",tet2,alcx->delcx->DeluanayTet[tet2].Neighbours[1],alcx->delcx->DeluanayTet[tet2].Neighbours[2],alcx->delcx->DeluanayTet[tet2].Neighbours[3],alcx->delcx->DeluanayTet[tet2].Neighbours[4]);
    }*/
    //fprintf(fp,"\n");
    /*for(uint i = 0;i<refinedCandidateTets.size();i++)
    {
        fprintf(fp,"%d\n",refinedCandidateTets[i].simplex);
    }*/
    //fclose(fp);
}

void Processor::ModifyFiltration()
{
    alcx->ModifyFiltration(refinedCandidateTrigs,refinedCandidateTets);
}

void Processor::UndoModify()
{
    alcx->UndoModifyFiltration(refinedCandidateTrigs,refinedCandidateTets);
}

int Processor::getNumberOfAtoms ()
{
    return vertexList.size ();
}

int Processor::getNumberOfSimplices ()
{
    return (vertexList.size () + alcx->delcx->DeluanayEdges.size () + alcx->delcx->DeluanayTrigs.size () + alcx->delcx->DeluanayTet.size ());
}
