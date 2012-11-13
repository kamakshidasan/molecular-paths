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
#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <QTableWidget>
#include <QString>
#include <QTableWidgetItem>
#include <QTableWidget>
#include <QLineEdit>
#include <QSpinBox>
#include <deluanaycomplex.h>
#include <pocket.h>
#include <mouth.h>
#include <alphacomplex.h>
#include <volume.h>
#include <vertex.h>
#include <filereader.h>
#include <rankmap.h>
#include <simplexmasterlistmap.h>
#include <proteinrenderer.h>

class Processor
{
    //private:
    public:
        std::vector<Vertex> vertexList;
        std::vector<int> PocketNMouths;
        AlphaComplex *alcx;
        PowerDiagram* powerDiagram;
        ProteinRenderer * proteinRenderer;
        FileReader *fr;
        Pocket *pocket;
        Mouth *mouth;
        Volume *volume;
        double scale;
        bool isRenderable;

        std::vector<SimplexMasterListMap> refinedCandidateTrigs;
        std::vector<SimplexMasterListMap> refinedCandidateTets;
        
        std::vector<RankMap> candidateTrigs;
        std::vector<RankMap> candidateTets;

        void MapPocketsToMouths();
        void FillTables();

    //public:
        QTableWidget *table1;
        QTableWidget *table2;

        QLineEdit *totVol;
        QLineEdit *totSurf;
        QLineEdit *pocVol;
        QLineEdit *pocSurf;
        QLineEdit *mouSurf;

        QSpinBox *pspin;

        Processor();
        ~Processor();

        bool IsRenderable();

        void read(const char *filename,double centre[],double *size, bool constantRadius, bool incrementRadius);
        int getMaxRank();
        int getMaxPersistence();
        int getRankForAlpha(double alphavalue);
        double getAlphaValue(int rank);
        void CalculateEverythingFor(int rank, int persistence);
        void CalculateRelevant(int rank, int persistence);
        void CalculateVolumes(int rank, int persistence);
        void Render(int persistence,bool alphaShape,bool allPockets,bool onlyPockets,bool onlyVoids,bool skinSurface,bool mouths,bool allindflag,int pnum,int rank,bool wireFrame,
                    bool alphaSkinSurface,bool alphaSkinWireFrame,bool smoothShading,bool skinWireFrame,bool pocketWireFrame, bool powerDiag,
                    bool cHull, bool cHullWF, bool cHullNorm, bool showPath, bool showSpaceFill);

        void ProcessEpsilonInterval(double epsilon,int lowRank,int highRank);
        void ModifyFiltration();
        void UndoModify();
        int getNumberOfAtoms();
        int getNumberOfSimplices();
};

#endif // PROCESSOR_H
