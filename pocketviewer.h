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
#ifndef POCKETVIEWER_H
#define POCKETVIEWER_H

#include "QGLViewer/qglviewer.h"
#include <processor.h>

class PocketViewer : public QGLViewer,public Processor
{
    Q_OBJECT
    public:
        PocketViewer(Processor *pr,QWidget* parent=NULL);
        void setRank(int rank);
        void setPersistence(int persistence);
        void setAlphaShape();
        void setWireFrame();
        void setAllPockets();
        void setMouths();
        void setOnlyPockets();
        void setOnlyVoids();

        void setPowerDiag();
        void setComplementSpacePD();
        void setOnlyInsideVerts();
        void setPruneIsolatedVerts();
        void setIntersectEdges();

        void setCHull();
        void setCHullWF();
        void setCHullNorm();

        void setPocketSkinSurface();
        void setPocketSkinWireFrame();
        void setAlphaSkinSurface();
        void setAlphaSkinWireFrame();
        void setSmoothShading(bool flag);
        void setAllOrIndivPocs(bool flag);
        void setCurrentPocNum(int num);
        void setPocketWireFrame();

    protected :
        virtual void draw();
        virtual void init();
        virtual QString helpString() const;

    private:
        Processor *m_processor;
        int rank;
        int persistence;
        int currentPocNum;
        bool alphaShape;
        bool wireFrame;
        bool allPockets;
        bool mouths;
        bool onlyPockets;
        bool onlyVoids;

        bool powerDiag;
        bool complementSpacePD;
        bool onlyInsideVerts;
        bool pruneIsolatedVerts;
        bool intersectEdges;

        bool cHull;
        bool cHullWF;
        bool cHullNorm;

        bool pocketSkinSurface;
        bool pocketSkinWireFrame;
        bool alphaSkinSurface;
        bool alphaSkinWireFrame;
        bool smoothShading;
        bool allOrIndivPocs;
        bool volumeEnabled;
        bool pocketWireframe;
};

#endif // POCKETVIEWER_H
