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

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include <vector>
#include "vertex.h"
#include <cmath>

class Tetrahedron
{
    public:
        int Corners[5];
        int Neighbours[5];
        int Nindex[5];
        int TetLink[5];
        int TetFlow[5];
        int Rho;
        double Size;
        int ufKey;
        int PocIndex;

        int AlphaStatus;
        int Status;
        int Index;
        int Orient;
        int Hull;
        int Depth;
        int posIndex;
        int Repeats;
        int MainRepeats;
        int Persistence;

        bool isValid;

        std::vector<int> BoundarySimplices;

        Tetrahedron();
        ~Tetrahedron();

        void Kill();
        void center4(std::vector<Vertex> &vertexList,std::vector<Vertex> &slist);
        Vertex powerVert(std::vector<Vertex> &vertexList);
};

#endif // TETRAHEDRON_H
