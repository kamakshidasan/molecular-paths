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

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vector>
#include <vector3.h>
#include <vertex.h>

class Triangle
{
    public:
        int Corners[4];
        int TrigLink[4];
        int ReverseLink1;
        int ReverseLink2;
        int nLink;
        int Rho;
        int Mu1;
        int Mu2;
        int AlphaStatus;
        int ufKey;
        int PocIndex;
        int Hull;
        int Entry;
        int trigCoef;
        bool isValid;
        bool IsAttached;
        bool IsMouth;
        double Size;
        int Repeats;
        int MainRepeats;
        int Persistence;
        double AlphaPersistence;
        Vector3 *Normal;
        std::vector<int> BoundarySimplices;

        Triangle();
        Triangle(unsigned int n[3]);
        ~Triangle();

        Vector3 getCentroid(std::vector <Vertex> &vertlist);
};

#endif // TRIANGLE_H
