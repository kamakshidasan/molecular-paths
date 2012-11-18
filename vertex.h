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
#ifndef VERTEX_H
#define VERTEX_H

#include <gmp.h>
#include <vector>
#include <vector3.h>

class Vertex
{
    public:
        double Radius;
        double BackRadius;
        double Coordinates[4];
        double NormCoordinates[4];
        double BackWeight;
        double Weight;
        double ranValue;
        double Size;
        double Coef;
        bool valid;

        int AlphaStatus;
        int Index;
        int Redinfo;
        int Mu2;
        int Hull;
        int Mu1;
        int Rho;
        int Repeats;
        int ufKey;
        mpz_t V[5];

        int selected;

        Vertex(){
        }

        Vertex(double x, double y, double z, double radius,int index,double scale);
        ~Vertex();

        static bool ranSortFunc(const Vertex& lhs,const Vertex& rhs);
        bool operator < (const Vertex& rhs ) const;

        Vector3 getCoordVector(){
            return Vector3(Coordinates[1], Coordinates[2], Coordinates[3]);
        }
};

#endif // VERTEX_H
