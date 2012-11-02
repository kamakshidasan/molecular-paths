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

#include "tetrahedron.h"

#define Sign(X)  (((X) > 0) ? 1 : (((X) < 0) ? -1 : 0))

typedef double Matrix[4][4];

double det4 (Matrix A)
{
        double t0, t1, t2, t3;
        t0 = A[0][0] * (A[1][1] * (A[2][2] * A[3][3] - A[3][2] * A[2][3]) -
                        A[2][1] * (A[1][2] * A[3][3] - A[3][2] * A[1][3]) +
                        A[3][1] * (A[1][2] * A[2][3] - A[2][2] * A[1][3]));

        t1 = A[1][0] * (A[0][1] * (A[2][2] * A[3][3] - A[3][2] * A[2][3]) -
                        A[2][1] * (A[0][2] * A[3][3] - A[3][2] * A[0][3]) +
                        A[3][1] * (A[0][2] * A[2][3] - A[2][2] * A[0][3]));

        t2 = A[2][0] * (A[0][1] * (A[1][2] * A[3][3] - A[3][2] * A[1][3]) -
                        A[1][1] * (A[0][2] * A[3][3] - A[3][2] * A[0][3]) +
                        A[3][1] * (A[0][2] * A[1][3] - A[1][2] * A[0][3]));

        t3 = A[3][0] * (A[0][1] * (A[1][2] * A[2][3] - A[2][2] * A[1][3]) -
                        A[1][1] * (A[0][2] * A[2][3] - A[2][2] * A[0][3]) +
                        A[2][1] * (A[0][2] * A[1][3] - A[1][2] * A[0][3]));

        return (t0 - t1 + t2 - t3);
}

Tetrahedron::Tetrahedron()
{
        Kill();
        for (int i = 0; i < 5; i++)
        {
                this->Corners[i] = 0;
                this->Neighbours[i] = 0;
                this->Nindex[i] = 0;
                this->TetLink[i] = -1;
                this->TetFlow[i] = 0;
        }
        this->AlphaStatus = -1;
        this->Index = -1;
        this->Orient = 0;
        this->Hull = -1;
        this->Depth = 0;
        this->Rho = 0;
        this->Size = -1.0;
        this->ufKey = -1;
        this->PocIndex = -1;
        this->posIndex = -1;
        this->Repeats = -1;
        this->MainRepeats = -1;
        this->Persistence = -1;
        this->isValid = false;
}


Tetrahedron::~Tetrahedron()
{
}

/*!
    \fn Tetrahedron::Kill()
 */
void Tetrahedron::Kill()
{
        this->Status = 0;
}


/*!
    \fn Tetrahedron::center4(std::vector<Vertex> &vertexList,std::vector<Vertex> &slist)
 */
void Tetrahedron::center4(std::vector<Vertex> &vertexList,std::vector<Vertex> &slist)
{
        Vertex vert = powerVert(vertexList);
        slist.push_back(vert);
}

Vertex Tetrahedron::powerVert(std::vector<Vertex> &vertexList){
    Matrix M;
    int i,j,k,l;
    double  i0, j0, k0, l0;
    double  i4, j4, k4, l4;
    double D0, Dx, Dy, Dz;
    double v[3] = {0,0,0},W[4];
    double dx2[4],dy2[4],dz2[4],wi2[4];

    i = this->Corners[1];
    j = this->Corners[2];
    k = this->Corners[3];
    l = this->Corners[4];

    i4 = Sign(vertexList[i].BackRadius) * vertexList[i].BackRadius * vertexList[i].BackRadius;
    j4 = Sign(vertexList[j].BackRadius) * vertexList[j].BackRadius * vertexList[j].BackRadius;
    k4 = Sign(vertexList[k].BackRadius) * vertexList[k].BackRadius * vertexList[k].BackRadius;
    l4 = Sign(vertexList[l].BackRadius) * vertexList[l].BackRadius * vertexList[l].BackRadius;

    i0 = (1.0/2.0)*(i4 - pow(vertexList[i].Coordinates[1],2.0) - pow(vertexList[i].Coordinates[2],2.0) - pow(vertexList[i].Coordinates[3],2.0));
    j0 = (1.0/2.0)*(j4 - pow(vertexList[j].Coordinates[1],2.0) - pow(vertexList[j].Coordinates[2],2.0) - pow(vertexList[j].Coordinates[3],2.0));
    k0 = (1.0/2.0)*(k4 - pow(vertexList[k].Coordinates[1],2.0) - pow(vertexList[k].Coordinates[2],2.0) - pow(vertexList[k].Coordinates[3],2.0));
    l0 = (1.0/2.0)*(l4 - pow(vertexList[l].Coordinates[1],2.0) - pow(vertexList[l].Coordinates[2],2.0) - pow(vertexList[l].Coordinates[3],2.0));

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = vertexList[l].Coordinates[1]; M[3][1] = vertexList[l].Coordinates[2]; M[3][2] = vertexList[l].Coordinates[3]; M[3][3] = 1.0;
    D0 = det4 (M);

    if (fabs(D0) < 0.00001)
            printf("center4 fabs(D0) = %lf",fabs(D0));

    M[0][0] = -i0; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = -j0; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = -k0; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = -l0; M[3][1] = vertexList[l].Coordinates[2]; M[3][2] = vertexList[l].Coordinates[3]; M[3][3] = 1.0;
    Dx = det4 (M);

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = -i0; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = -j0; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = -k0; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = vertexList[l].Coordinates[1]; M[3][1] = -l0; M[3][2] = vertexList[l].Coordinates[3]; M[3][3] = 1.0;
    Dy = det4 (M);

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = -i0; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = -j0; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = -k0; M[2][3] = 1.0;
    M[3][0] = vertexList[l].Coordinates[1]; M[3][1] = vertexList[l].Coordinates[2]; M[3][2] = -l0; M[3][3] = 1.0;
    Dz = det4 (M);

    //center
    v[0] = Dx/D0;  v[1] = Dy/D0;  v[2] = Dz/D0;

    dx2[0] = pow((v[0] - vertexList[i].Coordinates[1]),2.0);
    dy2[0] = pow((v[1] - vertexList[i].Coordinates[2]),2.0);
    dz2[0] = pow((v[2] - vertexList[i].Coordinates[3]),2.0);
    wi2[0] = pow(vertexList[i].Radius,2.0);

    dx2[1] = pow((v[0] - vertexList[j].Coordinates[1]),2.0);
    dy2[1] = pow((v[1] - vertexList[j].Coordinates[2]),2.0);
    dz2[1] = pow((v[2] - vertexList[j].Coordinates[3]),2.0);
    wi2[1] = pow(vertexList[j].Radius,2.0);

    dx2[2] = pow((v[0] - vertexList[k].Coordinates[1]),2.0);
    dy2[2] = pow((v[1] - vertexList[k].Coordinates[2]),2.0);
    dz2[2] = pow((v[2] - vertexList[k].Coordinates[3]),2.0);
    wi2[2] = pow(vertexList[k].Radius,2.0);

    dx2[3] = pow((v[0] - vertexList[l].Coordinates[1]),2.0);
    dy2[3] = pow((v[1] - vertexList[l].Coordinates[2]),2.0);
    dz2[3] = pow((v[2] - vertexList[l].Coordinates[3]),2.0);
    wi2[3] = pow(vertexList[l].Radius,2.0);

    for(int m = 0;m<4;m++)
    {
            W[m] = dx2[m] + dy2[m] + dz2[m] - wi2[m];
    }

    if(W[0]<0.0)
    {
            W[0] = 0.0;
    }
    W[0] = sqrt(W[0]);

    Vertex vert(v[0],v[1],v[2],W[0],-1,0.0);
    return vert;
}
