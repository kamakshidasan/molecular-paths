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

#include "metric.h"

#define Sign(X)  (((X) > 0) ? 1 : (((X) < 0) ? -1 : 0))

Metric::Metric(std::vector<Vertex> &vertexList,int vcount)
{
    PI = 4 * atan(1);
    Radius2.reserve (vcount);
    for (unsigned int i = 1; i < vertexList.size(); i++)
    {
            double r = vertexList[i].Radius;
            Radius2[i] = pow(r, 2.0);
    }
    size = new Size();
    sos = new Sos();
}

Metric::~Metric()
{
    Radius2.clear ();
}

double Metric::Det2(Matrix A)
{
    return A[0][0] * A[1][1] - A[1][0] * A[0][1];
}

double Metric::Det3(Matrix A)
{
    return A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) -
           A[1][0] * (A[0][1] * A[2][2] - A[2][1] * A[0][2]) +
           A[2][0] * (A[0][1] * A[1][2] - A[1][1] * A[0][2]);
}

double Metric::Det4(Matrix A)
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

    return t0 - t1 + t2 - t3;
}

int Metric::Ccw(std::vector<Vertex> vertexList,int i,int j,int k,int l)
{
    int a;
    sos->Minor4 (vertexList,&i,&j,&k,&l,&a);
    return a;
}

void Swap(int *a,int *b)
{
    int temp;

    temp = *a;
    *a = *b;
    *b = temp;
}

Vector3 Metric::Center2(std::vector<Vertex> &vertexList, int i, int j)
{
    Vector3 D,sS,sT,sum;
    double aux,lambda,i4,j4;
    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3::DiffVector (&D,&s,&t);
    Vector3::DotProduct (&D,&D,&aux);

    i4 = Sign(vertexList[i].Radius) * vertexList[i].Radius * vertexList[i].Radius;
    j4 = Sign(vertexList[j].Radius) * vertexList[j].Radius * vertexList[j].Radius;

    lambda = (1.0/2.0) - (i4 - j4)/(2 * aux);

    Vector3::Scale (&sS,&s,lambda);
    Vector3::Scale (&sT,&t,(1 - lambda));

    Vector3::Sum (&sum,&sS,&sT);
    return sum;
}

Vector3 Metric::Center3(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    Matrix M;
    double i0, j0, k0;
    double i4, j4, k4;
    Vector3 Y;
    double A1, A2, A3, A4;
    double D0, Dx, Dy, Dz;
    M[0][0] = vertexList[i].Coordinates[2]; M[0][1] = vertexList[i].Coordinates[3]; M[0][2] = 1.0;
    M[1][0] = vertexList[j].Coordinates[2]; M[1][1] = vertexList[j].Coordinates[3]; M[1][2] = 1.0;
    M[2][0] = vertexList[k].Coordinates[2]; M[2][1] = vertexList[k].Coordinates[3]; M[2][2] = 1.0;
    A1 = Det3 (M);

    M[0][0] = vertexList[i].Coordinates[3]; M[0][1] = vertexList[i].Coordinates[1]; M[0][2] = 1.0;
    M[1][0] = vertexList[j].Coordinates[3]; M[1][1] = vertexList[j].Coordinates[1]; M[1][2] = 1.0;
    M[2][0] = vertexList[k].Coordinates[3]; M[2][1] = vertexList[k].Coordinates[1]; M[2][2] = 1.0;
    A2 = Det3 (M);

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = 1.0;
    A3 = Det3 (M);

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = vertexList[i].Coordinates[3];
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = vertexList[j].Coordinates[3];
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = vertexList[k].Coordinates[3];
    A4 = Det3 (M);

    i4 = Sign(vertexList[i].Radius) * vertexList[i].Radius * vertexList[i].Radius;
    j4 = Sign(vertexList[j].Radius) * vertexList[j].Radius * vertexList[j].Radius;
    k4 = Sign(vertexList[k].Radius) * vertexList[k].Radius * vertexList[k].Radius;

    i0 = (1.0/2.0)*(i4 - pow(vertexList[i].Coordinates[1],2.0) - pow(vertexList[i].Coordinates[2],2.0) - pow(vertexList[i].Coordinates[3],2.0));
    j0 = (1.0/2.0)*(j4 - pow(vertexList[j].Coordinates[1],2.0) - pow(vertexList[j].Coordinates[2],2.0) - pow(vertexList[j].Coordinates[3],2.0));
    k0 = (1.0/2.0)*(k4 - pow(vertexList[k].Coordinates[1],2.0) - pow(vertexList[k].Coordinates[2],2.0) - pow(vertexList[k].Coordinates[3],2.0));

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = A1;                           M[3][1] = A2;                           M[3][2] = A3;                           M[3][3] = 0.0;

    D0 = Det4 (M);
    //if (fabs(D0) < EPSILON)
    //  volbl_warning(fabs(D0), "center3", "fabs(D0)");

    M[0][0] = -i0; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = -j0; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = -k0; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = A4;  M[3][1] = A2;                           M[3][2] = A3;                           M[3][3] = 0.0;
    Dx = Det4 (M);

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = -i0; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = -j0; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = -k0; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = A1;                           M[3][1] = A4;  M[3][2] = A3;                           M[3][3] = 0.0;
    Dy = Det4 (M);

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = -i0; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = -j0; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = -k0; M[2][3] = 1.0;
    M[3][0] = A1;                           M[3][1] = A2;                           M[3][2] = A4;  M[3][3] = 0.0;
    Dz = Det4 (M);

    Y.X = Dx/D0;  Y.Y = Dy/D0;  Y.Z = Dz/D0;
    return Y;
}

Vector3 Metric::Center4(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    Matrix M;
    Vector3 V;
    double i0, j0, k0, l0;
    double i4, j4, k4, l4;
    double D0, Dx, Dy, Dz;

    i4 = Sign(vertexList[i].Radius) * vertexList[i].Radius * vertexList[i].Radius;
    j4 = Sign(vertexList[j].Radius) * vertexList[j].Radius * vertexList[j].Radius;
    k4 = Sign(vertexList[k].Radius) * vertexList[k].Radius * vertexList[k].Radius;
    l4 = Sign(vertexList[l].Radius) * vertexList[l].Radius * vertexList[l].Radius;

    i0 = (1.0/2.0)*(i4 - pow(vertexList[i].Coordinates[1],2.0) - pow(vertexList[i].Coordinates[2],2.0) - pow(vertexList[i].Coordinates[3],2.0));
    j0 = (1.0/2.0)*(j4 - pow(vertexList[j].Coordinates[1],2.0) - pow(vertexList[j].Coordinates[2],2.0) - pow(vertexList[j].Coordinates[3],2.0));
    k0 = (1.0/2.0)*(k4 - pow(vertexList[k].Coordinates[1],2.0) - pow(vertexList[k].Coordinates[2],2.0) - pow(vertexList[k].Coordinates[3],2.0));
    l0 = (1.0/2.0)*(l4 - pow(vertexList[l].Coordinates[1],2.0) - pow(vertexList[l].Coordinates[2],2.0) - pow(vertexList[l].Coordinates[3],2.0));

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = vertexList[l].Coordinates[1]; M[3][1] = vertexList[l].Coordinates[2]; M[3][2] = vertexList[l].Coordinates[3]; M[3][3] = 1.0;
    D0 = Det4 (M);

    //if (fabs(D0) < EPSILON)
    //volbl_warning(fabs(D0), "center4", "fabs(D0)");

    M[0][0] = -i0; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = -j0; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = -k0; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = -l0; M[3][1] = vertexList[l].Coordinates[2]; M[3][2] = vertexList[l].Coordinates[3]; M[3][3] = 1.0;
    Dx = Det4 (M);

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = -i0; M[0][2] = vertexList[i].Coordinates[3]; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = -j0; M[1][2] = vertexList[j].Coordinates[3]; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = -k0; M[2][2] = vertexList[k].Coordinates[3]; M[2][3] = 1.0;
    M[3][0] = vertexList[l].Coordinates[1]; M[3][1] = -l0; M[3][2] = vertexList[l].Coordinates[3]; M[3][3] = 1.0;
    Dy = Det4 (M);

    M[0][0] = vertexList[i].Coordinates[1]; M[0][1] = vertexList[i].Coordinates[2]; M[0][2] = -i0; M[0][3] = 1.0;
    M[1][0] = vertexList[j].Coordinates[1]; M[1][1] = vertexList[j].Coordinates[2]; M[1][2] = -j0; M[1][3] = 1.0;
    M[2][0] = vertexList[k].Coordinates[1]; M[2][1] = vertexList[k].Coordinates[2]; M[2][2] = -k0; M[2][3] = 1.0;
    M[3][0] = vertexList[l].Coordinates[1]; M[3][1] = vertexList[l].Coordinates[2]; M[3][2] = -l0; M[3][3] = 1.0;
    Dz = Det4 (M);

    V.X = Dx/D0; V.Y = Dy/D0; V.Z = Dz/D0;

    return V;
}

Vector3 Metric::TriangleDual(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    Vector3 Y, N, D,D1,D2,diff1,diff2,res1,res;
    double xi, r, S1, S2, S3, aux1,temp1,temp2;

    Y = Center3 (vertexList,i,j,k);
    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);

    Vector3::DiffVector (&diff1,&t,&s);
    Vector3::DiffVector (&diff2,&u,&s);

    Vector3::CrossProduct (&N,&diff1,&diff2);

    Vector3::DiffVector (&D,&Y,&s);
    Vector3::DiffVector (&D1,&Y,&t);
    Vector3::DiffVector (&D2,&Y,&u);

    Vector3::DotProduct (&D,&N,&S1);
    Vector3::DotProduct (&N,&N,&S2);
    Vector3::DotProduct (&D,&D,&S3);

    Vector3::DotProduct (&D1,&D1,&temp1);
    Vector3::DotProduct (&D2,&D2,&temp2);

    r = vertexList[i].Radius;
    aux1 = S1*S1 - S3*S2 + r*r*S2;

    //aux1 should be non negative
    if(aux1 < 0.0) aux1 = 0.0;

    xi = (-S1 + sqrt(aux1)) / S2;

    //xi also should be non negative
    if(xi < 0.0) xi = 0.0;

    Vector3::Scale (&res1,&N,xi);
    Vector3::Sum (&res,&Y,&res1);

    return res;
}

double Metric::AngleDihedral(Vector3 s, Vector3 t, Vector3 u, Vector3 v)
{
    Vector3 Mu, Mv, Nu, Nv;
    Vector3 r1, r2;
    double spu, spv, aux, phi;

    Vector3::DiffVector (&r1,&u,&s);
    Vector3::DiffVector (&r2,&u,&t);

    Vector3::CrossProduct (&Mu,&r1,&r2);

    Vector3::DiffVector (&r1,&v,&s);
    Vector3::DiffVector (&r2,&v,&t);

    Vector3::CrossProduct (&Mv,&r1,&r2);

    Vector3::DotProduct (&Mu,&Mu,&spu);
    Vector3::DotProduct (&Mv,&Mv,&spv); 

    Vector3::Scale (&Nu,&Mu,(1.0/sqrt (spu)));
    Vector3::Scale (&Nv,&Mv,(1.0/sqrt (spv)));

    Vector3::DotProduct (&Nu,&Nv,&aux);

    /* Assert_always(aux >= -1.0 and aux <= 1.0) */
    if (aux < -1.0)
    {
        aux = -1.0;
    }
    if (aux > 1.0)
    {
        aux = 1.0;
    }

    phi = acos(aux) / (2.0 * PI);

    return phi;
}

double Metric::AngleSolid(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    double phi_t, phi_u, phi_v;

    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);
    Vector3 v(vertexList[l].Coordinates[1],vertexList[l].Coordinates[2],vertexList[l].Coordinates[3]);

    phi_t = AngleDihedral(s,t,v,u);
    phi_u = AngleDihedral(s,u,t,v);
    phi_v = AngleDihedral(s,v,u,t);

    return (1.0/2.0) * (phi_t + phi_u + phi_v) - (1.0/4.0);
}

double Metric::DiskRadius(std::vector<Vertex> &vertexList, int i, int j)
{
    double h, aux;
    h = CapHeight (vertexList,i,j);
    aux = h * (2.0 * vertexList[i].Radius - h);
    //aux should be non negative
    if(aux < 0.0) aux = 0.0;

    return (sqrt (aux));
}

double Metric::DiskLength(std::vector<Vertex> &vertexList, int i, int j)
{
    return (2.0 * PI * DiskRadius (vertexList,i,j));
}

double Metric::DiskArea(std::vector<Vertex> &vertexList, int i, int j)
{
    return ((1.0/2.0) * DiskRadius (vertexList,i,j) * DiskLength (vertexList,i,j));
}

double Metric::SegmentHeight(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    Vector3 Y3,Y2;
    double h,rho,dist;
    bool isAttached;

    Y3 = Center3 (vertexList,i,j,k);
    Y2 = Center2 (vertexList,i,j);

    rho = DiskRadius (vertexList,i,j);

    Vector3::Distance (&dist, &Y2, &Y3);
    size->CheckEdge (vertexList[i].V,vertexList[j].V,vertexList[k].V,&isAttached);
    if(isAttached)
    {
        h = rho + dist;
        if(h > 2 * rho)
        {
            //h can never be greater than 2 * rho
            h = 2 * rho;
        }
    }
    else
    {
        h = rho - dist;
        if(h < 0.0)
        {
            //h can never be less than 0
            h = 0;
        }
    }
    return h;
}

double Metric::SegmentAngle(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    Vector3 Pjk, Pkj;
    double delta,ang1,ang2;

    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);

    Pjk = TriangleDual (vertexList,i,j,k);
    Pkj = TriangleDual (vertexList,i,k,j);

    ang1 = AngleDihedral (s,t,u,Pjk);
    ang2 = AngleDihedral (s,t,u,Pkj);

    delta = fabs(ang1 - ang2);

    //if(delta < less epsilon) TODO

    return ang1 + ang2;
}

double Metric::SegmentLength(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    return (SegmentAngle (vertexList,i,j,k) * DiskLength (vertexList,i,j));
}

double Metric::SegmentArea(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    double H,S,T,dist;
    Vector3 Pjk, Pkj;

    S = (1.0/2.0) * SegmentLength (vertexList,i,j,k) * DiskRadius (vertexList,i,j);

    Pjk = TriangleDual (vertexList,i,j,k);
    Pkj = TriangleDual (vertexList,i,k,j);
    Vector3::Distance (&dist,&Pjk,&Pkj);

    H = DiskRadius (vertexList,i,j) - SegmentHeight (vertexList,i,j,k);
    T = (1.0/2.0) * H * dist;

    return S - T;
}

double Metric::Segment2Angle(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    Vector3 Pjl,Pkj;
    double a;

    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);
    Vector3 v(vertexList[l].Coordinates[1],vertexList[l].Coordinates[2],vertexList[l].Coordinates[3]);

    Pjl = TriangleDual (vertexList,i,j,l);
    Pkj = TriangleDual (vertexList,i,k,j);

    a = AngleDihedral (s,t,u,Pkj) + AngleDihedral (s,t,u,Pjl) - AngleDihedral (s,t,u,v);

    return a;
}

double Metric::Segment2Length(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    return Segment2Angle (vertexList,i,j,k,l) * DiskLength (vertexList,i,j);
}

double Metric::Segment2Area(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    Vector3 Pjk,Pkj,Pjl,Plj,Y;
    double h_k,h_l,r_ij,S,Tk,Tl,dist1,dist2;

    if(! Ccw (vertexList,i,j,k,l)) Swap (&k,&l);

    Pjk = TriangleDual (vertexList,i,j,k);
    Pkj = TriangleDual (vertexList,i,k,j);
    Pjl = TriangleDual (vertexList,i,j,l);
    Plj = TriangleDual (vertexList,i,l,j);

    Y = Center4 (vertexList,i,j,k,l);

    h_k = SegmentHeight (vertexList,i,j,k);
    h_l = SegmentLength (vertexList,i,j,l);

    r_ij = DiskRadius (vertexList,i,j);

    Vector3::Distance (&dist1,&Pkj,&Y);
    Vector3::Distance (&dist2,&Pjl,&Y);

    S = (1.0/2.0) * r_ij * Segment2Length (vertexList,i,j,k,l);
    Tk = (1.0/2.0) * (r_ij - h_k) * dist1;
    Tl = (1.0/2.0) * (r_ij - h_l) * dist2;

    return S - Tk - Tl;
}

double Metric::CapHeight(std::vector<Vertex> &vertexList, int i, int j)
{
    Vector3 Y;
    double dist;
    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);

    Y = Center2 (vertexList,i,j);

    Vector3::Distance (&dist,&s,&Y);
    if(size->CheckVertex (vertexList,i,j))
        return (vertexList[i].Radius + dist);
    else
        return (vertexList[i].Radius - dist);
}

double Metric::CapArea(std::vector<Vertex> &vertexList, int i, int j)
{
    return (2.0 * PI * vertexList[i].Radius * CapHeight (vertexList,i,j));
}

double Metric::CapVolume(std::vector<Vertex> &vertexList, int i, int j)
{
    double r, S, C;

    r = vertexList[i].Radius;
    S = (1.0/3.0) * r * CapArea (vertexList,i,j);
    C = (1.0/3.0) * (r - CapHeight (vertexList,i,j)) * DiskArea (vertexList,i,j);

    return S - C;
}

double Metric::Cap2Volume(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    double r, S2, Cj, Ck;

    r = vertexList[i].Radius;
    S2 = (1.0/3.0) * r * Cap2Area (vertexList,i,j,k);
    Cj = (1.0/3.0) * (r - CapHeight (vertexList,i,j)) * SegmentArea (vertexList,i,j,k);
    Ck = (1.0/3.0) * (r - CapHeight (vertexList,i,k)) * SegmentArea (vertexList,i,k,j);

    return S2 - Cj - Ck;
}

double Metric::Cap2Area(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    Vector3 Pjk,Pkj;
    double r, A1,A2,A3;
    double l_j, l_k, phi_jk, phi_kj;

    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);

    Pjk = TriangleDual (vertexList,i,j,k);
    Pkj = TriangleDual (vertexList,i,k,j);

    l_j = SegmentAngle (vertexList,i,j,k);
    l_k = SegmentAngle (vertexList,i,k,j);

    phi_jk = (1.0/2.0) - AngleDihedral (s,Pjk,t,u);
    phi_kj = (1.0/2.0) - AngleDihedral (s,Pkj,u,t);

    //if (fabs(phi_jk - phi_kj) > EPSILON) correction TODO

    r = vertexList[i].Radius;
    A1 = (1.0/2.0) * BallArea (i) * (phi_jk + phi_kj);
    A2 = 2.0 * PI * r * l_j * (r - CapHeight (vertexList,i,j));
    A3 = 2.0 * PI * r * l_k * (r - CapHeight (vertexList,i,k));

    return A1 - A2 - A3;
}

double Metric::Cap3Volume(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    double r, S3, Cj, Ck, Cl;

    r = vertexList[i].Radius;
    S3 = (1.0/3.0) * r * Cap3Area (vertexList,i,j,k,l);
    Cj = (1.0/3.0) * (r - CapHeight (vertexList,i,j)) * Segment2Area (vertexList,i,j,k,l);
    Ck = (1.0/3.0) * (r - CapHeight (vertexList,i,k)) * Segment2Area (vertexList,i,k,j,l);
    Cl = (1.0/3.0) * (r - CapHeight (vertexList,i,l)) * Segment2Area (vertexList,i,l,j,k);

    return S3 - Cj - Ck - Cl;
}

double Metric::Cap3Area(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{

    Vector3 Pkj, Plk, Pjl;
    double l_j, l_k, l_l;
    double phi_kj, phi_lk, phi_jl;
    double r, A1, A2, A3, A4;

    if (! Ccw (vertexList,i,j,k,l)) Swap (&k,&l);

    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);
    Vector3 v(vertexList[l].Coordinates[1],vertexList[l].Coordinates[2],vertexList[l].Coordinates[3]);

    Pkj = TriangleDual (vertexList,i, k, j);
    Plk = TriangleDual (vertexList,i, l, k);
    Pjl = TriangleDual (vertexList,i, j, l);

    l_j = Segment2Angle (vertexList,i,j,k,l);
    l_k = Segment2Angle (vertexList,i,k,l,j);
    l_l = Segment2Angle (vertexList,i,l,j,k);


    phi_kj = (1.0/2.0) - AngleDihedral (s,Pkj,u,t);
    phi_lk = (1.0/2.0) - AngleDihedral (s,Plk,v,u);
    phi_jl = (1.0/2.0) - AngleDihedral (s,Pjl,t,v);

    r = vertexList[i].Radius;
    A1 = (1.0/2.0) * BallArea (i) * (phi_kj + phi_lk + phi_jl - (1.0/2.0));
    A2 = 2.0 * PI * r * l_j * (r - CapHeight (vertexList,i,j));
    A3 = 2.0 * PI * r * l_k * (r - CapHeight (vertexList,i,k));
    A4 = 2.0 * PI * r * l_l * (r - CapHeight (vertexList,i,l));

    return A1 - A2 - A3 - A4;
}

double Metric::SectorArea(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    return AngleSolid (vertexList,i,j,k,l) * BallArea (i);
}

double Metric::SectorVolume(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    return AngleSolid (vertexList,i,j,k,l) * BallVolume (vertexList,i);
}

double Metric::WedgeLength(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);
    Vector3 v(vertexList[l].Coordinates[1],vertexList[l].Coordinates[2],vertexList[l].Coordinates[3]);

    return AngleDihedral (s,t,u,v) * DiskLength (vertexList,i,j);
}

double Metric::WedgeArea(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);
    Vector3 v(vertexList[l].Coordinates[1],vertexList[l].Coordinates[2],vertexList[l].Coordinates[3]);

    return AngleDihedral (s,t,u,v) * Ball2Area (vertexList,i,j);
}

double Metric::WedgeVolume(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    Vector3 s(vertexList[i].Coordinates[1],vertexList[i].Coordinates[2],vertexList[i].Coordinates[3]);
    Vector3 t(vertexList[j].Coordinates[1],vertexList[j].Coordinates[2],vertexList[j].Coordinates[3]);
    Vector3 u(vertexList[k].Coordinates[1],vertexList[k].Coordinates[2],vertexList[k].Coordinates[3]);
    Vector3 v(vertexList[l].Coordinates[1],vertexList[l].Coordinates[2],vertexList[l].Coordinates[3]);

    return AngleDihedral(s,t,u,v) * Ball2Volume (vertexList,i,j);
}

double Metric::PawnLength(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    return (1.0/2.0) * Ball3Length (vertexList,i,j,k);
}

double Metric::PawnArea(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    return (1.0/2.0) * Ball3Area (vertexList,i,j,k);
}

double Metric::PawnVolume(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    return (1.0/2.0) * Ball3Volume (vertexList,i,j,k);
}

double Metric::BallVolume(std::vector<Vertex> &vertexList, int idx)
{
    return (1.0/3.0) * vertexList[idx].Radius * BallArea (idx);
}

double Metric::BallArea(int idx)
{
    return 4.0 * PI * Radius2[idx];
}

double Metric::Ball2Length(std::vector<Vertex> &vertexList, int i, int j)
{
    return DiskLength (vertexList,i,j);
}

double Metric::Ball2Volume(std::vector<Vertex> &vertexList, int i, int j)
{
    return CapVolume (vertexList,i,j) + CapVolume (vertexList,j,i);
}

double Metric::Ball2Area(std::vector<Vertex> &vertexList, int i, int j)
{
    return CapArea (vertexList,i,j) + CapArea (vertexList,j,i);
}

double Metric::Ball3Length(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    return SegmentLength (vertexList,i,j,k) + SegmentLength (vertexList,i,k,j) + SegmentLength (vertexList,j,k,i);
}

double Metric::Ball3Volume(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    return Cap2Volume (vertexList,i,j,k) + Cap2Volume (vertexList,j,i,k) + Cap2Volume (vertexList,k,i,j);
}

double Metric::Ball3Area(std::vector<Vertex> &vertexList, int i, int j, int k)
{
    return Cap2Area (vertexList,i,j,k) + Cap2Area (vertexList,j,i,k) + Cap2Area (vertexList,k,i,j);
}

double Metric::Ball4Length(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    return Segment2Length (vertexList,i,j,k,l) + Segment2Length (vertexList,i,k,l,j) +
           Segment2Length (vertexList,i,l,j,k) + Segment2Length (vertexList,j,k,i,l) +
           Segment2Length (vertexList,j,l,k,i) + Segment2Length (vertexList,k,l,i,j);
}

double Metric::Ball4Volume(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    return Cap3Volume (vertexList,i,j,k,l) + Cap3Volume (vertexList,j,i,k,l) +
           Cap3Volume (vertexList,k,i,j,l) + Cap3Volume (vertexList,l,i,j,k);
}

double Metric::Ball4Area(std::vector<Vertex> &vertexList, int i, int j, int k, int l)
{
    return Cap3Area (vertexList,i,j,k,l) + Cap3Area (vertexList,j,i,k,l) +
           Cap3Area (vertexList,k,i,j,l) + Cap3Area (vertexList,l,i,j,k);
}

double Metric::TetraVolume(double a[], double b[], double c[], double d[])
{
    Matrix A;
    double vol;

    A[0][0] = a[1]; A[0][1] = a[2]; A[0][2] = a[3]; A[0][3] = 1.0;
    A[1][0] = b[1]; A[1][1] = b[2]; A[1][2] = b[3]; A[1][3] = 1.0;
    A[2][0] = c[1]; A[2][1] = c[2]; A[2][2] = c[3]; A[2][3] = 1.0;
    A[3][0] = d[1]; A[3][1] = d[2]; A[3][2] = d[3]; A[3][3] = 1.0;
    vol = fabs((1.0/6.0) * Det4 (A));

    return vol;
}

void Metric::Angle (double ang, Vector3 *T, Vector3 *u, Vector3 *v)
{
    Vector3 Tu, Tv;
    double res;
    Vector3::DiffVector (&Tu,T,u); Tu.Normalize ();
    Vector3::DiffVector (&Tv,T,v); Tv.Normalize ();
    Vector3::DotProduct (&Tu,&Tv,&res);

    ang = fabs (acos (res))/(2.0 * PI);

    if(ang > (1.0/2.0)) ang = 1.0 - ang;
}
