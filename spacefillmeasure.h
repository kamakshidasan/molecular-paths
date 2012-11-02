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

#ifndef SPACEFILLMEASURE_H
#define SPACEFILLMEASURE_H


#include <vector>
#include <cassert>
#include <deluanaycomplex.h>
#include <linktrig.h>
#include <vertex.h>
#include <vector3.h>

class SpaceFillMeasure
{
    private:
        double PI;

        void AngleDihed(double a[], double b[], double c[], double d[], double *ang, double *cos_ang);
        void SegmentCover(double a[], double p1[], double p2[], double d[], double wa, double wd, double *beta);
        void Tetra3Noder(double a[], double b[], double c[], double p[], double rab, double rac, double rbc, double ra, double rb, double rc, double *ang1, double *ang2, double *ang3, double *ang4, double *ang5, double *ang6, double cosine[], double sine[]);
        void TriangleDual(double a[], double b[], double c[], double center[], double *eps2, double ra2, double d1[], double d2[]);

    public:
        SpaceFillMeasure();
        ~SpaceFillMeasure();

        double Det41(double mat4[][5]);
        double TetraVolume(double a[], double b[], double c[], double d[]);
        void PrepareDeriv(DeluanayComplex *delcx, std::vector<int> & nlink_trig, std::vector<LinkTrig> & link_trig);
        void Distance2(std::vector<Vertex> & vertexList, int i, int j, double *distij2);

        void Center2(double a[], double b[], double ra2, double rb2, double rab2, double c[], double *lamda);
        void Center3(double a[], double b[], double c[], double i0, double j0, double k0, double y[]);
        void Center4(double a[], double b[], double c[], double d[], double i0, double j0, double k0, double l0, double y[]);

        void Tetra6Dihedral(double a[], double b[], double c[], double d[], double *ang1, double *ang2, double *ang3, double *ang4, double *ang5, double *ang6);

        void TwoSphereVol(double a[], double b[], double ra, double ra2, double rb, double rb2, double rab, double rab2, double *surfa, double *surfb, double *vola, double *volb, double dsurfa[][3], double dsurfb[][3],double dvola[][3], double dvolb[][3], int option);

        void ThreeVolDist(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double *surfa, double *surfb, double *surfc, double *vola, double *volb, double *volc,double dsurfa[][4], double dsurfb[][4], double dsurfc[][4], double dvola[][4], double dvolb[][4], double dvolc[][4],double dvolda[][4], double dvoldb[][4], double dvoldc[][4], double *eps, double pabc[], double pacb[], double angles[],double *sh_abc, double *sh_acb, double *sh_bca, int option);
        void ThreeVolDir(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double dvola[][4], double dvolb[][4], double dvolc[][4], double angles[], double pabc[], double pacb[],double *eps, double *sh_abc, double *sh_acb, double *sh_bca, int option);
        void ThreeSurfDir(double a[], double b[], double c[], double d1[], double d2[], int nlink, double pabc[], double pacb[], double eps, double ra, double rb, double rc,double ra2, double wa, double wb, double wc, double wd1, double wd2, double rab, double rac, double rbc, int flag_ab, int flag_ac, int flag_bc, double dsurfa[][4], double dsurfb[][4], double dsurfc[][4]);
        void ThreeSphereVol(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double *surfa, double *surfb, double *surfc, double *vola, double *volb, double *volc);

        void FourSphereVol(int a_index,int b_index,int c_index,int d_index,double a[], double b[], double c[], double d[], double ra, double rb, double rc, double rd, double ra2, double rb2, double rc2, double rd2, double rab, double rac, double rad, double rbc, double rbd, double rcd, double rab2, double rac2, double rad2, double rbc2, double rbd2, double rcd2, double wa, double wb, double wc, double wd, double eps1, double eps3, double eps5, double eps7, double shabc, double shacb, double shbca, double shabd, double shadb, double shbda, double shacd, double shadc, double shcda, double shbcd, double shbdc, double shcdb, double pacb[], double pabd[], double padc[], double pbcd[], double ang_abc[], double ang_abd[], double ang_acd[], double ang_bcd[], double *surfa, double *surfb, double *surfc, double *surfd, double *vola, double *volb, double *volc, double *vold, double dsurfa[][5], double dsurfb[][5], double dsurfc[][5], double dsurfd[][5], double dvola[][5], double dvolb[][5], double dvolc[][5], double dvold[][5], int option);
};

#endif // SPACEFILLMEASURE_H
