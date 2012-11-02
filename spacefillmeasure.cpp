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

#include "spacefillmeasure.h"

SpaceFillMeasure::SpaceFillMeasure()
{
        PI = 4 *atan(1.0);
}

SpaceFillMeasure::~SpaceFillMeasure()
{
}

/*!
    \fn SpaceFillMeasure::PrepareDeriv(DeluanayComplex *delcx, std::vector<int> & nlink_trig, std::vector<LinkTrig> & link_trig)
 */
void SpaceFillMeasure::PrepareDeriv(DeluanayComplex *delcx, std::vector<int> & nlink_trig, std::vector<LinkTrig> & link_trig)
{
        int i, j, k, l;
        uint idx;
        int trig1, trig2, trig3, trig4;

        for (idx = 1; idx < delcx->DeluanayTrigs.size(); idx++)
        {
                nlink_trig[idx] = 0;
                link_trig[idx].l1 = 0;
                link_trig[idx].l2 = 0;
        }

        for (idx = 1; idx < delcx->DeluanayTet.size(); idx++)
        {
                if (delcx->DeluanayTet[idx].AlphaStatus != 1) continue;

                i = delcx->DeluanayTet[idx].Corners[1];
                j = delcx->DeluanayTet[idx].Corners[2];
                k = delcx->DeluanayTet[idx].Corners[3];
                l = delcx->DeluanayTet[idx].Corners[4];

                trig1 = delcx->DeluanayTet[idx].TetLink[1];
                trig2 = delcx->DeluanayTet[idx].TetLink[2];
                trig3 = delcx->DeluanayTet[idx].TetLink[3];
                trig4 = delcx->DeluanayTet[idx].TetLink[4];

                nlink_trig[trig1]++;
                if(nlink_trig[trig1] == 1)
                {
                        link_trig[trig1].l1 = i;
                }
                else if(nlink_trig[trig1] == 2)
                {
                        link_trig[trig1].l2 = i;
                }

                nlink_trig[trig2]++;
                if(nlink_trig[trig2] == 1)
                {
                        link_trig[trig2].l1 = j;
                }
                else if(nlink_trig[trig2] == 2)
                {
                        link_trig[trig2].l2 = j;
                }

                nlink_trig[trig3]++;
                if(nlink_trig[trig3] == 1)
                {
                        link_trig[trig3].l1 = k;
                }
                else if(nlink_trig[trig3] == 2)
                {
                        link_trig[trig3].l2 = k;
                }

                nlink_trig[trig4]++;
                if(nlink_trig[trig4] == 1)
                {
                        link_trig[trig4].l1 = l;
                }
                else if(nlink_trig[trig4] == 2)
                {
                        link_trig[trig4].l2 = l;
                }
        }
}

/*!
    \fn SpaceFillMeasure::Distance2(std::vector<Vertex> & vertexList, int i, int j, double *distij2)
 */
void SpaceFillMeasure::Distance2(std::vector<Vertex> & vertexList, int i, int j, double *distij2)
{
        double val = 0.0;

        *distij2 = 0.0;
        for (int k = 1; k <= 3; k++)
        {
                val = vertexList[i].Coordinates[k] - vertexList[j].Coordinates[k];
                *distij2 += val * val;
        }
}

/*!
    \fn SpaceFillMeasure::Det41(double mat4[][5])
 */
double SpaceFillMeasure::Det41(double mat4[][5])
{
        double det4_1, a, b, c, d;
        double val1, val2, val3, val4, val5, val6;

        val1 = mat4[3][3] - mat4[4][3];
        val2 = mat4[2][3] - mat4[4][3];
        val3 = mat4[2][3] - mat4[3][3];
        val4 = mat4[1][3] - mat4[4][3];
        val5 = mat4[1][3] - mat4[3][3];
        val6 = mat4[1][3] - mat4[2][3];

        a = mat4[1][1] * (mat4[2][2] * val1 - mat4[3][2] * val2 + mat4[4][2] * val3);

        b = mat4[2][1] * (mat4[1][2] * val1 - mat4[3][2] * val4 + mat4[4][2] * val5);

        c = mat4[3][1] * (mat4[1][2] * val2 - mat4[2][2] * val4 + mat4[4][2] * val6);

        d = mat4[4][1] * (mat4[1][2] * val3 - mat4[2][2] * val5 + mat4[3][2] * val6);

        det4_1 = a - b + c - d;

        return (det4_1);
}

/*!
    \fn SpaceFillMeasure::Center2(double a[], double b[], double ra2, double rb2, double rab2, double c[], double *lamda)
 */
void SpaceFillMeasure::Center2(double a[], double b[], double ra2, double rb2, double rab2, double c[], double *lamda)
{
        double uml;
        int i;

        *lamda = 0.5 - (ra2 - rb2) / (2 * rab2);
        uml = 1.0 - *lamda;

        for (i = 1; i <= 3; i++)
        {
                c[i] = *lamda * a[i] + uml * b[i];
        }
}

/*!
    \fn SpaceFillMeasure::Center3(double a[], double b[], double c[], double i0, double j0, double k0, double y[])
 */
void SpaceFillMeasure::Center3(double a[], double b[], double c[], double i0, double j0, double k0, double y[])
{
        double a1, a2, a3, a4;
        double dx, dy, dz, d0;
        double val1, val2, val3;

        a1 = (b[2] - a[2]) * (c[3] - a[3]) - (c[2] - a[2]) * (b[3] - a[3]);
        a2 = (b[3] - a[3]) * (c[1] - a[1]) - (c[3] - a[3]) * (b[1] - a[1]);

        val1 = b[1] * c[2] - c[1] * b[2];
        val2 = a[1] * c[2] - c[1] * a[2];
        val3 = a[1] * b[2] - b[1] * a[2];

        a3 = val1 - val2 + val3;
        a4 = a[3] * val1 - b[3] * val2 + c[3] * val3;

        d0 = -a1 * a1 - a2 * a2 - a3 * a3;

        val1 = i0 * (b[3] - c[3]) - j0 * (a[3] - c[3]) + k0 * (a[3] - b[3]);
        val2 = i0 * (b[2] - c[2]) - j0 * (a[2] - c[2]) + k0 * (a[2] - b[2]);
        val3 = i0 * (b[1] - c[1]) - j0 * (a[1] - c[1]) + k0 * (a[1] - b[1]);

        dx = -a4 * a1 + a2 * val1 - a3 * val2;
        dy = -a1 * val1 - a4 * a2 + a3 * val3;
        dz = a1 * val2 - a2 * val3 - a4 * a3;

        y[1] = dx / d0;
        y[2] = dy / d0;
        y[3] = dz / d0;
}

/*!
    \fn SpaceFillMeasure::Center4(double a[], double b[], double c[], double d[], double i0, double j0, double k0, double l0, double y[])
 */
void SpaceFillMeasure::Center4(double a[], double b[], double c[], double d[], double i0, double j0, double k0, double l0, double y[])
{
        int i, j, im, im1;

        double detval[5];
        double mat[5][5];
        double mat4[5][5];

        for (i = 1; i <= 3; i++)
        {
                mat[1][i] = a[i];
                mat[2][i] = b[i];
                mat[3][i] = c[i];
                mat[4][i] = d[i];
        }

        for (i = 1; i <= 4; i++)
        {
                mat4[i][4] = 1;
        }

        for (im = 1; im <= 4; im++)
        {
                for (i = 1; i <= 3; i++)
                {
                        for (j = 1; j <= 4; j++)
                        {
                                mat4[j][i] = mat[j][i];
                        }
                }

                if (im > 1)
                {
                        im1 = im - 1;
                        mat4[1][im1] = i0;
                        mat4[2][im1] = j0;
                        mat4[3][im1] = k0;
                        mat4[4][im1] = l0;
                }
                detval[im] = Det41(mat4);
        }
        y[1] = detval[2] / detval[1];
        y[2] = detval[3] / detval[1];
        y[3] = detval[4] / detval[1];

        FILE *fp = fopen("center4","w");
        fprintf(fp,"%lf %lf %lf",y[1],y[2],y[3]);
        fclose (fp);

}

/*!
    \fn SpaceFillMeasure::Tetra3Noder(double a[], double b[], double c[], double p[], double rab, double rac, double rbc, double ra, double rb, double rc, double *ang1, double *ang2, double *ang3, double *ang4, double *ang5, double *ang6, double cosine[], double sine[])
 */
void SpaceFillMeasure::Tetra3Noder(double a[], double b[], double c[], double p[], double rab, double rac, double rbc, double ra, double rb, double rc, double *ang1, double *ang2, double *ang3, double *ang4, double *ang5, double *ang6, double cosine[], double sine[])
{
        double cos_ang1, cos_ang2, cos_ang3;
        double cos_ang4, cos_ang5 = 0.0, cos_ang6 = 0.0;
        double sin_ang1, sin_ang2, sin_ang3;
        double vol;
        double sum;
        double s1_2, s1;
        double s2_2, s2;
        double s3_2, s3;
        double s4_2, s4, sum_s_2;
/*
        Volume of the tetrahedron, based on all edge lengths (in fact, this is the square of the volume, multiplied by 288)
        The derivative are really (derivative/vol) and are computed from the fact that

        vol2 = 2*ra2*(4*rb2*rc2 - val_bc*val_bc) +
          val_ab*(-2*rc2*val_ab - val_bc*val_ac) -
          val_ac*(val_ab*val_bc + 2*rb2*val_ac)

        where vol2 = 288*vol^2
*/
        FILE *fp = fopen("tetnoder","w");
        vol = TetraVolume(a, b, c, p);
        fprintf(fp,"vol = %f\n",vol);
/*
        Surfaces s1,s2,s3,s4 of the four faces of the tetrahedron,
        (We use the fact that for a triangle T with side lengths a,b,c, then

        P =perimeter = (a+b+c)/2
        Surf^2 = p*(p-a)*(p-b)*(p-c)

        The four triangles considered are:

        Triangle	Surface

        T1 : BCP	s1
        T2 : ACP	s2
        T3 : ABP	s3
        T4 : ABC	s4
*/
        fprintf(fp,"ra = %f\n",ra);
        fprintf(fp,"rb = %f\n",rb);
        fprintf(fp,"rc = %f\n",rc);
        fprintf(fp,"rab = %f\n",rab);
        fprintf(fp,"rbc = %f\n",rbc);
        fprintf(fp,"rac = %f\n",rac);
        sum = (rb + rc + rbc) / 2.0;
        s1_2 = sum * (sum - rb) * (sum - rc) * (sum - rbc);
        s1 = sqrt(s1_2);

        sum = (ra + rc + rac) / 2.0;
        s2_2 = sum * (sum - ra) * (sum - rc) * (sum - rac);
        s2 = sqrt(s2_2);

        sum = (ra + rb + rab) / 2.0;
        s3_2 = sum * (sum - ra) * (sum - rb) * (sum - rab);
        s3 = sqrt(s3_2);

        sum = (rab + rac + rbc) / 2.0;
        s4_2 = sum * (sum - rab) * (sum - rac) * (sum - rbc);
        s4 = sqrt(s4_2);

        sum_s_2 = s1_2 + s2_2 + s3_2 + s4_2;
        fprintf(fp,"sum_s_2 = %f\n",sum_s_2);
/*
        Get all six dihedral angles

        ang1 = angle_dihed(a,b,c,p) = angle_dihed(T3,T4)
        ang2 = angle_dihed(a,c,b,p) = angle_dihed(T2,T4)
        ang3 = angle_dihed(b,c,a,p) = angle_dihed(T1,T4)
        ang4 = angle_dihed(a,p,b,c) = angle_dihed(T2,T3)
        ang5 = angle_dihed(b,p,a,c) = angle_dihed(T1,T3)
        ang6 = angle_dihed(c,p,a,b) = angle_dihed(T1,T2)
*/
        AngleDihed(a, b, c, p, ang1, &cos_ang1);
        AngleDihed(a, c, b, p, ang2, &cos_ang2);
        AngleDihed(b, c, a, p, ang3, &cos_ang3);
        AngleDihed(a, p, b, c, ang4, &cos_ang4);
        AngleDihed(b, p, a, c, ang5, &cos_ang5);
        AngleDihed(c, p, a, b, ang6, &cos_ang6);

        fprintf(fp,"ang1 = %f\n",*ang1);
        fprintf(fp,"ang2 = %f\n",*ang2);
        fprintf(fp,"ang3 = %f\n",*ang3);
        fprintf(fp,"ang4 = %f\n",*ang4);
        fprintf(fp,"ang5 = %f\n",*ang5);
        fprintf(fp,"ang6 = %f\n",*ang6);


        //Get 4 other dihedral angle using cosine rule ----> doesnt work

        /*val1 = 2.0 * s1 * s3 * cos_ang5;
        val2 = 2.0 * s1 * s2 * cos_ang6;

        c1 = sum_s_2 - 2 * s1_2;
        c2 = sum_s_2 - 2 * s2_2 - val1;
        c3 = sum_s_2 - 2 * s3_2 - val2;
        c4 = sum_s_2 - 2 * s4_2 - val1 - val2;

        cos_ang1 = (c1 + c2 - c3 - c4) / (4.0 * s3 * s4);
        cos_ang2 = (c1 - c2 + c3 - c4) / (4.0 * s2 * s4);
        cos_ang3 = (-c1 + c2 + c3 + c4) / (4.0 * s1 * s4);
        cos_ang4 = c4 / (2.0 * s2 * s3);

        *ang1 = acos(cos_ang1) / (2.0 * PI);
        *ang2 = acos(cos_ang2) / (2.0 * PI);
        *ang3 = acos(cos_ang3) / (2.0 * PI);
        *ang4 = acos(cos_ang4) / (2.0 * PI);

        fprintf(fp,"ang1 = %f\n",*ang1);
        fprintf(fp,"ang2 = %f\n",*ang2);
        fprintf(fp,"ang3 = %f\n",*ang3);
        fprintf(fp,"ang4 = %f\n",*ang4);*/

        sin_ang1 = 1.5 * rab * vol / (s3 * s4);
        sin_ang2 = 1.5 * rac * vol / (s2 * s4);
        sin_ang3 = 1.5 * rbc * vol / (s1 * s4);

        cosine[1] = cos_ang1;
        cosine[2] = cos_ang2;
        cosine[3] = cos_ang3;
        sine[1] = sin_ang1;
        sine[2] = sin_ang2;
        sine[3] = sin_ang3;

        fclose(fp);
}


/*!
    \fn SpaceFillMeasure::TetraVolume(double a[], double b[], double c[], double d[])
 */
double SpaceFillMeasure::TetraVolume(double a[], double b[], double c[], double d[])
{
        int i;
        double det, tetra_volume;
        double mat4[5][5];

        for (i = 1; i <= 4; i++)
        {
                mat4[i][4] = 1.0;
        }

        for (i = 1; i <= 3; i++)
        {
                mat4[1][i] = a[i];
                mat4[2][i] = b[i];
                mat4[3][i] = c[i];
                mat4[4][i] = d[i];
        }
        det = Det41(mat4);

        tetra_volume = fabs(det) / 6.0;
        return (tetra_volume);
}

/*!
    \fn SpaceFillMeasure::AngleDihed(double a[], double b[], double c[], double d[], double *ang, double *cos_ang)
 */
void SpaceFillMeasure::AngleDihed(double a[], double b[], double c[], double d[], double *ang, double *cos_ang)
{
        Vector3 *va = new Vector3(a[1], a[2], a[3]);
        Vector3 *vb = new Vector3(b[1], b[2], b[3]);
        Vector3 *vc = new Vector3(c[1], c[2], c[3]);
        Vector3 *vd = new Vector3(d[1], d[2], d[3]);

        Vector3 *u1 = new Vector3();
        Vector3 *u2 = new Vector3();
        Vector3 *m = new Vector3();
        Vector3 *n1 = new Vector3();
        Vector3 *n2 = new Vector3();

        Vector3::DiffVector(u1,va,vc);
        Vector3::DiffVector(u2,vb,vc);

        Vector3::CrossProduct(m,u1,u2);
        Vector3::Normalize(n1,m);

        Vector3::DiffVector(u1,va,vd);
        Vector3::DiffVector(u2,vb,vd);


        Vector3::CrossProduct(m,u1,u2);
        Vector3::Normalize(n2,m);

        Vector3::DotProduct(n1, n2,cos_ang);
        if(*cos_ang > 1.0)
        {
                *cos_ang = 1.0;
        }
        else if(*cos_ang < -1.0)
        {
                *cos_ang = -1.0;
        }
        *ang = acos(*cos_ang) / (2.0 * PI);

        delete va;delete vb;delete vc;delete vd;
        delete u1;delete u2;delete m;delete n1;delete n2;

}


/*!
    \fn SpaceFillMeasure::TriangleDual(double a[], double b[], double c[], double center[], double *eps2, double ra2, double d1[], double d2[])
 */
/*void SpaceFillMeasure::TriangleDual(double a[], double b[], double c[], double center[], double *eps2, double ra2, double d1[], double d2[])
{
        int i;
        double dn[4];
        double s2, s3, s1,eps,scale = 1e3;

        mpz_t a_mp[3],b_mp[3],c_mp[3],cen_mp[3];
        mpz_t u1[3],u2[3],u3[3],y[3],n[3],N[3];
        mpz_t temp1,temp2,temp3,s,S1,S2,S3,r2;
        mpz_init(temp1);mpz_init(temp2);mpz_init(temp3);mpz_init(s);
        mpz_init(S1);mpz_init(S2);mpz_init(S3);mpz_init(r2);

        for(int j = 0;j<3;j++)
        {
                mpz_init(a_mp[j]);mpz_set_d(a_mp[j],a[j+1]*scale);
                mpz_init(b_mp[j]);mpz_set_d(b_mp[j],b[j+1]*scale);
                mpz_init(c_mp[j]);mpz_set_d(c_mp[j],c[j+1]*scale);
                mpz_init(cen_mp[j]);mpz_set_d(cen_mp[j],center[j+1]*scale);
                mpz_init(u1[j]);mpz_sub(u1[j],b_mp[j],a_mp[j]);
                mpz_init(u2[j]);mpz_sub(u2[j],c_mp[j],a_mp[j]);
                mpz_init(u3[j]);mpz_sub(u3[j],cen_mp[j],a_mp[j]);
                mpz_init(n[j]);
                mpz_init(N[j]);
        }
        mpz_set_d(r2,ra2*scale*scale);

        mpz_mul(temp1,a_mp[1],b_mp[2]);
        mpz_mul(temp2,a_mp[2],b_mp[1]);
        mpz_sub(n[0],temp1,temp2);

        mpz_mul(temp1,a_mp[2],b_mp[0]);
        mpz_mul(temp2,a_mp[0],b_mp[2]);
        mpz_sub(n[1],temp1,temp2);

        mpz_mul(temp1,a_mp[0],b_mp[1]);
        mpz_mul(temp2,a_mp[1],b_mp[0]);
        mpz_sub(n[2],temp1,temp2);

        mpz_mul(temp3,n[0],n[0]);
        mpz_mul(temp2,n[1],n[1]);
        mpz_mul(temp1,n[2],n[2]);

        mpz_add(temp3,temp3,temp2);
        mpz_add(temp3,temp3,temp1);

        mpz_set(S2,temp3);

        double nval = mpz_get_d(S2);
        double n_1 = mpz_get_d(n[0])/sqrt(nval);
        double n_2 = mpz_get_d(n[1])/sqrt(nval);
        double n_3 = mpz_get_d(n[2])/sqrt(nval);

        mpz_set_d(N[0],n_1*scale);
        mpz_set_d(N[1],n_2*scale);
        mpz_set_d(N[2],n_3*scale);

        //s2 = mpz_get_d(temp3)/1e32;

        mpz_mul(temp1,u3[0],N[0]);
        mpz_mul(temp2,u3[1],N[1]);
        mpz_mul(temp3,u3[2],N[2]);

        mpz_add(temp3,temp3,temp2);
        mpz_add(temp3,temp3,temp1);

        mpz_set(S1,temp3);

        //s1 = mpz_get_d(temp3)/1e32;

        mpz_mul(temp1,u3[0],u3[0]);
        mpz_mul(temp2,u3[1],u3[1]);
        mpz_mul(temp3,u3[2],u3[2]);

        mpz_add(temp3,temp3,temp2);
        mpz_add(temp3,temp3,temp1);

        //s3 = mpz_get_d(temp3)/1e16;
        mpz_set(S3,temp3);

        mpz_mul(temp1,S1,S1);
        mpz_mul(temp2,S2,S3);
        mpz_mul(temp3,S2,r2);

        mpz_add(temp1,temp1,temp3);
        mpz_sub(temp1,temp1,temp2);

        mpz_sub(temp2,r2,S3);

        double value1 = mpz_get_d(temp2);
        double value2 = mpz_get_d(temp1);

        mpz_clear(temp1);mpz_clear(temp2);mpz_clear(temp3);mpz_clear(s);
        mpz_clear(S1);mpz_clear(S2);mpz_clear(S3);mpz_clear(r2);

// 	Vector3::DotProduct(N,N,&s2);
// 	Vector3::DotProduct(u3,N,&s1);
// 	Vector3::DotProduct(u3,u3,&s3);

        //(*eps2) = sqrt(ra2 - s3);
        //eps = *eps2 / sqrt(s2);

        eps = (-s1 + sqrt(s1 * s1 - s3 * s2 + ra2 * s2)) / s2;

// 	dn[1] = N->X;
// 	dn[2] = N->Y;
// 	dn[3] = N->Z;
//
// 	for (i = 1; i <= 3; i++)
// 	{
// 		d1[i] = center[i] + eps * dn[i];
// 		d2[i] = center[i] - eps * dn[i];
// 	}
//
// 	N->X = -N->X;
// 	N->Y = -N->Y;
// 	N->Z = -N->Z;
//
// 	Vector3::DotProduct(N,N,&s2);
// 	Vector3::DotProduct(u3,N,&s1);
// 	Vector3::DotProduct(u3,u3,&s3);
//
// 	eps = (-s1 + sqrt(s1 * s1 - s3 * s2 + ra2 * s2)) / s2;
//
// 	dn[1] = N->X;
// 	dn[2] = N->Y;
// 	dn[3] = N->Z;
//
// 	for (i = 1; i <= 3; i++)
// 	{
// 		d2[i] = center[i] + eps * dn[i];
// 		//d2[i] = center[i] - eps * dn[i];
// 	}
//
// 	delete s;delete t;delete u;delete y;delete n;delete u1;delete u2,delete u3;delete N;
}*/

void SpaceFillMeasure::TriangleDual(double a[], double b[], double c[], double center[], double *eps2, double ra2, double d1[], double d2[])
{
        int i;
        double dn[4];
        double s2, s3, s1,eps;

        Vector3 *s = new Vector3(a[1],a[2],a[3]);
        Vector3 *t = new Vector3(b[1],b[2],b[3]);
        Vector3 *u = new Vector3(c[1],c[2],c[3]);
        Vector3 *y = new Vector3(center[1],center[2],center[3]);

        Vector3 *u1 = new Vector3();
        Vector3 *u2 = new Vector3();
        Vector3 *u3 = new Vector3();

        Vector3::DiffVector(u1,t,s);
        Vector3::DiffVector(u2,u,s);
        Vector3::DiffVector(u3,y,s);

        Vector3 *n = new Vector3();
        Vector3 *N = new Vector3();

        Vector3::CrossProduct(n,u1,u2);
        Vector3::Normalize(N,n);

        Vector3::DotProduct(N,N,&s2);
        Vector3::DotProduct(u3,N,&s1);
        Vector3::DotProduct(u3,u3,&s3);


        if(s3>ra2)
        {
                eps = 0.0;
                *eps2 = 0.0;
        }
        else
        {
                //eps = (-s1 + sqrt(s1 * s1 - s3 * s2 + ra2 * s2)) / s2;
                //*eps2 = eps * eps;
            (*eps2) = sqrt(ra2 - s3);
            eps = *eps2 / sqrt(s2);
        }
        dn[1] = N->X;
        dn[2] = N->Y;
        dn[3] = N->Z;

        for (i = 1; i <= 3; i++)
        {
                d1[i] = center[i] + eps * dn[i];
                d2[i] = center[i] - eps * dn[i];
        }

        /*N->X = -N->X;
        N->Y = -N->Y;
        N->Z = -N->Z;

        Vector3::DotProduct(N,N,&s2);
        Vector3::DotProduct(u3,N,&s1);
        Vector3::DotProduct(u3,u3,&s3);
        if(s3>ra2)
        {
                eps = 0.0;
                *eps2 = 0.0;
        }
        else
        {
                eps = (-s1 + sqrt(s1 * s1 - s3 * s2 + ra2 * s2)) / s2;
                *eps2 = eps * eps;
        }
        dn[1] = N->X;
        dn[2] = N->Y;
        dn[3] = N->Z;

        for (i = 1; i <= 3; i++)
        {
                d2[i] = center[i] + eps * dn[i];
                //d2[i] = center[i] - eps * dn[i];
        }*/

        delete s;delete t;delete u;delete y;delete n;delete u1;delete u2,delete u3;delete N;
}

/*!
    \fn SpaceFillMeasure::SegmentCover(double a[], double p1[], double p2[], double d[], double wa, double wd, double *beta)
 */
void SpaceFillMeasure::SegmentCover(double a[], double p1[], double p2[], double d[], double wa, double wd, double *beta)
{
        double val1, val2, test;

        Vector3 *ad = new Vector3();
        Vector3 *p1p2 = new Vector3();

        Vector3 *op1 = new Vector3(p1[1], p1[2], p1[3]);
        Vector3 *op2 = new Vector3(p2[1], p2[2], p2[3]);

        Vector3::DiffVector(p1p2,op2,op1);

        delete op1;delete op2;

        op1 = new Vector3(a[1], a[2], a[3]);
        op2 = new Vector3(d[1], d[2], d[3]);

        Vector3::DiffVector(ad,op2,op1);

        delete op2;

        op2 = new Vector3(p2[1], p2[2], p2[3]);

        Vector3::DotProduct(ad, p1p2,&val1);
        Vector3::DotProduct(ad, op2,&val2);

        assert (val1 != 0);

        *beta = (wa - wd + val2) / val1;

        delete op1;delete op2;

        op1 = new Vector3(p1[1], p1[2], p1[3]);
        op2 = new Vector3(d[1], d[2], d[3]);

        Vector3::DotProduct(op1, op1, &val1);
        Vector3::DotProduct(op1, op2, &val2);

        test = val1 - 2.0 * (val2 - wd);

        if (test < 0)
        {
                *beta = 1.0 - *beta;
        }

        delete op1;delete op2;delete ad;delete p1p2;
}

/*!
    \fn SpaceFillMeasure::Tetra6Dihedral(double a[], double b[], double c[], double d[], double *ang1, double *ang2, double *ang3, double *ang4, double *ang5, double *ang6)
 */
void SpaceFillMeasure::Tetra6Dihedral(double a[], double b[], double c[], double d[], double *ang1, double *ang2, double *ang3, double *ang4, double *ang5, double *ang6)
{
        double cos_ang;

        Vector3 *va = new Vector3(a[1], a[2], a[3]);
        Vector3 *vb = new Vector3(b[1], b[2], b[3]);
        Vector3 *vc = new Vector3(c[1], c[2], c[3]);
        Vector3 *vd = new Vector3(d[1], d[2], d[3]);

        Vector3 *u_ab = new Vector3();
        Vector3 *u_ac = new Vector3();
        Vector3 *u_ad = new Vector3();
        Vector3 *u_bc = new Vector3();
        Vector3 *u_bd = new Vector3();
        Vector3 *u_abc = new Vector3();
        Vector3 *u_abd = new Vector3();
        Vector3 *u_acd = new Vector3();
        Vector3 *u_bcd = new Vector3();
        Vector3 *n_abc = new Vector3();
        Vector3 *n_abd = new Vector3();
        Vector3 *n_acd = new Vector3();
        Vector3 *n_bcd = new Vector3();

        Vector3::DiffVector(u_ab,vb, va);
        Vector3::DiffVector(u_ad,vd, va);
        Vector3::DiffVector(u_ac,vc, va);
        Vector3::DiffVector(u_bc,vc, vb);
        Vector3::DiffVector(u_bd,vd, vb);

        Vector3::CrossProduct(u_abc,u_ab, u_ac);
        Vector3::CrossProduct(u_abd,u_ab, u_ad);
        Vector3::CrossProduct(u_acd,u_ac, u_ad);
        Vector3::CrossProduct(u_bcd,u_bc, u_bd);

        Vector3::Normalize(n_abc,u_abc);
        Vector3::Normalize(n_abd,u_abd);
        Vector3::Normalize(n_acd,u_acd);
        Vector3::Normalize(n_bcd,u_bcd);

        Vector3::DotProduct(n_abc, n_abd,&cos_ang);
        if(cos_ang > 1.0)
        {
                cos_ang = 1.0;
        }
        else if(cos_ang < -1.0)
        {
                cos_ang = -1.0;
        }
        *ang1 = acos(cos_ang) / (2 * PI);

        Vector3::DotProduct(n_abc, n_acd,&cos_ang);
        if(cos_ang > 1.0)
        {
                cos_ang = 1.0;
        }
        else if(cos_ang < -1.0)
        {
                cos_ang = -1.0;
        }
        *ang2 = acos(-cos_ang) / (2 * PI);

        Vector3::DotProduct(n_abc, n_bcd,&cos_ang);
        if(cos_ang > 1.0)
        {
                cos_ang = 1.0;
        }
        else if(cos_ang < -1.0)
        {
                cos_ang = -1.0;
        }
        *ang3 = acos(cos_ang) / (2 * PI);

        Vector3::DotProduct(n_abd, n_acd,&cos_ang);
        if(cos_ang > 1.0)
        {
                cos_ang = 1.0;
        }
        else if(cos_ang < -1.0)
        {
                cos_ang = -1.0;
        }
        *ang4 = acos(cos_ang) / (2 * PI);

        Vector3::DotProduct(n_abd, n_bcd,&cos_ang);
        if(cos_ang > 1.0)
        {
                cos_ang = 1.0;
        }
        else if(cos_ang < -1.0)
        {
                cos_ang = -1.0;
        }
        *ang5 = acos(-cos_ang) / (2 * PI);

        Vector3::DotProduct(n_acd, n_bcd,&cos_ang);
        if(cos_ang > 1.0)
        {
                cos_ang = 1.0;
        }
        else if(cos_ang < -1.0)
        {
                cos_ang = -1.0;
        }
        *ang6 = acos(cos_ang) / (2 * PI);
}

/*!
    \fn SpaceFillMeasure::TwoSphereVol(double a[], double b[], double ra, double ra2, double rb, double rb2, double rab, double rab2, double *surfa, double *surfb, double *vola, double *volb, double dsurfa[][3], double dsurfb[][3],double dvola[][3], double dvolb[][3], int option)
 */
void SpaceFillMeasure::TwoSphereVol(double a[], double b[], double ra, double ra2, double rb, double rb2, double rab, double rab2, double *surfa, double *surfb, double *vola, double *volb, double dsurfa[][3], double dsurfb[][3],double dvola[][3], double dvolb[][3], int option)
{
        double vala, valb, lamda = 0.0, ha, hb, sa, ca, sb, cb;
        double dera, derb, der1, der2, coef1, coef2, coefa, coefb, Aab;

        double c[4];
        double u_ab[4];

        //Get "center" of the two spheres
        Center2(a, b, ra2, rb2, rab2, c, &lamda);

        valb = lamda * rab;
        vala = rab - valb;

        //Get height of the cap of sphere A occluded by sphere B
        ha = ra - vala;

        //same for sphere B
        hb = rb - valb;

        //Get surfaces of intersection
        *surfa = 2.0 * PI * ra * ha;
        *surfb = 2.0 * PI * rb * hb;

        //get volume

        Aab = PI * (ra2 - vala * vala);

        sa = ra * (*surfa);
        ca = vala * Aab;

        *vola = (sa - ca) / 3;

        sb = rb * (*surfb);
        cb = valb * Aab;

        *volb = (sb - cb) / 3;

        //Compute derivatives, if needed

        if (option == 0) return;

        for (int i = 1; i <= 3; i++)
        {
                u_ab[i] = (a[i] - b[i]) / rab;
        }

        dera = -lamda;
        derb = 1.0 - lamda;

        der1 = 2.0 * PI * ra * dera;
        der2 = 2.0 * PI * rb * derb;

        for(int i=1;i<=3;i++)
        {
                coef1 = der1 * u_ab[i];
                coef2 = der2 * u_ab[i];
                dsurfa[i][1] = coef1;
                dsurfa[i][2] = -coef1;
                dsurfb[i][1] = coef2;
                dsurfb[i][2] = -coef2;
        }

        coefa = Aab * lamda;
        coefb = Aab - coefa;

        for(int i=1;i<=3;i++)
        {
                coef1 = -coefa * u_ab[i];
                coef2 = -coefb * u_ab[i];
                dvola[i][1] = coef1;
                dvola[i][2] = -coef1;
                dvolb[i][1] = coef2;
                dvolb[i][2] = -coef2;
        }
}

/*!
    \fn SpaceFillMeasure::ThreeVolDist(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double *surfa, double *surfb, double *surfc, double *vola, double *volb, double *volc,double dsurfa[][4], double dsurfb[][4], double dsurfc[][4], double dvola[][4], double dvolb[][4], double dvolc[][4],double dvolda[][4], double dvoldb[][4], double dvoldc[][4], double *eps, double pabc[], double pacb[], double angles[],double *sh_abc, double *sh_acb, double *sh_bca, int option)
 */
void SpaceFillMeasure::ThreeVolDist(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double *surfa, double *surfb, double *surfc, double *vola, double *volb, double *volc,double dsurfa[][4], double dsurfb[][4], double dsurfc[][4], double dvola[][4], double dvolb[][4], double dvolc[][4],double dvolda[][4], double dvoldb[][4], double dvoldc[][4], double *eps, double pabc[], double pacb[], double angles[],double *sh_abc, double *sh_acb, double *sh_bca, int option)
{
        int i;

        double a1, a2, a3, s2, c1, c2;
        double seg_ang_abc = 0.0, seg_ang_acb = 0.0, seg_ang_bca = 0.0;
        double ang_abc, ang_acb, ang_bca;
        double cos_abc, cos_acb, cos_bca;
        double sin_abc, sin_acb, sin_bca;
        double ang_dih_abc = 0.0, ang_dih_cab = 0.0, ang_dih_bac = 0.0;
        double s_abc, s_acb, s_bca;
        double val1, val2, val3, l1 = 0.0, l2 = 0.0, l3 = 0.0;
        double val1b, val2b, val3b;
        double coef, coef1, coef2, coef3;
        double rho_ab2, rho_ac2, rho_bc2;
        double der_ab, der_ac, der_bc, diff_ab, diff_ac, diff_bc;
        double dsurfa_ab, dsurfa_ac;
        double dsurfb_ab, dsurfb_bc;
        double dsurfc_ac, dsurfc_bc;
        double coef_ab, coef_ac, coef_bc;
        double center[4];
        double c_ab[4];
        double c_ac[4];
        double c_bc[4];
        double u_ab[4];
        double u_ac[4];
        double u_bc[4];
        double Vab[4];
        double Vac[4];
        double Vbc[4];
        double cosine[4];
        double sine[4];

        Center2(a, b, ra2, rb2, rab2, c_ab, &l1);
        Center2(a, c, ra2, rc2, rac2, c_ac, &l2);
        Center2(b, c, rb2, rc2, rbc2, c_bc, &l3);

        val1 = l1 * rab;
        val2 = l2 * rac;
        val3 = l3 * rbc;

        val1b = rab - val1;
        val2b = rac - val2;
        val3b = rbc - val3;

        Center3(a, b, c, wa, wb, wc, center);

        TriangleDual(a, b, c, center, eps, ra2, pabc, pacb);

        Tetra3Noder(a, b, c, pabc, rab, rac, rbc, ra, rb, rc, &seg_ang_abc, &seg_ang_acb, &seg_ang_bca, &ang_dih_abc, &ang_dih_bac, &ang_dih_cab, cosine, sine);

        angles[1] = seg_ang_abc;
        angles[2] = seg_ang_acb;
        angles[3] = seg_ang_bca;

        FILE *fp = fopen("threevoldist.txt","w");

        fprintf(fp,"a(1) = %lf\n",a[1]);
        fprintf(fp,"a(2) = %lf\n",a[2]);
        fprintf(fp,"a(3) = %lf\n",a[3]);

        fprintf(fp,"b(1) = %lf\n",b[1]);
        fprintf(fp,"b(2) = %lf\n",b[2]);
        fprintf(fp,"b(3) = %lf\n",b[3]);

        fprintf(fp,"c(1) = %lf\n",c[1]);
        fprintf(fp,"c(2) = %lf\n",c[2]);
        fprintf(fp,"c(3) = %lf\n",c[3]);

        fprintf(fp,"center(1) = %lf\n",center[1]);
        fprintf(fp,"center(2) = %lf\n",center[2]);
        fprintf(fp,"center(3) = %lf\n",center[3]);
        fprintf(fp,"angles(1) = %lf\n",angles[1]);
        fprintf(fp,"angles(2) = %lf\n",angles[2]);
        fprintf(fp,"angles(3) = %lf\n",angles[3]);
        fprintf(fp,"ang_dih_abc = %lf\n",ang_dih_abc);
        fprintf(fp,"ang_dih_bac = %lf\n",ang_dih_bac);
        fprintf(fp,"ang_dih_cab = %lf\n",ang_dih_cab);

        fprintf(fp,"cos(1) = %lf\n",cosine[1]);
        fprintf(fp,"cos(2) = %lf\n",cosine[2]);
        fprintf(fp,"cos(3) = %lf\n",cosine[3]);
        fprintf(fp,"sine(1) = %lf\n",sine[1]);
        fprintf(fp,"sine(2) = %lf\n",sine[2]);
        fprintf(fp,"sine(3) = %lf\n",sine[3]);

        fprintf(fp,"eps = %lf\n",*eps);

        fclose(fp);

        a1 = ra * (1.0 - 2.0 * ang_dih_abc);
        a2 = 2.0 * seg_ang_abc * val1b;
        a3 = 2.0 * seg_ang_acb * val2b;

        *surfa = 2.0 * PI * ra * (a1 - a2 - a3);

        a1 = rb * (1.0 - 2.0 * ang_dih_bac);
        a2 = 2.0 * seg_ang_abc * val1;
        a3 = 2.0 * seg_ang_bca * val3b;

        *surfb = 2.0 * PI * rb * (a1 - a2 - a3);

        a1 = rc * (1.0 - 2.0 * ang_dih_cab);
        a2 = 2.0 * seg_ang_acb * val2;
        a3 = 2.0 * seg_ang_bca * val3;

        *surfc = 2.0 * PI * rc * (a1 - a2 - a3);

        ang_abc = 2.0 * PI * seg_ang_abc;
        ang_acb = 2.0 * PI * seg_ang_acb;
        ang_bca = 2.0 * PI * seg_ang_bca;

        cos_abc = cosine[1];
        sin_abc = sine[1];
        cos_acb = cosine[2];
        sin_acb = sine[2];
        cos_bca = cosine[3];
        sin_bca = sine[3];

        rho_ab2 = ra2 - val1b * val1b;
        rho_ac2 = ra2 - val2b * val2b;
        rho_bc2 = rb2 - val3b * val3b;

        s_abc = rho_ab2 * (ang_abc - sin_abc * cos_abc);
        s_acb = rho_ac2 * (ang_acb - sin_acb * cos_acb);
        s_bca = rho_bc2 * (ang_bca - sin_bca * cos_bca);

        s2 = ra * (*surfa);
        c1 = val1b * s_abc;
        c2 = val2b * s_acb;

        *vola = (s2 - c1 - c2) / 3.0;

        s2 = rb * (*surfb);
        c1 = val1 * s_abc;
        c2 = val3b * s_bca;

        *volb = (s2 - c1 - c2) / 3.0;

        s2 = rc * (*surfc);
        c1 = val2 * s_acb;
        c2 = val3 * s_bca;

        *volc = (s2 - c1 - c2) / 3.0;

        *sh_abc = *eps * cos_abc / sin_abc;
        *sh_acb = *eps * cos_acb / sin_acb;
        *sh_bca = *eps * cos_bca / sin_bca;

        if (option == 0) return;

        for (i = 1; i <= 3; i++)
        {
                u_ab[i] = (a[i] - b[i]) / rab;
                u_ac[i] = (a[i] - c[i]) / rac;
                u_bc[i] = (b[i] - c[i]) / rbc;
        }

        dsurfa_ab = -2.0 * ra * ang_abc * l1;
        dsurfa_ac = -2.0 * ra * ang_acb * l2;

        for (i = 1; i <= 3; i++)
        {
                diff_ab = dsurfa_ab * u_ab[i];
                diff_ac = dsurfa_ac * u_ac[i];
                dsurfa[i][1] = diff_ab + diff_ac;
                dsurfa[i][2] = -diff_ab;
                dsurfa[i][3] = -diff_ac;
        }

        dsurfb_ab = -2.0 * rb * ang_abc * (1.0 - l1);
        dsurfb_bc = -2.0 * rb * ang_bca * l3;

        for (i = 1; i <= 3; i++)
        {
                diff_ab = dsurfb_ab * u_ab[i];
                diff_bc = dsurfb_bc * u_bc[i];
                dsurfb[i][1] = diff_ab;
                dsurfb[i][2] = -diff_ab + diff_bc;
                dsurfb[i][3] = -diff_bc;
        }

        dsurfc_ac = -2.0 * rc * ang_acb * (1.0 - l2);
        dsurfc_bc = -2.0 * rc * ang_bca * (1.0 - l3);

        for (i = 1; i <= 3; i++)
        {
                diff_ac = dsurfc_ac * u_ac[i];
                diff_bc = dsurfc_bc * u_bc[i];
                dsurfc[i][1] = diff_ac;
                dsurfc[i][2] = diff_bc;
                dsurfc[i][3] = -diff_ac - diff_bc;
        }

        der_ab = -s_abc * l1;
        der_ac = -s_acb * l2;
        der_bc = 0.0;

        for (i = 1; i <= 3; i++)
        {
                coef_ab = der_ab * u_ab[i];
                coef_ac = der_ac * u_ac[i];
                dvola[i][1] = coef_ab + coef_ac;
                dvola[i][2] = -coef_ab;
                dvola[i][3] = -coef_ac;
        }

        der_ab = s_abc * (1.0 - l1);
        der_ac = 0.0;
        der_bc = -s_bca * l3;

        for (i = 1; i <= 3; i++)
        {
                coef_ab = der_ab * u_ab[i];
                coef_bc = der_bc * u_bc[i];
                dvolb[i][1] = -coef_ab;
                dvolb[i][2] = coef_ab + coef_bc;
                dvolb[i][3] = -coef_bc;
        }

        der_ab = 0.0;
        der_ac = s_acb * (1.0 - l2);
        der_bc = s_bca * (1.0 - l3);

        for (i = 1; i <= 3; i++)
        {
                coef_ac = der_ac * u_ac[i];
                coef_bc = der_bc * u_bc[i];
                dvolc[i][1] = -coef_ac;
                dvolc[i][2] = -coef_bc;
                dvolc[i][3] = coef_ac + coef_bc;
        }

        for (i = 1; i <= 3; i++)
        {
                Vab[i] = center[i] - c_ab[i];
                Vac[i] = center[i] - c_ac[i];
                Vbc[i] = center[i] - c_bc[i];
        }

        coef = -2.0 * pow(*eps,3.0) / 3.0;
        coef1 = coef / (*sh_abc * rab);
        coef2 = coef / (*sh_acb * rac);
        coef3 = coef / (*sh_bca * rbc);

        for (i = 1; i <= 3; i++)
        {
                coef_ab = coef1 * Vab[i];
                coef_ac = coef2 * Vac[i];
                coef_bc = coef3 * Vbc[i];
                dvolda[i][1] = coef_ab + coef_ac;
                dvolda[i][2] = -coef_ab;
                dvolda[i][3] = -coef_ac;
                dvoldb[i][1] = -coef_ab;
                dvoldb[i][2] = coef_ab + coef_bc;
                dvoldb[i][3] = -coef_bc;
                dvoldc[i][1] = -coef_ac;
                dvoldc[i][2] = -coef_bc;
                dvoldc[i][3] = coef_ac + coef_bc;
        }
}

/*!
    \fn SpaceFillMeasure::ThreeVolDir(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double dvola[][4], double dvolb[][4], double dvolc[][4], double angles[], double pabc[], double pacb[],double *eps, double *sh_abc, double *sh_acb, double *sh_bca, int option)
 */
void SpaceFillMeasure::ThreeVolDir(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double dvola[][4], double dvolb[][4], double dvolc[][4], double angles[], double pabc[], double pacb[],double *eps, double *sh_abc, double *sh_acb, double *sh_bca, int option)
{
        int i;

        double seg_ang_abc = 0.0, seg_ang_acb = 0.0, seg_ang_bca = 0.0;
        double cos_abc, cos_acb, cos_bca;
        double sin_abc, sin_acb, sin_bca;
        double ang_dih_abc = 0.0, ang_dih_cab = 0.0, ang_dih_bac = 0.0;
        double l1 = 0.0, l2 = 0.0, l3 = 0.0;
        double coef_ab, coef_ac, coef_bc;
        double coef, coef1, coef2, coef3;

        double center[4];
        double c_ab[4];
        double c_ac[4];
        double c_bc[4];
        double Vab[4];
        double Vac[4];
        double Vbc[4];
        double cosine[4];
        double sine[4];

        Center2(a, b, ra2, rb2, rab2, c_ab, &l1);
        Center2(a, c, ra2, rc2, rac2, c_ac, &l2);
        Center2(b, c, rb2, rc2, rbc2, c_bc, &l3);

        Center3(a, b, c, wa, wb, wc, center);

        TriangleDual(a, b, c, center, eps, ra2, pabc, pacb);

        Tetra3Noder(a, b, c, pabc, rab, rac, rbc, ra, rb, rc, &seg_ang_abc, &seg_ang_acb, &seg_ang_bca, &ang_dih_abc, &ang_dih_bac, &ang_dih_cab, cosine, sine);

        angles[1] = seg_ang_abc;
        angles[2] = seg_ang_acb;
        angles[3] = seg_ang_bca;

        cos_abc = cosine[1];
        sin_abc = sine[1];
        cos_acb = cosine[2];
        sin_acb = sine[2];
        cos_bca = cosine[3];
        sin_bca = sine[3];

        *sh_abc = *eps * cos_abc / sin_abc;
        *sh_acb = *eps * cos_acb / sin_acb;
        *sh_bca = *eps * cos_bca / sin_bca;

        FILE *fp = fopen("threevoldir.txt","w");
        fprintf(fp,"angles(1) = %lf\n",angles[1]);
        fprintf(fp,"angles(2) = %lf\n",angles[2]);
        fprintf(fp,"angles(3) = %lf\n",angles[3]);
        fprintf(fp,"ang_dih_abc = %lf\n",ang_dih_abc);
        fprintf(fp,"ang_dih_bac = %lf\n",ang_dih_bac);
        fprintf(fp,"ang_dih_cab = %lf\n",ang_dih_cab);

        fprintf(fp,"cos(1) = %lf\n",cosine[1]);
        fprintf(fp,"cos(2) = %lf\n",cosine[2]);
        fprintf(fp,"cos(3) = %lf\n",cosine[3]);
        fprintf(fp,"sine(1) = %lf\n",sine[1]);
        fprintf(fp,"sine(2) = %lf\n",sine[2]);
        fprintf(fp,"sine(3) = %lf\n",sine[3]);

        fprintf(fp,"eps = %lf\n",*eps);

        fclose(fp);

        if (option == 0) return;

        for (i = 1; i <= 3; i++)
        {
                Vab[i] = center[i] - c_ab[i];
                Vac[i] = center[i] - c_ac[i];
                Vbc[i] = center[i] - c_bc[i];
        }

        coef = -2.0 * pow(*eps, 3.0) / 3.0;
        coef1 = coef / (*sh_abc * rab);
        coef2 = coef / (*sh_acb * rac);
        coef3 = coef / (*sh_bca * rbc);

        for (i = 1; i <= 3; i++)
        {
                coef_ab = coef1 * Vab[i];
                coef_ac = coef2 * Vac[i];
                coef_bc = coef3 * Vbc[i];
                dvola[i][1] = coef_ab + coef_ac;
                dvola[i][2] = -coef_ab;
                dvola[i][3] = -coef_ac;
                dvolb[i][1] = -coef_ab;
                dvolb[i][2] = coef_ab + coef_bc;
                dvolb[i][3] = -coef_bc;
                dvolc[i][1] = -coef_ac;
                dvolc[i][2] = -coef_bc;
                dvolc[i][3] = coef_ac + coef_bc;
        }
}

/*!
    \fn SpaceFillMeasure::ThreeSurfDir(double a[], double b[], double c[], double d1[], double d2[], int nlink, double pabc[], double pacb[], double eps, double ra, double rb, double rc,double ra2, double wa, double wb, double wc, double wd1, double wd2, double rab, double rac, double rbc, int flag_ab, int flag_ac, int flag_bc, double dsurfa[][4], double dsurfb[][4], double dsurfc[][4])
 */
void SpaceFillMeasure::ThreeSurfDir(double a[], double b[], double c[], double d1[], double d2[], int nlink, double pabc[], double pacb[], double eps, double ra, double rb, double rc,double ra2, double wa, double wb, double wc, double wd1, double wd2, double rab, double rac, double rbc, int flag_ab, int flag_ac, int flag_bc, double dsurfa[][4], double dsurfb[][4], double dsurfc[][4])
{
        int i, j;
        int flag_a, flag_b, flag_c;

        double beta = 0, beta1 = 0, beta2 = 0;
        double coef2, coef3, coef_ab, coef_ac, coef_bc;
        double c_ab_ac, c_ab_bc, c_ac_bc;

        double center[4];
        double u_ab[4];
        double u_ac[4];
        double u_bc[4];
        double e_abc[4];
        double e_cab[4];
        double e_bca[4];
        double u_abc[4];
        double u_cab[4];
        double u_bca[4];

        for (i = 1; i <= 3; i++)
        {
                for (j = 1; j <= 3; j++)
                {
                        dsurfa[j][i] = 0;
                        dsurfb[j][i] = 0;
                        dsurfc[j][i] = 0;
                }
        }

        flag_a = flag_ab | flag_ac;
        flag_b = flag_ab | flag_bc;
        flag_c = flag_ac | flag_bc;

        if (nlink == 2)
        {
                Center3(a, b, c, wa, wb, wc, center);
                TriangleDual(a, b, c, center, &eps, ra2, pabc, pacb);
        }

        if (nlink == 1)
        {
                SegmentCover(a, pabc, pacb, d1, wa, wd1, &beta);
                beta = 1 - beta;
        }
        else if (nlink == 2)
        {
                SegmentCover(a, pabc, pacb, d1, wa, wd1, &beta1);
                SegmentCover(a, pabc, pacb, d2, wa, wd2, &beta2);
                beta = 1 - (beta1 + beta2);
        }
        else
        {
                beta = 1;
        }

        beta = beta * 2 * eps;

        if (flag_a == 0)
        {
                for (i = 1; i <= 3; i++)
                {
                        u_ab[i] = (a[i] - b[i]) / rab;
                        u_bc[i] = (b[i] - c[i]) / rbc;
                }

                Vector3 *v1 = new Vector3(u_ab[1], u_ab[2], u_ab[3]);
                Vector3 *v2 = new Vector3(u_bc[1], u_bc[2], u_bc[3]);

                Vector3::DotProduct(v1,v2,&c_ab_bc);

                for (i = 1; i <= 3; i++)
                {
                        e_bca[i] = -u_ab[i] + c_ab_bc * u_bc[i];
                }

                Vector3 *vec = new Vector3(e_bca[1], e_bca[2], e_bca[3]);

                vec->Normalize();

                u_bca[1] = vec->X; u_bca[2] = vec->Y; u_bca[3] = vec->Z;

                coef3 = beta * rb / rbc;

                for (i = 1; i <= 3; i++)
                {
                        coef_bc = coef3 * u_bca[i];
                        dsurfb[i][1] = 0;
                        dsurfb[i][2] = coef_bc;
                        dsurfb[i][3] = -coef_bc;
                }

                coef3 = beta * rc / rbc;

                for (i = 1; i <= 3; i++)
                {
                        coef_bc = coef3 * u_bca[i];
                        dsurfc[i][1] = 0;
                        dsurfc[i][2] = -coef_bc;
                        dsurfc[i][3] = coef_bc;
                }
                delete vec;delete v1;delete v2;
                return;
        }

        if (flag_b == 0)
        {
                for (i = 1; i <= 3; i++)
                {
                        u_ac[i] = (a[i] - c[i]) / rac;
                        u_bc[i] = (b[i] - c[i]) / rbc;
                }

                Vector3 *v1 = new Vector3(u_ac[1], u_ac[2], u_ac[3]);
                Vector3 *v2 = new Vector3(u_bc[1], u_bc[2], u_bc[3]);

                Vector3::DotProduct(v1,v2 ,&c_ac_bc);

                for (i = 1; i <= 3; i++)
                {
                        e_cab[i] = -u_bc[i] + c_ac_bc * u_ac[i];
                }

                Vector3 *vec = new Vector3(e_cab[1], e_cab[2], e_cab[3]);

                vec->Normalize();

                u_cab[1] = vec->X; u_cab[2] = vec->Y; u_cab[3] = vec->Z;

                coef3 = beta * ra / rac;

                for (i = 1; i <= 3; i++)
                {
                        coef_ac = coef3 * u_cab[i];
                        dsurfa[i][1] = coef_ac;
                        dsurfa[i][2] = 0;
                        dsurfa[i][3] = -coef_ac;
                }

                coef2 = beta * rc / rac;

                for (i = 1; i <= 3; i++)
                {
                        coef_ac = coef2 * u_cab[i];
                        dsurfc[i][1] = -coef_ac;
                        dsurfc[i][2] = 0;
                        dsurfc[i][3] = coef_ac;
                }
                delete vec;delete v1;delete v2;
                return;
        }

        if (flag_c == 0)
        {
                for (i = 1; i <= 3; i++)
                {
                        u_ab[i] = (a[i] - b[i]) / rab;
                        u_ac[i] = (a[i] - c[i]) / rac;
                }
                Vector3 *v1 = new Vector3(u_ab[1], u_ab[2], u_ab[3]);
                Vector3 *v2 = new Vector3(u_ac[1], u_ac[2], u_ac[3]);

                Vector3::DotProduct(v1,v2 ,&c_ab_ac);

                for (i = 1; i <= 3; i++)
                {
                        e_abc[i] = u_ac[i] - c_ab_ac * u_ab[i];
                }

                Vector3 *vec = new Vector3(e_abc[1], e_abc[2], e_abc[3]);

                vec->Normalize();

                u_abc[1] = vec->X; u_abc[2] = vec->Y; u_abc[3] = vec->Z;

                coef2 = beta * ra / rab;

                for (i = 1; i <= 3; i++)
                {
                        coef_ab = coef2 * u_abc[i];
                        dsurfa[i][1] = coef_ab;
                        dsurfa[i][2] = -coef_ab;
                        dsurfa[i][3] = 0;
                }

                coef2 = beta * rb / rab;

                for (i = 1; i <= 3; i++)
                {
                        coef_ab = coef2 * u_abc[i];
                        dsurfb[i][1] = -coef_ab;
                        dsurfb[i][2] = coef_ab;
                        dsurfb[i][3] = 0;
                }
                delete vec;delete v1;delete v2;
                return;
        }

        for (i = 1; i <= 3; i++)
        {
                u_ab[i] = (a[i] - b[i]) / rab;
                u_ac[i] = (a[i] - c[i]) / rac;
                u_bc[i] = (b[i] - c[i]) / rbc;
        }

        Vector3 *v1 = new Vector3(u_ab[1], u_ab[2], u_ab[3]);
        Vector3 *v2 = new Vector3(u_ac[1], u_ac[2], u_ac[3]);
        Vector3 *v3 = new Vector3(u_bc[1], u_bc[2], u_bc[3]);

        Vector3::DotProduct(v1,v2, &c_ab_ac);
        Vector3::DotProduct(v1,v3, &c_ab_bc);
        Vector3::DotProduct(v2,v3, &c_ac_bc);

        for (i = 1; i <= 3; i++)
        {
                e_abc[i] =  u_ac[i] - c_ab_ac * u_ab[i];
                e_bca[i] = -u_ab[i] + c_ab_bc * u_bc[i];
                e_cab[i] = -u_bc[i] + c_ac_bc * u_ac[i];
        }

        Vector3 *vec1 = new Vector3(e_abc[1], e_abc[2], e_abc[3]);
        vec1->Normalize();
        u_abc[1] = vec1->X; u_abc[2] = vec1->Y; u_abc[3] = vec1->Z;

        delete vec1;
        vec1 = new Vector3(e_bca[1], e_bca[2], e_bca[3]);
        vec1->Normalize();
        u_bca[1] = vec1->X; u_bca[2] = vec1->Y; u_bca[3] = vec1->Z;

        delete vec1;
        vec1 = new Vector3(e_cab[1], e_cab[2], e_cab[3]);
        vec1->Normalize();
        u_cab[1] = vec1->X; u_cab[2] = vec1->Y; u_cab[3] = vec1->Z;

        coef2 = flag_ab * beta * ra / rab;
        coef3 = flag_ac * beta * ra / rac;

        for (i = 1; i <= 3; i++)
        {
                coef_ab = coef2 * u_abc[i];
                coef_ac = coef3 * u_cab[i];
                dsurfa[i][1] = coef_ab + coef_ac;
                dsurfa[i][2] = -coef_ab;
                dsurfa[i][3] = -coef_ac;
        }

        coef2 = flag_ab * beta * rb / rab;
        coef3 = flag_bc * beta * rb / rbc;

        for (i = 1; i <= 3; i++)
        {
                coef_ab = coef2 * u_abc[i];
                coef_bc = coef3 * u_bca[i];
                dsurfb[i][1] = -coef_ab;
                dsurfb[i][2] = coef_ab + coef_bc;
                dsurfb[i][3] = -coef_bc;
        }

        coef2 = flag_ac * beta * rc / rac;
        coef3 = flag_bc * beta * rc / rbc;

        for (i = 1; i <= 3; i++)
        {
                coef_ac = coef2 * u_cab[i];
                coef_bc = coef3 * u_bca[i];
                dsurfc[i][1] = -coef_ac;
                dsurfc[i][2] = -coef_bc;
                dsurfc[i][3] = coef_ac + coef_bc;
        }

        delete v1;delete v2;delete v3;delete vec1;
}

/*!
    \fn SpaceFillMeasure::ThreeSphereVol(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double *surfa, double *surfb, double *surfc, double *vola, double *volb, double *volc)
 */
void SpaceFillMeasure::ThreeSphereVol(double a[], double b[], double c[], double ra, double rb, double rc, double ra2, double rb2, double rc2, double wa, double wb, double wc, double rab, double rac, double rbc, double rab2, double rac2, double rbc2, double *surfa, double *surfb, double *surfc, double *vola, double *volb, double *volc)
{
        double a1, a2, a3, s2, c1, c2, eps = 0.0;
        double seg_ang_abc = 0.0, seg_ang_acb = 0.0, seg_ang_bca = 0.0;
        double ang_abc, ang_acb, ang_bca;
        double cos_abc, cos_acb, cos_bca;
        double sin_abc, sin_acb, sin_bca;
        double ang_dih_abc = 0.0, ang_dih_cab = 0.0, ang_dih_bac = 0.0;
        double s_abc, s_acb, s_bca;
        double val1, val2, val3, l1 = 0.0, l2 = 0.0, l3 = 0.0;
        double val1b, val2b, val3b;
        double rho_ab2, rho_ac2, rho_bc2;

        double center[4];
        double c_ab[4];
        double c_ac[4];
        double c_bc[4];
        double pabc[4];
        double pacb[4];
        double cosine[4];
        double sine[4];
        double angles[4];

        Center2(a, b, ra2, rb2, rab2, c_ab, &l1);
        Center2(a, c, ra2, rc2, rac2, c_ac, &l2);
        Center2(b, c, rb2, rc2, rbc2, c_bc, &l3);

        val1 = l1 * rab;
        val2 = l2 * rac;
        val3 = l3 * rbc;

        val1b = rab - val1;
        val2b = rac - val2;
        val3b = rbc - val3;

        Center3(a, b, c, wa, wb, wc, center);

        TriangleDual(a, b, c, center, &eps, ra2, pabc, pacb);

        Tetra3Noder(a, b, c, pabc, rab, rac, rbc, ra, rb, rc, &seg_ang_abc, &seg_ang_acb, &seg_ang_bca, &ang_dih_abc, &ang_dih_bac, &ang_dih_cab, cosine, sine);

        angles[1] = seg_ang_abc;
        angles[2] = seg_ang_acb;
        angles[3] = seg_ang_bca;

        a1 = ra * (1.0 - 2.0 * ang_dih_abc);
        a2 = 2.0 * seg_ang_abc * val1b;
        a3 = 2.0 * seg_ang_acb * val2b;

        *surfa = 2.0 * PI * ra * (a1 - a2 - a3);

        a1 = rb * (1.0 - 2.0 * ang_dih_bac);
        a2 = 2.0 * seg_ang_abc * val1;
        a3 = 2.0 * seg_ang_bca * val3b;

        *surfb = 2.0 * PI * rb * (a1 - a2 - a3);

        a1 = rc * (1.0 - 2.0 * ang_dih_cab);
        a2 = 2.0 * seg_ang_acb * val2;
        a3 = 2.0 * seg_ang_bca * val3;

        *surfc = 2.0 * PI * rc * (a1 - a2 - a3);

        ang_abc = 2.0 * PI * seg_ang_abc;
        ang_acb = 2.0 * PI * seg_ang_acb;
        ang_bca = 2.0 * PI * seg_ang_bca;

        cos_abc = cosine[1];
        sin_abc = sine[1];
        cos_acb = cosine[2];
        sin_acb = sine[2];
        cos_bca = cosine[3];
        sin_bca = sine[3];

        rho_ab2 = ra2 - val1b * val1b;
        rho_ac2 = ra2 - val2b * val2b;
        rho_bc2 = rb2 - val3b * val3b;

        s_abc = rho_ab2 * (ang_abc - sin_abc * cos_abc);
        s_acb = rho_ac2 * (ang_acb - sin_acb * cos_acb);
        s_bca = rho_bc2 * (ang_bca - sin_bca * cos_bca);

        s2 = ra * (*surfa);
        c1 = val1b * s_abc;
        c2 = val2b * s_acb;

        *vola = (s2 - c1 - c2) / 3.0;

        s2 = rb * (*surfb);
        c1 = val1 * s_abc;
        c2 = val3b * s_bca;

        *volb = (s2 - c1 - c2) / 3.0;

        s2 = rc * (*surfc);
        c1 = val2 * s_acb;
        c2 = val3 * s_bca;

        *volc = (s2 - c1 - c2) / 3.0;
}

/*!
    \fn SpaceFillMeasure::FourSphereVol(double a[], double b[], double c[], double d[], double ra, double rb, double rc, double rd, double ra2, double rb2, double rc2, double rd2, double rab, double rac, double rad, double rbc, double rbd, double rcd, double rab2, double rac2, double rad2, double rbc2, double rbd2, double rcd2, double wa, double wb, double wc, double wd, double eps1, double eps3, double eps5, double eps7, double shabc, double shacb, double shbca, double shabd, double shadb, double shbda, double shacd, double shadc, double shcda, double shbcd, double shbdc, double shcdb, double pacb[], double pabd[], double padc[], double pbcd[], double ang_abc[], double ang_abd[], double ang_acd[], double ang_bcd[], double *surfa, double *surfb, double *surfc, double *surfd, double *vola, double *volb, double *volc, double *vold, double dsurfa[][5], double dsurfb[][5], double dsurfc[][5], double dsurfd[][5], double dvola[][5], double dvolb[][5], double dvolc[][5], double dvold[][5], int option)
 */
void SpaceFillMeasure::FourSphereVol(int a_index,int b_index,int c_index,int d_index,double a[], double b[], double c[], double d[], double ra, double rb, double rc, double rd, double ra2, double rb2, double rc2, double rd2, double rab, double rac, double rad, double rbc, double rbd, double rcd, double rab2, double rac2, double rad2, double rbc2, double rbd2, double rcd2, double wa, double wb, double wc, double wd, double eps1, double eps3, double eps5, double eps7, double shabc, double shacb, double shbca, double shabd, double shadb, double shbda, double shacd, double shadc, double shcda, double shbcd, double shbdc, double shcdb, double pacb[], double pabd[], double padc[], double pbcd[], double ang_abc[], double ang_abd[], double ang_acd[], double ang_bcd[], double *surfa, double *surfb, double *surfc, double *surfd, double *vola, double *volb, double *volc, double *vold, double dsurfa[][5], double dsurfb[][5], double dsurfc[][5], double dsurfd[][5], double dvola[][5], double dvolb[][5], double dvolc[][5], double dvold[][5], int option)
{
        int i;

        double val_ab, val_ac, val_ad, val_bc, val_bd, val_cd;
        double val2_ab, val2_ac, val2_ad, val2_bc, val2_bd, val2_cd;
        double ang1 = 0.0, ang2 = 0.0, ang3 = 0.0, ang4 = 0.0, ang5 = 0.0, ang6 = 0.0;
        double l_ab = 0.0, l_ac = 0.0, l_ad = 0.0, l_bc = 0.0, l_bd = 0.0, l_cd = 0.0;
        double dist1, dist3, dist5, dist7;
        double h1, h3, h5, h7;
        double dab2, dac2, dad2, dbc2, dbd2, dcd2;
        double s1, t1, t2, coef;
        double der_ab, der_ac, der_ad, der_bc, der_bd, der_cd;
        double diff_ab, diff_ac, diff_ad, diff_bc, diff_bd, diff_cd;
        double cap_ab, cap_ac, cap_ad, cap_bc, cap_bd, cap_cd;
        double sin_ab, sin_ac, sin_ad, sin_bc, sin_bd, sin_cd;
        double cos_ab, cos_ac, cos_ad, cos_bc, cos_bd, cos_cd;

        double c_ab[4]; double c_ac[4]; double c_ad[4];
        double c_bc[4]; double c_bd[4]; double c_cd[4];

        double c_abcd[4];

        double u_ab[4]; double u_ac[4]; double u_ad[4];
        double u_bc[4]; double u_bd[4]; double u_cd[4];

        double v_abc[4]; double v_abd[4]; double v_acb[4];
        double v_acd[4]; double v_adb[4]; double v_adc[4];

        double v_bca[4]; double v_bcd[4]; double v_bda[4];
        double v_bdc[4]; double v_cda[4]; double v_cdb[4];

        double vect_ab[4]; double vect_ac[4]; double vect_ad[4];
        double vect_bc[4]; double vect_bd[4]; double vect_cd[4];

        double cof_ab[4]; double cof_ac[4]; double cof_ad[4];
        double cof_bc[4]; double cof_bd[4]; double cof_cd[4];

        Tetra6Dihedral(a, b, c, d, &ang1, &ang2, &ang4, &ang3, &ang5, &ang6);

        Center2(a, b, ra2, rb2, rab2, c_ab, &l_ab);
        Center2(a, c, ra2, rc2, rac2, c_ac, &l_ac);
        Center2(a, d, ra2, rd2, rad2, c_ad, &l_ad);
        Center2(b, c, rb2, rc2, rbc2, c_bc, &l_bc);
        Center2(b, d, rb2, rd2, rbd2, c_bd, &l_bd);
        Center2(c, d, rc2, rd2, rcd2, c_cd, &l_cd);

        val_ab = l_ab * rab;
        val_ac = l_ac * rac;
        val_ad = l_ad * rad;
        val_bc = l_bc * rbc;
        val_bd = l_bd * rbd;
        val_cd = l_cd * rcd;

        val2_ab = rab - val_ab;
        val2_ac = rac - val_ac;
        val2_ad = rad - val_ad;
        val2_bc = rbc - val_bc;
        val2_bd = rbd - val_bd;
        val2_cd = rcd - val_cd;

        *surfa = (-1.0/2.0) * ra + ang1 * val2_ab + ang2 * val2_ac + ang3 * val2_ad;
        *surfa = 2.0* PI * ra * (*surfa);

        *surfb = (-1.0/2.0) * rb + ang1 * val_ab + ang5 * val2_bd + ang4 * val2_bc;
        *surfb = 2.0 * PI * rb * (*surfb);

        *surfc = (-1.0/2.0) * rc + ang2 * val_ac + ang4 * val_bc + ang6 * val2_cd;
        *surfc = 2.0 * PI * rc * (*surfc);

        *surfd = (-1.0/2.0) * rd + ang3 * val_ad + ang6 * val_cd + ang5 * val_bd;
        *surfd = 2.0 * PI * rd * (*surfd);

        //compute volume
        Center4(a, b, c, d, wa, wb, wc, wd, c_abcd);

        dab2 = ra2 - val2_ab * val2_ab;
        dac2 = ra2 - val2_ac * val2_ac;
        dad2 = ra2 - val2_ad * val2_ad;
        dbc2 = rb2 - val2_bc * val2_bc;
        dbd2 = rb2 - val2_bd * val2_bd;
        dcd2 = rc2 - val2_cd * val2_cd;

        dist1 = 0;
        dist3 = 0;
        dist5 = 0;
        dist7 = 0;

        for (i = 1; i <= 3; i++)
        {
                dist1 = dist1 + pow((pacb[i] - c_abcd[i]), 2.0);
                dist3 = dist3 + pow((pabd[i] - c_abcd[i]), 2.0);
                dist5 = dist5 + pow((padc[i] - c_abcd[i]), 2.0);
                dist7 = dist7 + pow((pbcd[i] - c_abcd[i]), 2.0);
        }

        dist1 = sqrt(dist1);
        dist3 = sqrt(dist3);
        dist5 = sqrt(dist5);
        dist7 = sqrt(dist7);

        h1 = dist1 - eps1;
        h3 = dist3 - eps3;
        h5 = dist5 - eps5;
        h7 = dist7 - eps7;

        s1 = -2.0 * PI * dab2 * ang1;
        t1 = shabc * h1;
        t2 = shabd * h3;

        cap_ab = s1 - t1 - t2;

        s1 = -2.0 * PI * dac2 * ang2;
        t1 = shacd * h5;
        t2 = shacb * h1;

        cap_ac = s1 - t1 - t2;

        s1 = -2.0 * PI * dad2 * ang3;
        t1 = shadb * h3;
        t2 = shadc * h5;

        cap_ad = s1 - t1 - t2;

        s1 = -2.0 * PI * dbc2 * ang4;
        t1 = shbca * h1;
        t2 = shbcd * h7;

        cap_bc = s1 - t1 - t2;

        s1 = -2.0 * PI * dbd2 * ang5;
        t1 = shbdc * h7;
        t2 = shbda * h3;

        cap_bd = s1 - t1 - t2;

        s1 = -2.0 * PI * dcd2 * ang6;
        t1 = shcda * h5;
        t2 = shcdb * h7;

        cap_cd = s1 - t1 - t2;

        FILE *fp = fopen("fourvol.txt","w");

        fprintf(fp,"i = %d\n",a_index);
        fprintf(fp,"j = %d\n",b_index);
        fprintf(fp,"k = %d\n",c_index);
        fprintf(fp,"l = %d\n",d_index);

        fprintf(fp,"a(1) = %lf\n",a[1]);
        fprintf(fp,"a(2) = %lf\n",a[2]);
        fprintf(fp,"a(3) = %lf\n",a[3]);

        fprintf(fp,"b(1) = %lf\n",b[1]);
        fprintf(fp,"b(2) = %lf\n",b[2]);
        fprintf(fp,"b(3) = %lf\n",b[3]);

        fprintf(fp,"c(1) = %lf\n",c[1]);
        fprintf(fp,"c(2) = %lf\n",c[2]);
        fprintf(fp,"c(3) = %lf\n",c[3]);

        fprintf(fp,"d(1) = %lf\n",d[1]);
        fprintf(fp,"d(2) = %lf\n",d[2]);
        fprintf(fp,"d(3) = %lf\n",d[3]);

        fprintf(fp,"pacb(1) = %lf\n",pacb[1]);
        fprintf(fp,"pacb(2) = %lf\n",pacb[2]);
        fprintf(fp,"pacb(3) = %lf\n",pacb[3]);

        fprintf(fp,"pabd(1) = %lf\n",pabd[1]);
        fprintf(fp,"pabd(2) = %lf\n",pabd[2]);
        fprintf(fp,"pabd(3) = %lf\n",pabd[3]);

        fprintf(fp,"padc(1) = %lf\n",padc[1]);
        fprintf(fp,"padc(2) = %lf\n",padc[2]);
        fprintf(fp,"padc(3) = %lf\n",padc[3]);

        fprintf(fp,"pbcd(1) = %lf\n",pbcd[1]);
        fprintf(fp,"pbcd(2) = %lf\n",pbcd[2]);
        fprintf(fp,"pbcd(3) = %lf\n",pbcd[3]);

        fprintf(fp,"center(1) = %lf\n",c_abcd[1]);
        fprintf(fp,"center(2) = %lf\n",c_abcd[2]);
        fprintf(fp,"center(3) = %lf\n",c_abcd[3]);

        fprintf(fp,"dist1 = %lf\n",dist1);
        fprintf(fp,"dist3 = %lf\n",dist3);
        fprintf(fp,"dist5 = %lf\n",dist5);
        fprintf(fp,"dist7 = %lf\n",dist7);

        fprintf(fp,"h1 = %lf\n",h1);
        fprintf(fp,"h3 = %lf\n",h3);
        fprintf(fp,"h5 = %lf\n",h5);
        fprintf(fp,"h7 = %lf\n",h7);

        fprintf(fp,"center(3) = %lf\n",c_abcd[3]);
        fprintf(fp,"center(1) = %lf\n",c_abcd[1]);
        fprintf(fp,"center(2) = %lf\n",c_abcd[2]);
        fprintf(fp,"center(3) = %lf\n",c_abcd[3]);
        fprintf(fp,"center(1) = %lf\n",c_abcd[1]);
        fprintf(fp,"center(2) = %lf\n",c_abcd[2]);

        fprintf(fp,"cap_ab = %lf\n",cap_ab);
        fprintf(fp,"cap_ac = %lf\n",cap_ac);
        fprintf(fp,"cap_ad = %lf\n",cap_ad);
        fprintf(fp,"cap_bc = %lf\n",cap_bc);
        fprintf(fp,"cap_bd = %lf\n",cap_bd);
        fprintf(fp,"cap_cd = %lf\n",cap_cd);

        fclose(fp);

        *vola = 2.0 * ra * (*surfa) - val2_ab * cap_ab - val2_ac * cap_ac - val2_ad * cap_ad;
        *vola = *vola / 6.0;

        *volb = 2.0 * rb * (*surfb) - val_ab * cap_ab - val2_bd * cap_bd - val2_bc * cap_bc;
        *volb = *volb / 6.0;

        *volc = 2.0 * rc * (*surfc) - val_ac * cap_ac - val_bc * cap_bc - val2_cd * cap_cd;
        *volc = *volc / 6.0;

        *vold = 2.0 * rd * (*surfd) - val_ad * cap_ad - val_bd * cap_bd - val_cd * cap_cd;
        *vold = *vold / 6.0;

        if (option == 0) return;

        for (i = 1; i <= 3; i++)
        {
                u_ab[i] = (a[i] - b[i]) / rab;
                u_ac[i] = (a[i] - c[i]) / rac;
                u_ad[i] = (a[i] - d[i]) / rad;
                u_bc[i] = (b[i] - c[i]) / rbc;
                u_bd[i] = (b[i] - d[i]) / rbd;
                u_cd[i] = (c[i] - d[i]) / rcd;
        }

        coef = 2.0 * PI * ra;
        der_ab = coef * ang1 * l_ab;
        der_ac = coef * ang2 * l_ac;
        der_ad = coef * ang3 * l_ad;

        for (i = 1; i <= 3; i++)
        {
                diff_ab = der_ab * u_ab[i];
                diff_ac = der_ac * u_ac[i];
                diff_ad = der_ad * u_ad[i];
                dsurfa[i][1] = diff_ab + diff_ac + diff_ad;
                dsurfa[i][2] = -diff_ab;
                dsurfa[i][3] = -diff_ac;
                dsurfa[i][4] = -diff_ad;
        }

        coef = 2.0 * PI * rb;
        der_ab = coef * ang1 * (1.0 - l_ab);
        der_bc = coef * ang4 * l_bc;
        der_bd = coef * ang5 * l_bd;

        for (i = 1; i <= 3; i++)
        {
                diff_ab = der_ab * u_ab[i];
                diff_bc = der_bc * u_bc[i];
                diff_bd = der_bd * u_bd[i];
                dsurfb[i][1] = diff_ab;
                dsurfb[i][2] = -diff_ab + diff_bc + diff_bd;
                dsurfb[i][3] = -diff_bc;
                dsurfb[i][4] = -diff_bd;
        }

        coef = 2.0 * PI * rc;
        der_ac = coef * ang2 * (1.0 - l_ac);
        der_bc = coef * ang4 * (1.0 - l_bc);
        der_cd = coef * ang6 * l_cd;

        for (i = 1; i <= 3; i++)
        {
                diff_ac = der_ac * u_ac[i];
                diff_bc = der_bc * u_bc[i];
                diff_cd = der_cd * u_cd[i];
                dsurfc[i][1] = diff_ac;
                dsurfc[i][2] = diff_bc;
                dsurfc[i][3] = -diff_ac - diff_bc + diff_cd;
                dsurfc[i][4] = -diff_cd;
        }

        coef = 2.0 * PI * rd;
        der_ad = coef * ang3 * (1.0 - l_ad);
        der_bd = coef * ang5 * (1.0 - l_bd);
        der_cd = coef * ang6 * (1.0 - l_cd);

        for (i = 1; i <= 3; i++)
        {
                diff_ad = der_ad * u_ad[i];
                diff_bd = der_bd * u_bd[i];
                diff_cd = der_cd * u_cd[i];
                dsurfd[i][1] = diff_ad;
                dsurfd[i][2] = diff_bd;
                dsurfd[i][3] = diff_cd;
                dsurfd[i][4] = -diff_ad - diff_bd - diff_cd;
        }

            //get volume derivatives

        for (i = 1; i <= 3; i++)
        {
                v_abc[i] = (c_abcd[i] + pacb[i] - 2.0 * c_ab[i]) / 2.0;
                v_abd[i] = (c_abcd[i] + pabd[i] - 2.0 * c_ab[i]) / 2.0;
                v_acb[i] = (c_abcd[i] + pacb[i] - 2.0 * c_ac[i]) / 2.0;
                v_acd[i] = (c_abcd[i] + padc[i] - 2.0 * c_ac[i]) / 2.0;
                v_adb[i] = (c_abcd[i] + pabd[i] - 2.0 * c_ad[i]) / 2.0;
                v_adc[i] = (c_abcd[i] + padc[i] - 2.0 * c_ad[i]) / 2.0;
                v_bca[i] = (c_abcd[i] + pacb[i] - 2.0 * c_bc[i]) / 2.0;
                v_bcd[i] = (c_abcd[i] + pbcd[i] - 2.0 * c_bc[i]) / 2.0;
                v_bda[i] = (c_abcd[i] + pabd[i] - 2.0 * c_bd[i]) / 2.0;
                v_bdc[i] = (c_abcd[i] + pbcd[i] - 2.0 * c_bd[i]) / 2.0;
                v_cda[i] = (c_abcd[i] + padc[i] - 2.0 * c_cd[i]) / 2.0;
                v_cdb[i] = (c_abcd[i] + pbcd[i] - 2.0 * c_cd[i]) / 2.0;
                vect_ab[i] = (pacb[i] + pabd[i] - 2.0 * c_ab[i]) / 2.0;
                vect_ac[i] = (pacb[i] + padc[i] - 2.0 * c_ac[i]) / 2.0;
                vect_ad[i] = (pabd[i] + padc[i] - 2.0 * c_ad[i]) / 2.0;
                vect_bc[i] = (pacb[i] + pbcd[i] - 2.0 * c_bc[i]) / 2.0;
                vect_bd[i] = (pabd[i] + pbcd[i] - 2.0 * c_bd[i]) / 2.0;
                vect_cd[i] = (padc[i] + pbcd[i] - 2.0 * c_cd[i]) / 2.0;
        }

        sin_ab = sin(PI * (ang_abc[1] + ang_abd[1] - ang1));
        sin_ac = sin(PI * (ang_abc[2] + ang_acd[1] - ang2));
        sin_ad = sin(PI * (ang_abd[2] + ang_acd[2] - ang3));
        sin_bc = sin(PI * (ang_abc[3] + ang_bcd[1] - ang4));
        sin_bd = sin(PI * (ang_abd[3] + ang_bcd[2] - ang5));
        sin_cd = sin(PI * (ang_acd[3] + ang_bcd[3] - ang6));

        cos_ab = cos(PI * (ang_abc[1] + ang_abd[1] - ang1));
        cos_ac = cos(PI * (ang_abc[2] + ang_acd[1] - ang2));
        cos_ad = cos(PI * (ang_abd[2] + ang_acd[2] - ang3));
        cos_bc = cos(PI * (ang_abc[3] + ang_bcd[1] - ang4));
        cos_bd = cos(PI * (ang_abd[3] + ang_bcd[2] - ang5));
        cos_cd = cos(PI * (ang_acd[3] + ang_bcd[3] - ang6));

        for (i = 1; i <= 3; i++)
        {
                cof_ab[i] = (shabc * dist1 * v_abc[i] + shabd * dist3 * v_abd[i] - 2.0 * dab2 * sin_ab * vect_ab[i] / cos_ab) / (3.0 * rab);
                cof_ac[i] = (shacb * dist1 * v_acb[i] + shacd * dist5 * v_acd[i] - 2.0 * dac2 * sin_ac * vect_ac[i] / cos_ac) / (3.0 * rac);
                cof_ad[i] = (shadb * dist3 * v_adb[i] + shadc * dist5 * v_adc[i] - 2.0 * dad2 * sin_ad * vect_ad[i] / cos_ad) / (3.0 * rad);
                cof_bc[i] = (shbca * dist1 * v_bca[i] + shbcd * dist7 * v_bcd[i] - 2.0 * dbc2 * sin_bc * vect_bc[i] / cos_bc) / (3.0 * rbc);
                cof_bd[i] = (shbda * dist3 * v_bda[i] + shbdc * dist7 * v_bdc[i] - 2.0 * dbd2 * sin_bd * vect_bd[i] / cos_bd) / (3.0 * rbd);
                cof_cd[i] = (shcda * dist5 * v_cda[i] + shcdb * dist7 * v_cdb[i] - 2.0 * dcd2 * sin_cd * vect_cd[i] / cos_cd) / (3.0 * rcd);
        }

        der_ab = -0.5 * l_ab * cap_ab;
        der_ac = -0.5 * l_ac * cap_ac;
        der_ad = -0.5 * l_ad * cap_ad;

        for (i = 1; i <= 3; i++)
        {
                diff_ab = der_ab * u_ab[i];
                diff_ac = der_ac * u_ac[i];
                diff_ad = der_ad * u_ad[i];
                dvola[i][1] = diff_ab + diff_ac + diff_ad;
                dvola[i][2] = -diff_ab;
                dvola[i][3] = -diff_ac;
                dvola[i][4] = -diff_ad;
        }

        for (i = 1; i <= 3; i++)
        {
                dvola[i][1] = dvola[i][1] + cof_ab[i] + cof_ac[i] + cof_ad[i];
                dvola[i][2] = dvola[i][2] - cof_ab[i];
                dvola[i][3] = dvola[i][3] - cof_ac[i];
                dvola[i][4] = dvola[i][4] - cof_ad[i];
        }

        der_ab = -0.5 * cap_ab * (1.0 - l_ab);
        der_bc = -0.5 * cap_bc * l_bc;
        der_bd = -0.5 * cap_bd * l_bd;

        for (i = 1; i <= 3; i++)
        {
                diff_ab = der_ab * u_ab[i];
                diff_bc = der_bc * u_bc[i];
                diff_bd = der_bd * u_bd[i];
                dvolb[i][1] = diff_ab;
                dvolb[i][2] = -diff_ab + diff_bc + diff_bd;
                dvolb[i][3] = -diff_bc;
                dvolb[i][4] = -diff_bd;
        }

        for (i = 1; i <= 3; i++)
        {
                dvolb[i][1] = dvolb[i][1] - cof_ab[i];
                dvolb[i][2] = dvolb[i][2] + cof_ab[i] + cof_bc[i] + cof_bd[i];
                dvolb[i][3] = dvolb[i][3] - cof_bc[i];
                dvolb[i][4] = dvolb[i][4] - cof_bd[i];
        }

        der_ac = -0.5 * cap_ac * (1.0 - l_ac);
        der_bc = -0.5 * cap_bc * (1.0 - l_bc);
        der_cd = -0.5 * cap_cd * l_cd;

        for (i = 1; i <= 3; i++)
        {
                diff_ac = der_ac * u_ac[i];
                diff_bc = der_bc * u_bc[i];
                diff_cd = der_cd * u_cd[i];
                dvolc[i][1] = diff_ac;
                dvolc[i][2] = diff_bc;
                dvolc[i][3] = -diff_ac - diff_bc + diff_cd;
                dvolc[i][4] = -diff_cd;
        }

        for (i = 1; i <= 3; i++)
        {
                dvolc[i][1] = dvolc[i][1] - cof_ac[i];
                dvolc[i][2] = dvolc[i][2] - cof_bc[i];
                dvolc[i][3] = dvolc[i][3] + cof_ac[i] + cof_bc[i] + cof_cd[i];
                dvolc[i][4] = dvolc[i][4] - cof_cd[i];
        }

        der_ad = -0.5 * cap_ad * (1.0 - l_ad);
        der_bd = -0.5 * cap_bd * (1.0 - l_bd);
        der_cd = -0.5 * cap_cd * (1.0 - l_cd);


        for (i = 1; i <= 3; i++)
        {
                diff_ad = der_ad * u_ad[i];
                diff_bd = der_bd * u_bd[i];
                diff_cd = der_cd * u_cd[i];
                dvold[i][1] = diff_ad;
                dvold[i][2] = diff_bd;
                dvold[i][3] = diff_cd;
                dvold[i][4] = -diff_ad - diff_bd - diff_cd;
        }

        for (i = 1; i <= 3; i++)
        {
                dvold[i][1] = dvold[i][1] - cof_ad[i];
                dvold[i][2] = dvold[i][2] - cof_bd[i];
                dvold[i][3] = dvold[i][3] - cof_cd[i];
                dvold[i][4] = dvold[i][4] + cof_ad[i] + cof_bd[i] + cof_cd[i];
        }
}

