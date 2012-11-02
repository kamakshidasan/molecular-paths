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

#include "size.h"

Size::Size()
{
}

Size::~Size()
{
}

/*!
    \fn Size::int ISort3(int *a, int *b, int *c)
 */
int Size::ISort3(int *a, int *b, int *c)
{
        int swaps = 0, aux;
        if (*a > *b)
        {
                aux = *a;
                *a = *b;
                *b = aux;
                swaps++;
        }
        if (*b > *c)
        {
                aux = *b;
                *b = *c;
                *c = aux;
                swaps++;
                if (*a > *b)
                {
                        aux = *a;
                        *a = *b;
                        *b = aux;
                        swaps++;
                }
        }
        return (swaps);
}

/*!
    \fn Size::ISort4(int *a, int *b, int *c, int *d)
 */
int Size::ISort4(int *a, int *b, int *c, int *d)
{
        int swaps = 0, aux;
        assert ((*a <= *b) && (*b <= *c));
        if (*d < *c)
        {
                if (*d < *a)
                {
                        aux = *c;			//swap (*c, *d);
                        *c = *d;
                        *d = aux;
                        swaps++;
                        aux = *b;			//swap (*b, *c);
                        *b = *c;
                        *c = aux;
                        swaps++;
                        aux = *a;			//swap (*a, *b);
                        *a = *b;
                        *b = aux;
                        swaps++;
                }
                else if (*d < *b)
                {
                        aux = *c;			//swap (*c, *d);
                        *c = *d;
                        *d = aux;
                        swaps++;
                        aux = *b;			//swap (*b, *c);
                        *b = *c;
                        *c = aux;
                        swaps++;
                }
                else
                {
                        aux = *c;			//swap (*c, *d);
                        *c = *d;
                        *d = aux;
                        swaps++;
                }
        }
        return (swaps);
}

/*!
    \fn Size::AlfRatioCompare(mpz_t* a, mpz_t* b,mpz_t* c, mpz_t* d)
 */
int Size::AlfRatioCompare(mpz_t* a, mpz_t* b,mpz_t* c, mpz_t* d)
{
        int ret;
        mpz_t temp;
        mpz_t num;
        mpz_t den;

        if ((mpz_cmp (*a, *c) == 0) && (mpz_cmp (*b, *d) == 0))		//trivially equal
        {
                return (0);
        }
        else if (mpz_sgn (*b) == 0)					//a/b == infinity does not occur, or???
        {
                // Assert_always (FALSE);
                if (mpz_sgn (*d) == 0)					//both are special ratios
                {

                        mpz_init(temp);
                        mpz_sub (temp, *a, *c);
                        ret = mpz_sgn(temp);
                        mpz_clear(temp);
                        return (ret);
                }
                else
                {
                        switch (mpz_sgn (*a))   			//only a/b is special
                        {
                                case -1: return (-1);			//-infinity < anything
                                case  1: return ( 1);			//infinity > anything
                                case  0: return (-mpz_sgn (*c));	//compare 0 to something !=0
                                default: return (-99);			//for lint only
                        }
                }
        }
        else if (mpz_sgn (*d) == 0)  					//c/d == infinity does not occur, or???
        {
                //   Assert_always (FALSE);
                return (-AlfRatioCompare (c, d, a, b));  		//swap ratios!
        }
        else								 //normal case: a/b versus c/d
        {
                mpz_init(num);mpz_init(den);
                mpz_init(temp);
                mpz_mul (num, *a, *d);
                mpz_mul (den, *b, *c);
                mpz_sub (temp,num,den);  				//ad - bc
                ret = mpz_sgn(temp);
                mpz_clear(temp);
                mpz_clear(num);mpz_clear(den);
                return (ret);
        }
}

/*!
    \fn Size::VertSize(std::vector<Vertex> & vertexList, int i, mpz_t *p, mpz_t *q)
 */
void Size::VertSize(std::vector<Vertex> & vertexList, int i, mpz_t *p, mpz_t *q)
{
/*
        Local GMP variables
*/
        mpz_t num, den,temp1;
/*
        Initialise local GMP variables
*/
        mpz_init (num); mpz_init (den);mpz_init (temp1);
/*
        calculate numerator
*/
        mpz_mul(temp1,vertexList[i].V[1],vertexList[i].V[1]);
        mpz_sub(num,vertexList[i].V[4],temp1);
        mpz_mul(temp1,vertexList[i].V[2],vertexList[i].V[2]);
        mpz_sub(num,num,temp1);
        mpz_mul(temp1,vertexList[i].V[3],vertexList[i].V[3]);
        mpz_sub(num,num,temp1);
/*
        set denominator to 1
*/
        mpz_set_ui(den,1);
        //set p and q to num and den respectively
        mpz_set(*p,num); mpz_set(*q,den);
/*
        Clear local GMP variables
*/
        mpz_clear (num); mpz_clear (den);mpz_clear (temp1);
}

/*!
    \fn Size::CheckVertex(std::vector<Vertex> & vertexList, int a, int b)
 */
int Size::CheckVertex(std::vector<Vertex> & vertexList, int a, int b)
{
        int i,testa,coef;

        mpz_t temp[3];
        mpz_t T[3];
        mpz_t d0,d1;
        mpz_t num,den;
/*
        Initialise local GMP variables
*/
        for(i=0;i<3;i++)
        {
                mpz_init(T[i]);
                mpz_init(temp[i]);
        }

        mpz_init(d0);mpz_init(d1);
        mpz_init(num);mpz_init(den);
/*
        perform hidden test
*/
        for(i=0;i<3;i++)
        {
                mpz_sub(temp[i],vertexList[a].V[i+1],vertexList[b].V[i+1]);
                mpz_mul(temp[i],temp[i],vertexList[a].V[i+1]);
                mpz_add(d0,d0,temp[i]);
        }
        coef = 2;
        mpz_mul_si(d1,d0,coef);
        mpz_add(d1,d1,vertexList[b].V[4]);
        mpz_sub(d1,d1,vertexList[a].V[4]);

        switch(mpz_sgn(d1))
        {
                case  1:testa=false;
                break;
                case  0:testa=true+true;
                break;
                case -1:testa=true;
                break;
                default:break;
        }
/*
        clear local GMP variables
*/
        for(i=0;i<3;i++)
        {
                mpz_clear(T[i]);
                mpz_clear(temp[i]);
        }

        mpz_clear(d0);mpz_clear(d1);
        mpz_clear(num);mpz_clear(den);

        return testa;
}

/*!
    \fn Size::EdgeSize(std::vector<Vertex> & vertexList, int a, int b, mpz_t *p, mpz_t *q)
 */
void Size::EdgeSize(std::vector<Vertex> & vertexList, int a, int b, mpz_t *p, mpz_t *q)
{
        int i,j,coef;
        int co_ord[5];
/*
        Local GMP variables
*/
        mpz_t num, den;
        mpz_t temp1,temp2,temp3;
        mpz_t val;
        mpz_t d0,d1,d2,d3,d4,m123;
        mpz_t results[5][5];
/*
        Initialise local GMP variables
*/
        for (i = 0; i < 5; i++)
        {
                co_ord[i] = 0;
        }
        for (i=0; i < 5; i++)
        {
                for(j = 0; j < 5; j++)
                {
                        mpz_init(results[i][j]);
                }
        }

        mpz_init(d0);mpz_init(d1);mpz_init(d2);
        mpz_init(d3);mpz_init(d4);mpz_init(m123);

        mpz_init(temp1); mpz_init(temp2); mpz_init(temp3);
        mpz_init(val);
        mpz_init (num); mpz_init (den);
/*
        First, determine proper ordering for the coordinates.
*/
        co_ord[4] = 4;
        if (mpz_cmp(vertexList[a].V[1],vertexList[b].V[1]))
        {
                co_ord[1] = 1; co_ord[2] = 2; co_ord[3] = 3;
        }
        else if (mpz_cmp(vertexList[a].V[2],vertexList[b].V[2]))
        {
                co_ord[1] = 2; co_ord[2] = 3; co_ord[3] = 1;
        }
        else if (mpz_cmp(vertexList[a].V[3],vertexList[b].V[3]))
        {
                co_ord[1] = 3; co_ord[2] = 1; co_ord[3] = 2;
        }
/*
        Compute and store all the 2 x 2 subdeterminants needed.
*/
        for(i=1;i<5;i++)						//M_{a,b,co_ord[i],0}
        {
                mpz_sub(results[0][i],vertexList[a].V[co_ord[i]], vertexList[b].V[co_ord[i]]);
        }

        for(i=1;i<4;i++)						//M_{a,b,co_ord[i],co_ord[j]}
        {
                for(j = i+1;j<5;j++)
                {
                        mpz_mul(temp1,vertexList[a].V[co_ord[i]],vertexList[b].V[co_ord[j]]);
                        mpz_mul(temp2,vertexList[a].V[co_ord[j]],vertexList[b].V[co_ord[i]]);
                        mpz_sub(results[i][j],temp1,temp2);
                }
        }
/*
        Compute d0 = -2 * M{a,b,co_ord[1],0} *
        [ M{a,b,co_ord[1],0} ^ 2 +
        M{a,b,co_ord[2],0} ^ 2 +
        M{a,b,co_ord[3],0} ^ 2
        ]
*/
        coef = -2.0;
        mpz_mul(temp1,results[0][1],results[0][1]);
        mpz_mul(temp2,results[0][2],results[0][2]);
        mpz_mul(temp3,results[0][3],results[0][3]);
        mpz_add(temp1,temp1,temp2);
        mpz_add(temp1,temp1,temp3);
        mpz_mul(temp1,temp1,results[0][1]);
        mpz_mul_si(d0,temp1,coef);
/*
        Compute d1 = M{a,b,co_ord[1],0} *
        [ 2 * {
        M{a,b,co_ord[3],0} * M{a,b,co_ord[1],co_ord[3]} +
        M{a,b,co_ord[2],0} * M{a,b,co_ord[1],co_ord[2]}
} -
        M{a,b,co_ord[1],0} * M{a,b,co_ord[4],0}
        ]
*/
        coef = 2.0;
        mpz_mul(temp1,results[0][3],results[1][3]);
        mpz_mul(temp2,results[0][2],results[1][2]);
        mpz_add(temp1,temp1,temp2);
        mpz_mul_si(temp1,temp1,coef);
        mpz_mul(temp2,results[0][1],results[0][4]);
        mpz_sub(temp1,temp1,temp2);
        mpz_mul(d1,temp1,results[0][1]);
/*
        Compute d2 = -2 * M{a,b,co_ord[1],co_ord[2]} *
        [ M{a,b,co_ord[3],0} ^2 + M{a,b,co_ord[1],0} ^ 2] -
        M{a,b,co_ord[2],0} *
        [ M{a,b,co_ord[1],0} * M{a,b,co_ord[4],0}
        -2 * M{a,b,co_ord[1],co_ord[3]} * M{a,b,co_ord[3],0}
        ]
*/
        coef = -2.0;
        mpz_mul(temp1,results[0][1],results[0][1]);
        mpz_mul(temp2,results[0][3],results[0][3]);
        mpz_add(temp1,temp1,temp2);
        mpz_mul(temp1,temp1,results[1][2]);
        mpz_mul_si(temp1,temp1,coef);

        mpz_mul(temp2,results[0][1],results[0][4]);
        mpz_mul(temp3,results[0][3],results[1][3]);
        mpz_mul_si(temp3,temp3,coef);
        mpz_add(temp2,temp2,temp3);
        mpz_mul(temp2,temp2,results[0][2]);

        mpz_sub(d2,temp1,temp2);
/*
        Compute d3 = -2 * M{a,b,co_ord[1],co_ord[3]} *
        [ M{a,b,co_ord[2],0} ^2 + M{a,b,co_ord[1],0} ^ 2]
        - M{a,b,co_ord[3],0} *
        [ M{a,b,co_ord[1],0} * M{a,b,co_ord[4],0}
        -2 * M{a,b,co_ord[1],co_ord[2]} * M{a,b,co_ord[2],0}
        ]
*/
        coef = -2.0;
        mpz_mul(temp1,results[0][1],results[0][1]);
        mpz_mul(temp2,results[0][2],results[0][2]);
        mpz_add(temp1,temp1,temp2);
        mpz_mul(temp1,temp1,results[1][3]);
        mpz_mul_si(temp1,temp1,coef);

        mpz_mul(temp2,results[0][1],results[0][4]);
        mpz_mul(temp3,results[0][2],results[1][2]);
        mpz_mul_si(temp3,temp3,coef);
        mpz_add(temp2,temp2,temp3);
        mpz_mul(temp2,temp2,results[0][3]);

        mpz_sub(d3,temp1,temp2);
/*
        Compute d4 =  2 * M{a,b,co_ord[1],0} *
        [ M{a,b,co_ord[1],0} * M{a,b,co_ord[1],co_ord[4]} +
        M{a,b,co_ord[2],0} * M{a,b,co_ord[2],co_ord[4]} +
        M{a,b,co_ord[3],0} * M{a,b,co_ord[3],co_ord[4]}
        ] +
        4 *
        [ M{a,b,co_ord[2],co_ord[3]} *
        { M{a,b,co_ord[3],0} * M{a,b,co_ord[1],co_ord[2]} -
        M{a,b,co_ord[2],0} * M{a,b,co_ord[1],co_ord[3]}
} -
        M{a,b,co_ord[1],0} *
        { M{a,b,co_ord[1],co_ord[3]} ^ 2 + M{a,b,co_ord[1],co_ord[2]} ^ 2 }
        ]

*/
        coef = 2.0;
        mpz_mul(val,results[0][1],results[1][4]);
        mpz_mul(temp1,results[0][2],results[2][4]);
        mpz_mul(temp2,results[0][3],results[3][4]);
        mpz_add(val,val,temp2);
        mpz_add(val,val,temp1);
        mpz_mul(val,val,results[0][1]);
        mpz_mul_si(val,val,coef);

        mpz_mul(temp1,results[0][3],results[1][2]);
        mpz_mul(temp2,results[0][2],results[1][3]);
        mpz_sub(temp1,temp1,temp2);
        mpz_mul(temp1,temp1,results[2][3]);

        mpz_mul(temp2,results[1][2],results[1][2]);
        mpz_mul(temp3,results[1][3],results[1][3]);
        mpz_add(temp2,temp2,temp3);
        mpz_mul(temp2,temp2,results[0][1]);

        coef = 4.0;
        mpz_sub(temp1,temp1,temp2);
        mpz_mul_si(temp1,temp1,coef);

        mpz_add(d4,val,temp1);
/*
        Compute top of fraction = d1 ^ 2 + d2 ^ 2 + d3 ^ 2 - d0 * d4
*/
        mpz_mul(temp1,d1,d1);
        mpz_mul(temp2,d2,d2);
        mpz_mul(temp3,d3,d3);
        mpz_add(temp1,temp1,temp2);
        mpz_add(temp1,temp1,temp3);
        mpz_mul(val,d0,d4);
        mpz_sub(num,temp1,val);
/*
        Compute bottom of fraction = d0 ^ 2
*/
        mpz_mul(den,d0,d0);

        //set p and q to num and den respectively
        mpz_set(*p,num); mpz_set(*q,den);

/*
        Clear local GMP variables
*/
        for (i=0; i < 5; i++)
        {
                for(j = 0; j < 5; j++)
                {
                        mpz_clear(results[i][j]);
                }
        }
        mpz_clear(d0);mpz_clear(d1);mpz_clear(d2);
        mpz_clear(d3);mpz_clear(d4);mpz_clear(m123);

        mpz_clear(val);
        mpz_clear (num); mpz_clear (den);
        mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
}

/*!
    \fn Size::CheckEdge(mpz_t *p_mp,mpz_t *q_mp,mpz_t *r_mp,bool *isAttached)
 */
void Size::CheckEdge(mpz_t *p_mp,mpz_t *q_mp,mpz_t *r_mp,bool *isAttached)
{
        int i,j,coef;
        int co_ord[4];
/*
        Local GMP variables
*/
        mpz_t temp[3];
        mpz_t T[9];
        mpz_t d0,d1;
        mpz_t results[5][5];
        mpz_t gamma,lambda;
/*
        Initialise local GMP variables
*/
        for(i=0;i<3;i++)
        {
                mpz_init(T[3*i]);
                mpz_init(T[3*i+1]);
                mpz_init(T[3*i+2]);
                mpz_init(temp[i]);
        }
        for (i=0; i < 5; i++)
        {
                for(j = 0; j < 5; j++)
                {
                        mpz_init(results[i][j]);
                }
        }
        mpz_init(d0);mpz_init(d1);
        mpz_init(gamma);mpz_init(lambda);
/*
        This is the "hidden1" part
*/
/*
        First, determine proper ordering for the coordinates.
*/
        if (mpz_cmp(p_mp[1],q_mp[1]))
        {
                co_ord[1] = 1; co_ord[2] = 2; co_ord[3] = 3;
        }
        else if (mpz_cmp(p_mp[2],q_mp[2]))
        {
                co_ord[1] = 2; co_ord[2] = 3; co_ord[3] = 1;
        }
        else if (mpz_cmp(p_mp[3],q_mp[3]))
        {
                co_ord[1] = 3; co_ord[2] = 1; co_ord[3] = 2;
        }
/*
        Allocate and pre-compute determinants used more than once
*/
        for(i = 0;i<3;i++)		//M{p,q,co_ord[i+1],0}
        {
                mpz_sub(results[i][0],p_mp[co_ord[i+1]],q_mp[co_ord[i+1]]);
        }

        mpz_mul(temp[0],p_mp[co_ord[1]],q_mp[co_ord[2]]);
        mpz_mul(temp[1],p_mp[co_ord[2]],q_mp[co_ord[1]]);
        mpz_sub(results[3][0],temp[0],temp[1]);			//M{p,q,co_ord[1],co_ord[2]}

        mpz_mul(temp[0],p_mp[co_ord[1]],q_mp[co_ord[3]]);
        mpz_mul(temp[1],p_mp[co_ord[3]],q_mp[co_ord[1]]);
        mpz_sub(results[4][0],temp[0],temp[1]);			//M{p,q,co_ord[1],co_ord[3]}

/*
        Compute Det(Gamma) = M{p,q,co_ord[1],0} *
        [
        M{p,q,co_ord[1],0} ^ 2 +
        M{p,q,co_ord[2],0} ^ 2 +
        M{p,q,co_ord[3],0} ^ 2
        ]
*/
        mpz_mul(T[0],results[0][0],results[0][0]);
        mpz_mul(T[1],results[1][0],results[1][0]);
        mpz_mul(T[2],results[2][0],results[2][0]);
        mpz_add(T[0],T[0],T[1]);
        mpz_add(T[0],T[0],T[2]);
        mpz_mul(gamma,T[0],results[0][0]);
/*
        Compute Det(Lambda) =
        M{p,q,co_ord[1],0} *
        [ M{p,q,co_ord[1],0} * M{p,q,r,co_ord[1],co_ord[4],0} +
        M{p,q,co_ord[2],0} * M{p,q,r,co_ord[2],co_ord[4],0} +
        M{p,q,co_ord[3],0} * M{p,q,r,co_ord[3],co_ord[4],0} -
        2 * {
        M{p,q,co_ord[1],co_ord[3]} * M{p,q,r,co_ord[1],co_ord[3],0} +
        M{p,q,co_ord[1],co_ord[2]} * M{p,q,r,co_ord[1],co_ord[2],0}
}
        ]
        + 2 * M{p,q,r,co_ord[2],co_ord[3],0} *
        [ M{p,q,co_ord[1],co_ord[3]} * M{p,q,co_ord[3],0} -
        M{p,q,co_ord[1],co_ord[3]} * M{p,q,co_ord[2],0}
        ]
*/
        for(i=0;i<3;i++)
        {
                mpz_mul(T[i],p_mp[co_ord[i+1]],q_mp[4]);
                mpz_mul(temp[0],q_mp[co_ord[i+1]],r_mp[4]);
                mpz_mul(temp[1],r_mp[co_ord[i+1]],p_mp[4]);
                mpz_add(T[i],T[i],temp[0]);
                mpz_add(T[i],T[i],temp[1]);

                mpz_mul(temp[0],p_mp[co_ord[i+1]],r_mp[4]);
                mpz_mul(temp[1],q_mp[co_ord[i+1]],p_mp[4]);
                mpz_mul(temp[2],r_mp[co_ord[i+1]],q_mp[4]);
                mpz_add(temp[0],temp[0],temp[1]);
                mpz_add(temp[0],temp[0],temp[2]);

                mpz_sub(T[i],T[i],temp[0]);			//M{p,q,r,co_ord[i+1],co_ord[4],0}
        }

        mpz_mul(T[0],T[0],results[0][0]);
        mpz_mul(T[1],T[1],results[1][0]);
        mpz_mul(T[2],T[2],results[2][0]);
        mpz_add(T[0],T[0],T[1]);
        mpz_add(T[0],T[0],T[2]);

        for(i=1;i<3;i++)
        {
                mpz_mul(T[i],p_mp[co_ord[1]],q_mp[co_ord[i+1]]);
                mpz_mul(temp[0],q_mp[co_ord[1]],r_mp[co_ord[i+1]]);
                mpz_mul(temp[1],r_mp[co_ord[1]],p_mp[co_ord[i+1]]);
                mpz_add(T[i],T[i],temp[0]);
                mpz_add(T[i],T[i],temp[1]);

                mpz_mul(temp[0],p_mp[co_ord[1]],r_mp[co_ord[i+1]]);
                mpz_mul(temp[1],q_mp[co_ord[1]],p_mp[co_ord[i+1]]);
                mpz_mul(temp[2],r_mp[co_ord[1]],q_mp[co_ord[i+1]]);
                mpz_add(temp[0],temp[0],temp[1]);
                mpz_add(temp[0],temp[0],temp[2]);

                mpz_sub(T[i],T[i],temp[0]);			//M{p,q,r,co_ord[1],co_ord[2,3],0}
                mpz_mul(T[i],T[i],results[i+2][0]);

        }

        mpz_add(T[2],T[2],T[1]);
        coef = 2;
        mpz_mul_si(T[2],T[2],coef);

        mpz_sub(T[0],T[0],T[2]);
        mpz_mul(T[0],T[0],results[0][0]);
        //
        mpz_mul(T[2],results[3][0],results[2][0]);
        mpz_mul(T[1],results[4][0],results[1][0]);
        mpz_sub(T[2],T[2],T[1]);
        //
        mpz_mul(T[1],p_mp[co_ord[2]],q_mp[co_ord[3]]);
        mpz_mul(temp[0],q_mp[co_ord[2]],r_mp[co_ord[3]]);
        mpz_mul(temp[1],r_mp[co_ord[2]],p_mp[co_ord[3]]);
        mpz_add(T[1],T[1],temp[0]);
        mpz_add(T[1],T[1],temp[1]);

        mpz_mul(temp[0],p_mp[co_ord[2]],r_mp[co_ord[3]]);
        mpz_mul(temp[1],q_mp[co_ord[2]],p_mp[co_ord[3]]);
        mpz_mul(temp[2],r_mp[co_ord[2]],q_mp[co_ord[3]]);
        mpz_add(temp[0],temp[0],temp[1]);
        mpz_add(temp[0],temp[0],temp[2]);

        mpz_sub(T[1],T[1],temp[0]);				//M{p,q,r,co_ord[2],co_ord[3],0}
        //
        mpz_mul_si(T[1],T[1],coef);

        mpz_mul(T[2],T[2],T[1]);
        mpz_add(lambda,T[0],T[2]);

        mpz_mul(gamma,gamma,lambda);

        switch(mpz_sgn(gamma))
        {
                case 1:*isAttached = true;
                break;
                case 0:*isAttached = true;
                break;
                case -1:*isAttached = false;
                break;
                default:break;
        }
/*
        Clear local GMP variables
*/
        for(i=0;i<3;i++)
        {
                mpz_clear(T[3*i]);
                mpz_clear(T[3*i+1]);
                mpz_clear(T[3*i+2]);
                mpz_clear(temp[i]);
        }
        for (i=0; i < 5; i++)
        {
                for(j = 0; j < 5; j++)
                {
                        mpz_init(results[i][j]);
                }
        }
        mpz_clear(d0);mpz_clear(d1);
        mpz_clear(gamma);mpz_clear(lambda);
}

/*!
    \fn Size::TrigSize(std::vector<Vertex> & vertexList, int a, int b, int c, mpz_t *p, mpz_t *q)
 */
void Size::TrigSize(std::vector<Vertex> & vertexList, int a, int b, int c, mpz_t *p, mpz_t *q)
{
        int i,j,coef;
/*
        Local GMP variables
*/
        mpz_t num, den;
        mpz_t temp1,temp2,temp3;
        mpz_t val;
        mpz_t d0,d1,d2,d3,d4,m123,m124,m134,m234;
        mpz_t results[5][5];
/*
        Initialise local GMP variables
*/
        for (i=0; i < 5; i++)
        {
                for(j = 0; j < 5; j++)
                {
                        mpz_init(results[i][j]);
                }
        }
        mpz_init(d0);mpz_init(d1);mpz_init(d2);mpz_init(d3);mpz_init(d4);
        mpz_init(m123);mpz_init(m124);mpz_init(m134);mpz_init(m234);

        mpz_init (num); mpz_init (den);
        mpz_init(temp1); mpz_init(temp2); mpz_init(temp3);
        mpz_init(val);
/*
        Determinant usage
        Used once:
        124,134,234

        Used more than once:
        120,130,140,230,240,340,123
*/
/*
        Allocate and pre-compute determinants used more than once
*/
        for(i = 0; i < 3; i++)				//M_{a,b,c,i+1,j+2,0} checked
        {
                for(j = i; j < 3; j++)
                {
                        mpz_mul(temp1,vertexList[a].V[i+1],vertexList[b].V[j+2]);
                        mpz_mul(temp2,vertexList[b].V[i+1],vertexList[c].V[j+2]);
                        mpz_add(temp1,temp1,temp2);
                        mpz_mul(temp2,vertexList[c].V[i+1],vertexList[a].V[j+2]);
                        mpz_add(temp1,temp1,temp2);

                        mpz_mul(temp2,vertexList[a].V[i+1],vertexList[c].V[j+2]);
                        mpz_mul(temp3,vertexList[b].V[i+1],vertexList[a].V[j+2]);
                        mpz_add(temp2,temp2,temp3);
                        mpz_mul(temp3,vertexList[c].V[i+1],vertexList[b].V[j+2]);
                        mpz_add(temp2,temp2,temp3);

                        mpz_sub(results[i][j],temp1,temp2);
                }
        }
        //M_{a,b,c,1,2,3}
        mpz_mul(val,vertexList[b].V[2],vertexList[c].V[3]);
        mpz_mul(temp1,vertexList[c].V[2],vertexList[b].V[3]);
        mpz_sub(val,val,temp1);
        mpz_mul(val,val,vertexList[a].V[1]);

        mpz_mul(temp1,vertexList[b].V[1],vertexList[c].V[3]);
        mpz_mul(temp2,vertexList[c].V[1],vertexList[b].V[3]);
        mpz_sub(temp1,temp1,temp2);
        mpz_mul(temp1,temp1,vertexList[a].V[2]);

        mpz_mul(temp2,vertexList[b].V[1],vertexList[c].V[2]);
        mpz_mul(temp3,vertexList[c].V[1],vertexList[b].V[2]);
        mpz_sub(temp2,temp2,temp3);
        mpz_mul(temp2,temp2,vertexList[a].V[3]);

        mpz_sub(val,val,temp1);
        mpz_add(m123,val,temp2);

        for(i = 1;i<3;i++)		//checked
        {
                for(j= i+1;j<4;j++)
                {
                        mpz_mul(val,vertexList[b].V[j],vertexList[c].V[4]);
                        mpz_mul(temp1,vertexList[c].V[j],vertexList[b].V[4]);
                        mpz_sub(val,val,temp1);
                        mpz_mul(val,val,vertexList[a].V[i]);

                        mpz_mul(temp1,vertexList[b].V[i],vertexList[c].V[4]);
                        mpz_mul(temp2,vertexList[c].V[i],vertexList[b].V[4]);
                        mpz_sub(temp1,temp1,temp2);
                        mpz_mul(temp1,temp1,vertexList[a].V[j]);

                        mpz_mul(temp2,vertexList[b].V[i],vertexList[c].V[j]);
                        mpz_mul(temp3,vertexList[c].V[i],vertexList[b].V[j]);
                        mpz_sub(temp2,temp2,temp3);
                        mpz_mul(temp2,temp2,vertexList[a].V[4]);

                        mpz_sub(val,val,temp1);
                        mpz_add(results[4][i+j-2],val,temp2);
                }
        }
        mpz_set(m124,results[4][1]);			//M_{a,b,c,1,2,4}
        mpz_set(m134,results[4][2]);			//M_{a,b,c,1,3,4}
        mpz_set(m234,results[4][3]);			//M_{a,b,c,2,3,4}
/*
        Compute d0 = 4 * [ M_{1,2,0} ^ 2 + M_{1,3,0} ^ 2 +M_{2,3,0} ^ 2 ]
*/
        coef = 4.0;
        mpz_mul(temp1,results[0][0],results[0][0]);	//M_{1,2,0}
        mpz_mul(temp2,results[0][1],results[0][1]);	//M_{1,3,0}
        mpz_mul(temp3,results[1][1],results[1][1]);	//M_{2,3,0}
        mpz_add(temp1,temp1,temp2);
        mpz_add(temp1,temp1,temp3);
        mpz_mul_si(d0,temp1,coef);
/*
        Compute d1 = -2 * [
        M_{1,3,0} * M_{3,4,0} +
        M_{1,2,0} * M_{2,4,0} -
        2 * M_{1,2,3} * M_{2,3,0}
        ]
*/
        coef = -2.0;
        mpz_mul(temp1,results[0][1],results[2][2]);	//M_{1,3,0} * M_{3,4,0}
        mpz_mul(temp2,results[0][0],results[1][2]);	//M_{1,2,0} * M_{2,4,0}
        mpz_mul(temp3,m123,results[1][1]);		//M_{1,2,3} * M_{2,3,0}
        mpz_mul_si(temp3,temp3,coef);
        mpz_add(temp1,temp1,temp2);
        mpz_add(temp1,temp1,temp3);
        mpz_mul_si(d1,temp1,coef);
/*
        Compute d2 =  2 * [
        M_{1,2,0} * M_{1,4,0} -
        M_{2,3,0} * M_{3,4,0} -
        2 * M_{1,2,3} * M_{1,3,0}
        ]
*/
        coef = -2.0;
        mpz_mul(temp1,results[0][0],results[0][2]);	//M_{1,2,0} * M_{1,4,0}
        mpz_mul(temp2,results[1][1],results[2][2]);	//M_{2,3,0} * M_{3,4,0}
        mpz_sub(temp1,temp1,temp2);
        mpz_mul(temp3,m123,results[0][1]);		//M_{1,2,3} * M_{1,3,0}
        mpz_mul_si(temp3,temp3,coef);
        mpz_add(temp1,temp1,temp3);
        coef = 2.0;
        mpz_mul_si(d2,temp1,coef);
/*
        Compute d3 =  2 * [
        M_{2,3,0} * M_{2,4,0} +
        M_{1,3,0} * M_{1,4,0} +
        2 * M_{1,2,3} * M_{1,2,0}
        ]
*/
        coef = 2.0;
        mpz_mul(temp1,results[1][1],results[1][2]);     //M_{2,3,0} * M_{2,4,0}
        mpz_mul(temp2,results[0][1],results[0][2]);     //M_{1,3,0} * M_{1,4,0}
        mpz_mul(temp3,m123,results[0][0]);     		//M_{1,2,3} * M_{1,2,0}
        mpz_mul_si(temp3,temp3,coef);
        mpz_add(temp1,temp1,temp2);
        mpz_add(temp1,temp1,temp3);
        mpz_mul_si(d3,temp1,coef);
/*
        Compute d4 = -4 * [
        M_{1,2,0} * M_{1,2,4} +
        M_{1,3,0} * M_{1,3,4} +
        M_{2,3,0} * M_{2,3,4} -
        2 * M_{1,2,3} ^ 2
        ]
*/
        coef = -2.0;
        mpz_mul(temp1,results[0][0],m124);     //M_{1,2,0} * M_{1,2,4}
        mpz_mul(temp2,results[0][1],m134);     //M_{1,3,0} * M_{1,3,4}
        mpz_mul(temp3,results[1][1],m234);     //M_{2,3,0} * M_{2,3,4}
        mpz_mul(val,m123,m123);
        mpz_mul_si(val,val,coef);
        mpz_add(temp1,temp1,temp2);
        mpz_add(temp1,temp1,temp3);
        mpz_add(temp1,temp1,val);
        coef = -4.0;
        mpz_mul_si(d4,temp1,coef);
/*
        Compute top of fraction = d1 ^ 2 + d2 ^ 2 + d3 ^ 2 - d0 * d4
*/
        mpz_mul(temp1,d1,d1);
        mpz_mul(temp2,d2,d2);
        mpz_mul(temp3,d3,d3);
        mpz_mul(val,d0,d4);
        mpz_add(temp1,temp1,temp2);
        mpz_add(temp1,temp1,temp3);
        mpz_sub(num,temp1,val);
/*
        Compute bottom of fraction = d0 ^ 2
*/
        mpz_mul(den,d0,d0);

        mpz_set(*p,num); mpz_set(*q,den);
/*
        Clear local GMP variables
*/
        for (i=0; i < 5; i++)
        {
                for(j = 0; j < 5; j++)
                {
                        mpz_clear(results[i][j]);
                }
        }
        mpz_clear(d0);mpz_clear(d1);mpz_clear(d2);mpz_clear(d3);mpz_clear(d4);
        mpz_clear(m123);mpz_clear(m124);mpz_clear(m134);mpz_clear(m234);

        mpz_clear(val);
        mpz_clear (num); mpz_clear (den);
        mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
}

/*!
    \fn Size::CheckTriangle(mpz_t *p_mp,mpz_t *q_mp,mpz_t *r_mp,mpz_t *s_mp, bool *isAttached)
 */
void Size::CheckTriangle(mpz_t *p_mp,mpz_t *q_mp,mpz_t *r_mp,mpz_t *s_mp, bool *isAttached)
{
        int i,j,coef;
/*
        Local GMP variables
*/
        mpz_t temp[3];
        mpz_t T[9];
        mpz_t d0,d1;
        mpz_t results[5][5];
        mpz_t gamma,lambda,m123;
        mpz_t rs4,rsi,rsj,rsij,rsi4,rsj4;
/*
        Initialise local GMP variables
*/
        for(i=0;i<3;i++)
        {
                mpz_init(T[3*i]);
                mpz_init(T[3*i+1]);
                mpz_init(T[3*i+2]);
                mpz_init(temp[i]);
        }
        for (i=0; i < 5; i++)
        {
                for(j = 0; j < 5; j++)
                {
                        mpz_init(results[i][j]);
                }
        }
        mpz_init(d0);mpz_init(d1);
        mpz_init(gamma);mpz_init(lambda);mpz_init(m123);
        mpz_init(rs4);mpz_init(rsi);mpz_init(rsj);mpz_init(rsi4);mpz_init(rsj4);mpz_init(rsij);
/*
        Check if triangle is attached
*/
/*
        Compute Det(Gamma) = M{2,3,0} ^ 2 + M{1,3,0} ^ 2 + M{1,2,0} ^ 2
*/
        for(i=1;i<3;i++)
        {
                for(j=i+1;j<4;j++)
                {
                        mpz_mul(d0,p_mp[i],q_mp[j]);
                        mpz_mul(temp[0],q_mp[i],r_mp[j]);
                        mpz_mul(temp[1],r_mp[i],p_mp[j]);
                        mpz_add(d0,d0,temp[0]);
                        mpz_add(d0,d0,temp[1]);

                        mpz_mul(temp[0],p_mp[i],r_mp[j]);
                        mpz_mul(temp[1],q_mp[i],p_mp[j]);
                        mpz_mul(temp[2],r_mp[i],q_mp[j]);
                        mpz_add(temp[0],temp[0],temp[1]);
                        mpz_add(temp[0],temp[0],temp[2]);

                        mpz_sub(T[i+j],d0,temp[0]);
                }
        }

        mpz_set(results[0][0],T[5]);		//M{p,q,r,2,3,0}
        mpz_set(results[1][0],T[4]);		//M{p,q,r,1,3,0}
        mpz_set(results[2][0],T[3]);		//M{p,q,r,1,2,0}

        for(i=0;i<3;i++)
        {
                mpz_mul(T[i],results[i][0],results[i][0]);
        }
        mpz_add(T[0],T[0],T[1]);
        mpz_add(gamma,T[0],T[2]);
/*
        Compute Det(Lambda)
        M{2,3,0} * M{2,3,4,0} +
        M{1,3,0} * M{1,3,4,0} +
        M{1,2,0} * M{1,2,4,0} -
        2 * M{1,2,3} * M{1,2,3,0}
*/
        mpz_sub(rs4,r_mp[4],s_mp[4]);				//rs4 = r[4] - s[4]

        for(i=1;i<3;i++)
        {
                for(j=i+1;j<4;j++)
                {

                        mpz_sub(rsi,r_mp[i],s_mp[i]);		//rsi = r[i] - s[i]
                        mpz_sub(rsj,r_mp[j],s_mp[j]);		//rsj = r[j] - s[j]

                        mpz_mul(rsi4,r_mp[i],s_mp[4]);
                        mpz_mul(T[0],s_mp[i],r_mp[4]);
                        mpz_sub(rsi4,rsi4,T[0]);		//rsi4 = r[i]s[4] - s[i]r[4]

                        mpz_mul(rsj4,r_mp[j],s_mp[4]);
                        mpz_mul(T[0],s_mp[j],r_mp[4]);
                        mpz_sub(rsj4,rsj4,T[0]);		//rsj4 = r[j]s[4] - s[j]r[4]

                        mpz_mul(rsij,r_mp[i],s_mp[j]);
                        mpz_mul(T[0],s_mp[i],r_mp[j]);
                        mpz_sub(rsij,rsij,T[0]);		//rsij = r[i]s[j] - s[i]r[j]

                        //1st term
                        mpz_mul(T[0],q_mp[j],rs4);
                        mpz_mul(T[1],q_mp[4],rsj);
                        mpz_add(T[0],T[0],rsj4);
                        mpz_sub(T[0],T[0],T[1]);
                        mpz_mul(T[0],T[0],p_mp[i]);

                        //2nd term
                        mpz_mul(T[1],q_mp[i],rs4);
                        mpz_mul(T[2],q_mp[4],rsi);
                        mpz_add(T[1],T[1],rsi4);
                        mpz_sub(T[1],T[1],T[2]);
                        mpz_mul(T[1],T[1],p_mp[j]);

                        //3rd term
                        mpz_mul(T[2],q_mp[i],rsj);
                        mpz_mul(T[3],q_mp[j],rsi);
                        mpz_add(T[2],T[2],rsij);
                        mpz_sub(T[2],T[2],T[3]);
                        mpz_mul(T[2],T[2],p_mp[4]);

                        //4th term
                        mpz_mul(T[3],q_mp[i],rsj4);
                        mpz_mul(T[4],q_mp[j],rsi4);
                        mpz_sub(T[3],T[3],T[4]);
                        mpz_mul(T[4],q_mp[4],rsij);
                        mpz_add(T[3],T[3],T[4]);

                        //result = 1st - 2nd + 3rd - 4th
                        mpz_add(T[0],T[0],T[2]);
                        mpz_add(T[1],T[1],T[3]);
                        mpz_sub(T[i+j+3],T[0],T[1]);
                }
        }

        mpz_mul(T[0],T[8],results[0][0]);		//M{p,q,r,2,3,0} * M{p,q,r,s,2,3,4,0}
        mpz_mul(T[1],T[7],results[1][0]);		//M{p,q,r,1,3,0} * M{p,q,r,s,1,3,4,0}
        mpz_mul(T[2],T[6],results[2][0]);		//M{p,q,r,1,2,0} * M{p,q,r,s,1,2,4,0}

        mpz_add(T[0],T[0],T[1]);
        mpz_add(T[0],T[0],T[2]);

        mpz_set(T[8],T[0]);				//store T[0]

        coef = 2;

        //M_{a,b,c,d,1,2,3,0}
        //use variables as follows rs4 = rs3 i = 1 j = 2 i.e rsi = rs1 etc

        mpz_sub(rs4,r_mp[3],s_mp[3]);

        mpz_sub(rsi,r_mp[1],s_mp[1]);		//rsi = r[1] - s[1]
        mpz_sub(rsj,r_mp[2],s_mp[2]);		//rsj = r[2] - s[2]

        mpz_mul(rsi4,r_mp[1],s_mp[3]);
        mpz_mul(T[0],s_mp[1],r_mp[3]);
        mpz_sub(rsi4,rsi4,T[0]);		//rsi4 = r[1]s[3] - s[1]r[3]

        mpz_mul(rsj4,r_mp[2],s_mp[3]);
        mpz_mul(T[0],s_mp[2],r_mp[3]);
        mpz_sub(rsj4,rsj4,T[0]);		//rsj4 = r[2]s[3] - s[2]r[3]

        mpz_mul(rsij,r_mp[1],s_mp[2]);
        mpz_mul(T[0],s_mp[1],r_mp[2]);
        mpz_sub(rsij,rsij,T[0]);		//rsij = r[1]s[2] - s[1]r[2]

        //1st term
        mpz_mul(T[0],q_mp[2],rs4);
        mpz_mul(T[1],q_mp[3],rsj);
        mpz_add(T[0],T[0],rsj4);
        mpz_sub(T[0],T[0],T[1]);
        mpz_mul(T[0],T[0],p_mp[1]);

        //2nd term
        mpz_mul(T[1],q_mp[1],rs4);
        mpz_mul(T[2],q_mp[3],rsi);
        mpz_add(T[1],T[1],rsi4);
        mpz_sub(T[1],T[1],T[2]);
        mpz_mul(T[1],T[1],p_mp[2]);

        //3rd term
        mpz_mul(T[2],q_mp[1],rsj);
        mpz_mul(T[3],q_mp[2],rsi);
        mpz_add(T[2],T[2],rsij);
        mpz_sub(T[2],T[2],T[3]);
        mpz_mul(T[2],T[2],p_mp[3]);

        //4th term
        mpz_mul(T[3],q_mp[1],rsj4);
        mpz_mul(T[4],q_mp[2],rsi4);
        mpz_sub(T[3],T[3],T[4]);
        mpz_mul(T[4],q_mp[3],rsij);
        mpz_add(T[3],T[3],T[4]);

        //result = 1st - 2nd + 3rd - 4th
        mpz_add(T[0],T[0],T[2]);
        mpz_add(T[1],T[1],T[3]);
        mpz_sub(T[0],T[0],T[1]);

        //M_{a,b,c,1,2,3}
        mpz_mul(d0,q_mp[2],r_mp[3]);
        mpz_mul(temp[0],r_mp[2],q_mp[3]);
        mpz_sub(d0,d0,temp[0]);
        mpz_mul(d0,d0,p_mp[1]);

        mpz_mul(temp[0],q_mp[1],r_mp[3]);
        mpz_mul(temp[1],r_mp[1],q_mp[3]);
        mpz_sub(temp[0],temp[0],temp[1]);
        mpz_mul(temp[0],temp[0],p_mp[2]);

        mpz_mul(temp[1],q_mp[1],r_mp[2]);
        mpz_mul(temp[2],r_mp[1],q_mp[2]);
        mpz_sub(temp[1],temp[1],temp[2]);
        mpz_mul(temp[1],temp[1],p_mp[3]);

        mpz_sub(d0,d0,temp[0]);
        mpz_add(m123,d0,temp[1]);			//M{p,q,r,1,2,3}

        mpz_mul(T[0],T[0],m123);			//M{p,q,r,1,2,3} * M{p,q,r,s,1,2,3,0}
        mpz_mul_si(T[0],T[0],coef);

        mpz_sub(lambda,T[8],T[0]);

        mpz_mul(gamma,gamma,lambda);

        switch(mpz_sgn(gamma))
        {
                case 1:*isAttached = true;
                break;
                case 0:*isAttached = true;
                break;
                case -1:*isAttached = false;
                break;
                default:break;
        }

/*
        Clear local GMP variables
*/
        for(i=0;i<3;i++)
        {
                mpz_clear(T[3*i]);
                mpz_clear(T[3*i+1]);
                mpz_clear(T[3*i+2]);
                mpz_clear(temp[i]);
        }
        for (i=0; i < 5; i++)
        {
                for(j = 0; j < 5; j++)
                {
                        mpz_init(results[i][j]);
                }
        }
        mpz_clear(d0);mpz_clear(d1);
        mpz_clear(gamma);mpz_clear(lambda);mpz_clear(m123);
        mpz_clear(rs4);mpz_clear(rsi);mpz_clear(rsj);mpz_clear(rsi4);mpz_clear(rsj4);mpz_clear(rsij);
}

/*!
    \fn Size::TetraSize(DeluanayComplex *delcx, std::vector<Vertex> & vertexList, int a, int b, int c, int d, mpz_t *p, mpz_t *q, int index)
 */
void Size::TetraSize(DeluanayComplex *delcx, std::vector<Vertex> & vertexList, int a, int b, int c, int d, mpz_t *p, mpz_t *q, int index)
{	int i,j,k,coef;
        int tri[4],pair[6];
/*
        Local GMP variables
*/
        mpz_t num, den;
        mpz_t temp1,temp2,temp3;
        mpz_t val1,val2,val3,val4;
        mpz_t Dabc,Dabd,Dacd,Dbcd;
        mpz_t Det1,Det2,Det3,Det4,Dabcd;

        mpz_t Sab[4],Sac[4],Sad[4],Sbc[4],Sbd[4],Scd[4];
        mpz_t Dab[4],Dac[4],Dad[4],Dbc[4],Dbd[4],Dcd[4];
        mpz_t Sa[4],Sb[4],Sc[4],Sd[4];
        mpz_t Sam1[4],Sbm1[4],Scm1[4],Sdm1[4];
        mpz_t Ta[4],Tb[4],Tc[4],Td[4];
        mpz_t Tam1[4],Tbm1[4],Tcm1[4],Tdm1[4];
        mpz_t Ua[4],Ub[4],Uc[4],Ud[4];
        mpz_t Deter[4];
/*
        Initialise local GMP variables
*/
        for (i = 0; i < 4; i++)
        {
                tri[i] = 0;
                if(2*i<6) {pair[2*i] = pair[2*i+1] = 0;}
        }

        mpz_init (num); mpz_init (den);

        mpz_init(temp1); mpz_init(temp2); mpz_init(temp3);
        mpz_init(val1); mpz_init(val2); mpz_init(val3); mpz_init(val4);

        for (i=0; i < 4; i++)
        {
                mpz_init(Sab[i]); mpz_init(Sac[i]); mpz_init(Sad[i]);
                mpz_init(Sbc[i]); mpz_init(Sbd[i]); mpz_init(Scd[i]);
                mpz_init(Dab[i]); mpz_init(Dac[i]); mpz_init(Dad[i]);
                mpz_init(Dbc[i]); mpz_init(Dbd[i]); mpz_init(Dcd[i]);
                mpz_init(Sa[i]); mpz_init(Sb[i]); mpz_init(Sc[i]);
                mpz_init(Sd[i]);
                mpz_init(Sam1[i]); mpz_init(Sbm1[i]); mpz_init(Scm1[i]);
                mpz_init(Sdm1[i]);
                mpz_init(Ta[i]); mpz_init(Tb[i]); mpz_init(Tc[i]);
                mpz_init(Td[i]);
                mpz_init(Tam1[i]); mpz_init(Tbm1[i]); mpz_init(Tcm1[i]);
                mpz_init(Tdm1[i]);
                mpz_init(Ua[i]); mpz_init(Ub[i]); mpz_init(Uc[i]);
                mpz_init(Ud[i]);
                mpz_init(Deter[i]);
        }

        mpz_init(Dabc); mpz_init(Dabd); mpz_init(Dacd); mpz_init(Dbcd);
        mpz_init(Det1); mpz_init(Det2); mpz_init(Det3); mpz_init(Dabcd);
        mpz_init(Det4);
/*
        1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
                                                            | n(i)  n(j) |
        for all i in [1,2] and all j in [i+1,3]
*/
        for (i=1;  i<3; i++)
        {
                for (j=i+1; j<4 ; j++)
                {
                        k=i+j-2;
                        mpz_mul(temp1,vertexList[a].V[j],vertexList[b].V[i]);
                        mpz_mul(temp2,vertexList[a].V[i],vertexList[b].V[j]);
                        mpz_sub(Sab[k],temp2,temp1);
                        mpz_mul(temp1,vertexList[a].V[j],vertexList[c].V[i]);
                        mpz_mul(temp2,vertexList[a].V[i],vertexList[c].V[j]);
                        mpz_sub(Sac[k],temp2,temp1);
                        mpz_mul(temp1,vertexList[a].V[j],vertexList[d].V[i]);
                        mpz_mul(temp2,vertexList[a].V[i],vertexList[d].V[j]);
                        mpz_sub(Sad[k],temp2,temp1);
                        mpz_mul(temp1,vertexList[b].V[j],vertexList[c].V[i]);
                        mpz_mul(temp2,vertexList[b].V[i],vertexList[c].V[j]);
                        mpz_sub(Sbc[k],temp2,temp1);
                        mpz_mul(temp1,vertexList[b].V[j],vertexList[d].V[i]);
                        mpz_mul(temp2,vertexList[b].V[i],vertexList[d].V[j]);
                        mpz_sub(Sbd[k],temp2,temp1);
                        mpz_mul(temp1,vertexList[c].V[j],vertexList[d].V[i]);
                        mpz_mul(temp2,vertexList[c].V[i],vertexList[d].V[j]);
                        mpz_sub(Scd[k],temp2,temp1);
                }
        }
/*
        Now compute all Minors
        Sq(i+j-2) = M(m,n,p,i,j,0) = Det| m(i) m(j) 1 |
                                        | n(i) n(j) 1 |
                                        | p(i) p(j) 1 |

        and all Minors
        Det(i+j-2) = M(m,n,p,q,i,j,4,0) =    Det| m(i) m(j) m(4) 1 |
                                                | n(i) n(j) n(4) 1 |
                                                | p(i) p(j) p(4) 1 |
                                                | q(i) q(j) q(4) 1 |

        m,n,p,q are the four vertices of the tetrahedron, i and j correspond
        to two of the coordinates of the vertices, and m(4) refers to the
        "weight" of vertices m
*/
        for (i=1; i<4; i++)
        {
                mpz_sub(temp1,Scd[i],Sbd[i]);
                mpz_add(Sa[i],temp1,Sbc[i]);
                mpz_mul(temp2,Sa[i],vertexList[a].V[4]);
                mpz_sub(temp1,Scd[i],Sad[i]);
                mpz_add(Sb[i],temp1,Sac[i]);
                mpz_mul(temp3,Sb[i],vertexList[b].V[4]);
                mpz_sub(temp2,temp2,temp3);
                mpz_sub(temp1,Sbd[i],Sad[i]);
                mpz_add(Sc[i],temp1,Sab[i]);
                mpz_mul(temp3,Sc[i],vertexList[c].V[4]);
                mpz_add(temp2,temp2,temp3);
                mpz_sub(temp1,Sbc[i],Sac[i]);
                mpz_add(Sd[i],temp1,Sab[i]);
                mpz_mul(temp3,Sd[i],vertexList[d].V[4]);
                mpz_sub(Deter[i],temp2,temp3);
                mpz_neg(Sam1[i],Sa[i]);
                mpz_neg(Sbm1[i],Sb[i]);
                mpz_neg(Scm1[i],Sc[i]);
                mpz_neg(Sdm1[i],Sd[i]);
        }
/*
        Now compute the determinant needed to compute the radius of the
        circumsphere of the tetrahedron :

        Det1 = Minor(a,b,c,d,4,2,3,0)
        Det2 = Minor(a,b,c,d,1,3,4,0)
        Det3 = Minor(a,b,c,d,1,2,4,0)
        Det4 = Minor(a,b,c,d,1,2,3,0)
*/
        mpz_set(Det1,Deter[3]);
        mpz_set(Det2,Deter[2]);
        mpz_set(Det3,Deter[1]);

        mpz_mul(temp1,vertexList[a].V[1],Sa[3]);
        mpz_mul(temp2,vertexList[b].V[1],Sb[3]);
        mpz_sub(temp3,temp1,temp2);
        mpz_mul(temp1,vertexList[c].V[1],Sc[3]);
        mpz_mul(temp2,vertexList[d].V[1],Sd[3]);
        mpz_sub(temp1,temp1,temp2);
        mpz_add(Det4,temp1,temp3);
/*
        Now compute all minors:
        Dmnp = Minor(m,n,p,1,2,3) = Det | m(1) m(2) m(3) |
                                        | n(1) n(2) n(3) |
                                        | p(1) p(2) p(3) |
*/
        mpz_mul(temp1,vertexList[a].V[1],Sbc[3]);
        mpz_mul(temp2,vertexList[b].V[1],Sac[3]);
        mpz_sub(temp3,temp1,temp2);
        mpz_mul(temp1,vertexList[c].V[1],Sab[3]);
        mpz_add(Dabc,temp3,temp1);

        mpz_mul(temp1,vertexList[a].V[1],Sbd[3]);
        mpz_mul(temp2,vertexList[b].V[1],Sad[3]);
        mpz_sub(temp3,temp1,temp2);
        mpz_mul(temp1,vertexList[d].V[1],Sab[3]);
        mpz_add(Dabd,temp3,temp1);

        mpz_mul(temp1,vertexList[a].V[1],Scd[3]);
        mpz_mul(temp2,vertexList[c].V[1],Sad[3]);
        mpz_sub(temp3,temp1,temp2);
        mpz_mul(temp1,vertexList[d].V[1],Sac[3]);
        mpz_add(Dacd,temp3,temp1);

        mpz_mul(temp1,vertexList[b].V[1],Scd[3]);
        mpz_mul(temp2,vertexList[c].V[1],Sbd[3]);
        mpz_sub(temp3,temp1,temp2);
        mpz_mul(temp1,vertexList[d].V[1],Sbc[3]);
        mpz_add(Dbcd,temp3,temp1);
/*
        We also need :
        Det =Det| m(1) m(2) m(3) m(4) |
                | n(1) n(2) n(3) n(4) |
                | p(1) p(2) p(3) p(4) |
                | q(1) q(2) q(3) q(4) |
*/
        mpz_mul(temp1,vertexList[a].V[4],Dbcd);
        mpz_mul(temp2,vertexList[b].V[4],Dacd);
        mpz_sub(temp3,temp2,temp1);
        mpz_mul(temp1,vertexList[c].V[4],Dabd);
        mpz_mul(temp2,vertexList[d].V[4],Dabc);
        mpz_sub(temp1,temp2,temp1);
        mpz_add(Dabcd,temp3,temp1);
/*
        The radius of the circumsphere of the weighted tetrahedron is then:
        r_t = (Det1*Det1 + Det2*Det2 + Det3*Det3 + 4*Det4*Dabcd)/(4*Det4*Det4)
*/

        mpz_mul(temp1,Det4,Det4); coef=4; mpz_mul_si(den,temp1,coef);

        mpz_mul(temp1,Det1,Det1); mpz_mul(temp2,Det2,Det2);
        mpz_add(temp1,temp1,temp2); mpz_mul(temp2,Det3,Det3);
        mpz_add(temp1,temp1,temp2); mpz_mul(temp2,Det4,Dabcd);
        mpz_mul_si(temp2,temp2,coef); mpz_add(num,temp1,temp2);
/*
        set p and q to num and den respctively
*/
        mpz_set(*p,num); mpz_set(*q,den);
/*
        sorting is important :D
*/
        ISort3(&a,&b,&c);
        ISort4(&a,&b,&c,&d);
/*
        Get triangles (ijk), (ijl), (ikl) and (jkl)
        (in that order)

        triangle(1) = tetra_link(4,index)
        triangle(2) = tetra_link(3,index)
        triangle(3) = tetra_link(2,index)
        triangle(4) = tetra_link(1,index)

        Get edges (ij),(ik),(il),(jk),(jl),(kl)
        (in that order)

        pair(1) = trig_link(3,triangle(1))
        pair(2) = trig_link(2,triangle(1))
        pair(3) = trig_link(2,triangle(2))
        pair(4) = trig_link(1,triangle(1))
        pair(5) = trig_link(1,triangle(2))
        pair(6) = trig_link(1,triangle(3))
*/
//        Vertex* vA = &vertexList[a];
//        Vertex* vB = &vertexList[b];
//        Vertex* vC = &vertexList[c];
//        double* cA = vA->Coordinates;
//        double* cB = vB->Coordinates;
//        double* cC = vC->Coordinates;
        bool isattached = false;
        tri[0] = delcx->DeluanayTet[index].TetLink[4];                   //abc wrt d
        if (!delcx->DeluanayTrigs[tri[0]].IsAttached)
        {
                isattached = false;
                CheckTriangle(vertexList[a].V, vertexList[b].V, vertexList[c].V, vertexList[d].V, &isattached);
                delcx->DeluanayTrigs[tri[0]].IsAttached = isattached;
        }

        tri[1] = delcx->DeluanayTet[index].TetLink[3];           	//abd wrt c
        if (!delcx->DeluanayTrigs[tri[1]].IsAttached)
        {
                isattached = false;
                CheckTriangle(vertexList[a].V, vertexList[b].V, vertexList[d].V, vertexList[c].V, &isattached);
                delcx->DeluanayTrigs[tri[1]].IsAttached = isattached;
        }

        tri[2] = delcx->DeluanayTet[index].TetLink[2];                   //acd wrt b
        if (!delcx->DeluanayTrigs[tri[2]].IsAttached)
        {
                isattached = false;
                CheckTriangle(vertexList[a].V, vertexList[c].V, vertexList[d].V, vertexList[b].V, &isattached);
                delcx->DeluanayTrigs[tri[2]].IsAttached = isattached;
        }

        tri[3] = delcx->DeluanayTet[index].TetLink[1];                   //bcd wrt a
        if (!delcx->DeluanayTrigs[tri[3]].IsAttached)
        {
                isattached = false;
                CheckTriangle(vertexList[b].V, vertexList[c].V, vertexList[d].V, vertexList[a].V, &isattached);
                delcx->DeluanayTrigs[tri[3]].IsAttached = isattached;
        }

        pair[0] = delcx->DeluanayTrigs[tri[0]].TrigLink[3];
        if (!delcx->DeluanayEdges[pair[0]].IsAttached)                       //ab wrt c
        {
                isattached = false;
                CheckEdge(vertexList[a].V, vertexList[b].V, vertexList[c].V, &isattached);
                delcx->DeluanayEdges[pair[0]].IsAttached = isattached;
        }
        if (!delcx->DeluanayEdges[pair[0]].IsAttached)                       //ab wrt d
        {
                isattached = false;
                CheckEdge(vertexList[a].V, vertexList[b].V, vertexList[d].V, &isattached);
                delcx->DeluanayEdges[pair[0]].IsAttached = isattached;
        }

        pair[1] = delcx->DeluanayTrigs[tri[0]].TrigLink[2];
        if (!delcx->DeluanayEdges[pair[1]].IsAttached)               		//ac wrt b
        {
                isattached = false;
                CheckEdge(vertexList[a].V, vertexList[c].V, vertexList[b].V, &isattached);
                delcx->DeluanayEdges[pair[1]].IsAttached = isattached;
        }
        if (!delcx->DeluanayEdges[pair[1]].IsAttached)               		//ac wrt d
        {
                isattached = false;
                CheckEdge(vertexList[a].V, vertexList[c].V, vertexList[d].V, &isattached);
                delcx->DeluanayEdges[pair[1]].IsAttached = isattached;
        }

        pair[2] = delcx->DeluanayTrigs[tri[1]].TrigLink[2];
        if (!delcx->DeluanayEdges[pair[2]].IsAttached)		                //ad wrt b
        {
                isattached = false;
                CheckEdge(vertexList[a].V, vertexList[d].V, vertexList[b].V, &isattached);
                delcx->DeluanayEdges[pair[2]].IsAttached = isattached;
        }
        if (!delcx->DeluanayEdges[pair[2]].IsAttached)		                //ad wrt c
        {
                isattached = false;
                CheckEdge(vertexList[a].V, vertexList[d].V, vertexList[c].V, &isattached);
                delcx->DeluanayEdges[pair[2]].IsAttached = isattached;
        }

        pair[3] = delcx->DeluanayTrigs[tri[0]].TrigLink[1];
        if (!delcx->DeluanayEdges[pair[3]].IsAttached)		                //bc wrt a
        {
                isattached = false;
                CheckEdge(vertexList[b].V, vertexList[c].V, vertexList[a].V, &isattached);
                delcx->DeluanayEdges[pair[3]].IsAttached = isattached;
        }
        if (!delcx->DeluanayEdges[pair[3]].IsAttached)		                //bc wrt d
        {
                isattached = false;
                CheckEdge(vertexList[b].V, vertexList[c].V, vertexList[d].V, &isattached);
                delcx->DeluanayEdges[pair[3]].IsAttached = isattached;
        }

        pair[4] = delcx->DeluanayTrigs[tri[1]].TrigLink[1];
        if (!delcx->DeluanayEdges[pair[4]].IsAttached)		                //bd wrt a
        {
                isattached = false;
                CheckEdge(vertexList[b].V, vertexList[d].V, vertexList[a].V, &isattached);
                delcx->DeluanayEdges[pair[4]].IsAttached = isattached;
        }
        if (!delcx->DeluanayEdges[pair[4]].IsAttached)		                //bc wrt c
        {
                isattached = false;
                CheckEdge(vertexList[b].V, vertexList[d].V, vertexList[c].V, &isattached);
                delcx->DeluanayEdges[pair[4]].IsAttached = isattached;
        }

        pair[5] = delcx->DeluanayTrigs[tri[2]].TrigLink[1];
        if (!delcx->DeluanayEdges[pair[5]].IsAttached)		                //cd wrt b
        {
                isattached = false;
                CheckEdge(vertexList[c].V, vertexList[d].V, vertexList[b].V, &isattached);
                delcx->DeluanayEdges[pair[5]].IsAttached = isattached;
        }
        if (!delcx->DeluanayEdges[pair[5]].IsAttached)		                //cd wrt a
        {
                isattached = false;
                CheckEdge(vertexList[c].V, vertexList[d].V, vertexList[a].V, &isattached);
                delcx->DeluanayEdges[pair[5]].IsAttached = isattached;
        }
/*
        Clear local GMP variables
*/
        mpz_clear(val1); mpz_clear(val2); mpz_clear(val3); mpz_clear(val4);

        mpz_clear (num); mpz_clear (den);

        mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);

        for (i=0; i < 4; i++)
        {
                mpz_clear(Sab[i]); mpz_clear(Sac[i]); mpz_clear(Sad[i]);
                mpz_clear(Sbc[i]); mpz_clear(Sbd[i]); mpz_clear(Scd[i]);
                mpz_clear(Dab[i]); mpz_clear(Dac[i]); mpz_clear(Dad[i]);
                mpz_clear(Dbc[i]); mpz_clear(Dbd[i]); mpz_clear(Dcd[i]);
                mpz_clear(Sa[i]); mpz_clear(Sb[i]); mpz_clear(Sc[i]);
                mpz_clear(Sd[i]);
                mpz_clear(Sam1[i]); mpz_clear(Sbm1[i]); mpz_clear(Scm1[i]);
                mpz_clear(Sdm1[i]);
                mpz_clear(Ta[i]); mpz_clear(Tb[i]); mpz_clear(Tc[i]);
                mpz_clear(Td[i]);
                mpz_clear(Tam1[i]); mpz_clear(Tbm1[i]); mpz_clear(Tcm1[i]);
                mpz_clear(Tdm1[i]);
                mpz_clear(Ua[i]); mpz_clear(Ub[i]); mpz_clear(Uc[i]);
                mpz_clear(Ud[i]);
                mpz_clear(Deter[i]);
        }

        mpz_clear(Dabc); mpz_clear(Dabd); mpz_clear(Dacd); mpz_clear(Dbcd);
        mpz_clear(Det1); mpz_clear(Det2); mpz_clear(Det3); mpz_clear(Dabcd);
        mpz_clear(Det4);
}
