/*
 *  Eigenfunction.h
 *
 *  Wormulator
 *  Copyright (C) 2011 Brian C. Ross
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or any later version.

 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef AnalyticsHeader
#define AnalyticsHeader

#include "LinMath.h"
#include "userfn.h"

typedef struct {
    ccInt variables_num;
    ccInt terms_num;
} PolyType;

typedef struct {
    double Coef;
    ccInt p_power;
    ccInt Ki_power;     // for alignment, there should be an even number of terms (double must be on even long boundary)
} OnePolyTerm;

typedef struct {
    PolyType Header;
    OnePolyTerm One;
} OneTermPoly;

typedef struct {
    PolyType Header;
    OnePolyTerm One;
    OnePolyTerm Two;
} TwoTermPoly;

typedef struct {
    PolyType Header;
    OnePolyTerm One;
    OnePolyTerm Two;
    OnePolyTerm Three;
} ThreeTermPoly;

typedef struct {
    PolyType Header;
    OnePolyTerm One;
    OnePolyTerm Two;
    OnePolyTerm Three;
    OnePolyTerm Four;
} FourTermPoly;

typedef struct {
    cdouble root;
    cdouble uncertainty;
    int multiplicity;
    ccBool finished;           // used by the root-solver
} ComplexRootType;

extern int SetKVals(int argc, char **argv);
extern int call_MakeTwoPolyTables(int, char **);
extern int MakeTwoPolyTables(int);
extern int AddPolyProduct(linkedlist *, linkedlist *, linkedlist *, linkedlist *, double *);
extern int GetWFactors(int, char **);
extern int call_CalcG(int, char **);
extern int CalcG(int, vector_3, vector_3, vector_3, vector_3, vector_3, double, int, int, double, double, unsigned char, ccBool, ccBool, linkedlist[2], cdouble *, cdouble *, cdouble *);
extern int Gllmj(linkedlist[2],int, int, int, int, double, int, double, cdouble *);
extern double AngleBetween(vector_3, vector_3, vector_3);
extern int Renorm_3v(vector_3 *);
extern void GramSchmidt(vector_3 *, vector_3 *);
extern int InverseLaplace(cdouble *, int, cdouble *, int, double, cdouble *);
extern int AddToRootList(linkedlist *, linkedlist *, linkedlist *, linkedlist *, cdouble *, int, int);
extern int call_Wigner(int, char **);
extern double Wigner(int, int, int, double);
extern void GetAngles(vector_3 *, vector_3 *, vector_3 *, vector_3 *, vector_3 *, double *, double *);
extern void GetOneAngle(vector_3 *, vector_3 *, vector_3 *, vector_3 *, double *);
extern inline void fillVector(vector_3 *, double, double, double);
extern void RotateVector(vector_3 *, vector_3 *, vector_3 *, vector_3 *);
extern void p_Poly_Product(double *, double *, double *, int, int);
extern inline void Make_p_Poly(double *, double *, double);
extern inline cdouble Eval_p_Poly(double *, cdouble, int);

extern int CopyPoly(linkedlist *, linkedlist *);
extern int AddPoly(linkedlist *, linkedlist *);
extern int MultiplyPoly(linkedlist *, linkedlist *, linkedlist *);
extern void ScalarMultiplyPoly(linkedlist *, double *);
extern int ReducePoly_iK(linkedlist *, linkedlist *, int, double);
extern void iReducePoly(linkedlist *, int);
extern void ComplexConjugatePoly(linkedlist *, int);
extern cdouble EvaluatePoly(linkedlist *, cdouble *, int *, int);
extern void SimplifyPoly(linkedlist *);
extern int ResizePoly(linkedlist *, int, int);
extern int InitPoly(linkedlist *, int, int, char *);
extern inline int PolyElementIndex(linkedlist *, int);
extern double *PolyElement(linkedlist *, int);
extern inline int PolyVariables(linkedlist *);
extern int PolyTerms(linkedlist *);

extern int *PolyCoef(double *, int);
extern inline linkedlist *NumPolyCoef(linkedlist *, int, int, int);
extern inline linkedlist *DenomPolyCoef(linkedlist *, int, int, int);

extern double IntegerPow(double, int);
extern double LogFactorial(int);
extern double Factorial(int);

extern double K_max, K_step;
extern linkedlist *Wlm, *w_plus, *w_minus, *LogFactorialList, *FactorialList;
extern linkedlist *Wlm_roots, *w_plus_roots, *w_minus_roots;

#endif
