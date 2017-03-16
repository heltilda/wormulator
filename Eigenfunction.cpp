/*
 *  Eigenfunction.cpp - calculates polymer end statistics using the method of Spakowitz (2005)
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

#include <math.h>
#include <iostream>
//#include <stdio.h>
#include "LinMath.h"
//#include "extrnl.h"
#include "lnklst.h"
#include "intrpt.h"
#include "Eigenfunction.h"
#include "cpoly.h"
#include "Interpolation.h"

// a, b are alpha and beta in Spakowitz's paper
#define a_lmj_sq(l, m, j) (((double) (((l)-(m))*((l)+(m))*((l)-(j))*((l)+(j))))/((l)*(l)*(4*(l)*(l)-1)))
#define b_lmj(l, m, j) (((double) ((m)*(j)))/((l)*((l)+1)))

// TriangleSum() computes x_min + (x_min+1) + ... + x_max
#define TriangleSum(x_min, x_max) (((x_max)-(x_min)+1)*((x_max)+(x_min))/2)

double K_max, K_step;
linkedlist *Wlm, *w_plus, *w_minus;
linkedlist *Wlm_roots, *w_plus_roots, *w_minus_roots, *LogFactorialList, *FactorialList;



// Sets the global variables K_max and K_step

int SetKVals(int argc, char **argv)
{
    arg_info *ArgInfo = (arg_info *) *(argv+argc);
    
    const int ArgTypes[] = { double_type, double_type };
    const int ArgIndices[] = { 1, 1 };
    
    if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "SetKVals") != passed)  return 1;
    
    getArgs(argc, argv, byValue(&K_step), byValue(&K_max));
    
    if ((argc != 2) || (ArgInfo->argType != double_type) || ((ArgInfo+1)->argType < double_type)
                || (ArgInfo->argIndices != 1) || ((ArgInfo+1)->argIndices != 1))
        {  printf("usage: call(\"SetKVals\", (doubles) K_step, K_max)\n");  return 1;  }
    
    return passed;
}



// Next two routines:  initialize the continued fractions

int call_MakeTwoPolyTables(int argc, char **argv)
{
    int l_max, c1, WlmIndices, mtpt_retrn;
    arg_info *ArgInfo = (arg_info *) *(argv+argc);
    linkedlist *loopLL;
    
    const int ArgTypes[] = { string_type };
    const int ArgIndices[] = { -1 };
    
    if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "MakeTwoPolyTables") != passed)  return 1;
    
    
        // l_max wasn't passed explicitly; compute it along with number of polynomials in each w-, w+, W
    
    l_max = 0;
    WlmIndices = 2;
    while (ArgInfo->argIndices > 3*WlmIndices)     {
        l_max++;
        WlmIndices += 2*TriangleSum(1, l_max+1);        // factor of 2 for num/denom
        if (ArgInfo->argIndices < 3*WlmIndices)        // factor of 3 for w-, w+, W
            {  printf("usage: call(\"MakeTwoPolyTables\", string_array { Wlm, W_plus, W_minus })\n");  return 1;  }
    };
    
    
        // init num/denom of each fraction to a single polynomial in two variables (K and i) having zero terms
    
    for (c1 = 0; c1 < ArgInfo->argIndices; c1++)   {
        loopLL = ((linkedlist *) *argv) + c1;
        if ( loopLL->memory != 0 )  deleteLinkedList(loopLL);
        mtpt_retrn = newLinkedList(loopLL, sizeof(PolyType), 1, 0, ccTrue);
        if (mtpt_retrn != passed)  return 10;
        if (c1 < 3*WlmIndices)  ((PolyType *) element(loopLL, 1))->variables_num = 2;
        else  ((PolyType *) element(loopLL, 1))->variables_num = 2;
    }
    
    Wlm = (linkedlist *) argv[0];
    w_plus = Wlm + WlmIndices;
    w_minus = w_plus + WlmIndices;
    
    mtpt_retrn = MakeTwoPolyTables(l_max);
    if (mtpt_retrn != passed)  return 10;
    
    for (c1 = 0; c1 < ArgInfo->argIndices; c1++)   {
        SimplifyPoly(((linkedlist *) *argv) + c1);
    }
    
    return passed;
}

int MakeTwoPolyTables(int l_max)
{
    linkedlist *w_plus_num, *w_plus_denom, *w_minus_num, *w_minus_denom;
    linkedlist w_denominator, KsqSubPoly, ZeroPoly, OnePoly, ConstPoly, tempPoly;
    int poly_retrn, l, m, j, n_l, l_cutoff;
    
    
        // templates for the nums/denoms explicitly written in the Spakowitz paper
    
    const OneTermPoly Zero = { { 2, 1 }, { 0, 0, 0 } };
    const OneTermPoly One = { { 2, 1 }, { 1, 0, 0 } };
    const ThreeTermPoly w_denom_init = { { 2, 3 }, { 1, 1, 0 }, { 0, 0, 0 }, { 0, 0, 1 } };
    const FourTermPoly w_denom_init_plus_Ksq = { { 2, 4 }, { 1, 1, 0 }, { 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 2 } };
    const OneTermPoly Ksp_init = { { 2, 1 }, { 1, 0, 2 } };
    const OnePolyTerm p_offset = { 1, 1, 0 };
    
    
        // copy the 'const' template polynomials over to the part of memory that we are supposed to be writing to
    
    w_denominator.memory = 0;
    KsqSubPoly.memory = 0;
    ZeroPoly.memory = 0;
    OnePoly.memory = 0;
    tempPoly.memory = 0;
    
    poly_retrn = newLinkedList(&w_denominator, sizeof(ThreeTermPoly), 1, 0, ccTrue);
    if (poly_retrn == passed)  poly_retrn = newLinkedList(&KsqSubPoly, sizeof(OneTermPoly), 1, 0, ccTrue);
    if (poly_retrn == passed)  poly_retrn = newLinkedList(&ZeroPoly, sizeof(OneTermPoly), 1, 0, ccTrue);
    if (poly_retrn == passed)  poly_retrn = newLinkedList(&OnePoly, sizeof(OneTermPoly), 1, 0, ccTrue);
    if (poly_retrn == passed)  poly_retrn = newLinkedList(&ConstPoly, sizeof(OneTermPoly), 1, 0, ccTrue);
    if (poly_retrn == passed)  poly_retrn = newLinkedList(&tempPoly, sizeof(OneTermPoly), 1, 0, ccTrue);
    if (poly_retrn != passed)  return poly_retrn;
    
    setElements(&w_denominator, 1, w_denominator.elementNum, (void *) &w_denom_init);
    setElements(&KsqSubPoly, 1, KsqSubPoly.elementNum, (void *) &Ksp_init);
    setElements(&ZeroPoly, 1, ZeroPoly.elementNum, (void *) &Zero);
    setElements(&OnePoly, 1, OnePoly.elementNum, (void *) &One);
    setElements(&ConstPoly, 1, OnePoly.elementNum, (void *) &One);
    
    
        // now fill the continued fractions
    
    for (m = 0; m <= l_max; m++)  {
    for (j = 0; j <= m; j++)    {
        l_cutoff = m;
        
        
            // w- polynomials
        
        poly_retrn = InitPoly(NumPolyCoef(w_minus, l_cutoff, m, j), 2, 1, (char *) &(One.One));
        if (poly_retrn == passed)  poly_retrn = InitPoly(DenomPolyCoef(w_minus, l_cutoff, m, j), 2, 3, (char *) &(w_denom_init.One));
        if (poly_retrn != passed)  return poly_retrn;
        
            // for numerical reasons, subtract l_min^2 from all Laplace variables
        
        *PolyElement(DenomPolyCoef(w_minus, l_cutoff, m, j), 2) = l_cutoff*(l_cutoff+1) - l_cutoff*l_cutoff;
        if (l_cutoff != 0)  *PolyElement(DenomPolyCoef(w_minus, l_cutoff, m, j), 3) = -b_lmj(l_cutoff, m, j);
        else  *PolyElement(DenomPolyCoef(w_minus, l_cutoff, m, j), 3) = 0;
        
        n_l = 1;                    // order+1 of the preceding l-1 polynomial
        for (l = l_cutoff+1; l <= l_max; l++)   {
            *PolyElement(&KsqSubPoly, 1) = -a_lmj_sq(l, m, j);
            *PolyElement(&w_denominator, 2) = l*(l+1) - l_cutoff*l_cutoff;
            *PolyElement(&w_denominator, 3) = -b_lmj(l, m, j);
            
            poly_retrn = CopyPoly(DenomPolyCoef(w_minus, l-1, m, j), NumPolyCoef(w_minus, l, m, j));
            if (poly_retrn == passed)  poly_retrn = MultiplyPoly(DenomPolyCoef(w_minus, l-1, m, j), &w_denominator, DenomPolyCoef(w_minus, l, m, j));
            if (poly_retrn == passed)  poly_retrn = AddPolyProduct(DenomPolyCoef(w_minus, l, m, j), NumPolyCoef(w_minus, l-1, m, j), &KsqSubPoly, &tempPoly, 0);
            if (poly_retrn != passed)  return poly_retrn;
            
            SimplifyPoly(NumPolyCoef(w_minus, l, m, j));
            SimplifyPoly(DenomPolyCoef(w_minus, l, m, j));
            n_l++;
        }
        
        
            // w+ polynomials
        
        poly_retrn = InitPoly(NumPolyCoef(w_plus, l_max, m, j), 2, 1, (char *) &(One.One));
        if (poly_retrn == passed)  poly_retrn = InitPoly(DenomPolyCoef(w_plus, l_max, m, j), 2, 4, (char *) &(w_denom_init_plus_Ksq.One));
        if (poly_retrn != passed)  return poly_retrn;
        
         *PolyElement(DenomPolyCoef(w_plus, l_max, m, j), 2) = l_max*(l_max+1) - l_cutoff*l_cutoff;
        if (l_max != 0)  *PolyElement(DenomPolyCoef(w_plus, l_max, m, j), 3) = -b_lmj(l_max, m, j);
        else  *PolyElement(DenomPolyCoef(w_plus, l_max, m, j), 3) = 0;
        *PolyElement(DenomPolyCoef(w_plus, l_max, m, j), 4) = -a_lmj_sq(l_max+1, m, j)/((l_max+1)*(l_max+2));
        
        n_l = 1;        // order+1 of the preceding l+1 polynomial
        for (l = l_max-1; l >= l_cutoff; l--)   {
            *PolyElement(&KsqSubPoly, 1) = -a_lmj_sq(l+1, m, j);
            *PolyElement(&w_denominator, 2) = l*(l+1) - l_cutoff*l_cutoff;
            if (l != 0)  *PolyElement(&w_denominator, 3) = -b_lmj(l, m, j);
            else  *PolyElement(&w_denominator, 3) = 0;
            
            poly_retrn = CopyPoly(DenomPolyCoef(w_plus, l+1, m, j), NumPolyCoef(w_plus, l, m, j));
            if (poly_retrn == passed)  poly_retrn = MultiplyPoly(DenomPolyCoef(w_plus, l+1, m, j), &w_denominator, DenomPolyCoef(w_plus, l, m, j));
            if (poly_retrn == passed)  poly_retrn = AddPolyProduct(DenomPolyCoef(w_plus, l, m, j), NumPolyCoef(w_plus, l+1, m, j), &KsqSubPoly, &tempPoly, 0);
            if (poly_retrn != passed)  return poly_retrn;
            
            SimplifyPoly(NumPolyCoef(w_plus, l, m, j));
            SimplifyPoly(DenomPolyCoef(w_plus, l, m, j));
            n_l++;
        }
        
        
            // W polynomials
        
        for (l = l_cutoff; l <= l_max; l++) {
            if (l == l_max)     {
                w_plus_num = &ConstPoly;
                *PolyElement(&ConstPoly, 1) = 1./((l+1)*(l+2));
                w_plus_denom = &OnePoly;            }
            else        {
                w_plus_num = NumPolyCoef(w_plus, l+1, m, j);
                w_plus_denom = DenomPolyCoef(w_plus, l+1, m, j);        }
            
            if (l == l_cutoff)      {
                w_minus_num = &ZeroPoly;
                w_minus_denom = &OnePoly;       }
            else        {
                w_minus_num = NumPolyCoef(w_minus, l-1, m, j);
                w_minus_denom = DenomPolyCoef(w_minus, l-1, m, j);      }
            
            AddPolyProduct(NumPolyCoef(Wlm, l, m, j), w_minus_denom, w_plus_denom, &tempPoly, 0);
            AddPolyProduct(DenomPolyCoef(Wlm, l, m, j), w_minus_denom, w_plus_denom, &tempPoly, (double *) &p_offset);
            *PolyElement(&w_denominator, 2) = l*(l+1) - l_cutoff*l_cutoff;
            AddPolyProduct(DenomPolyCoef(Wlm, l, m, j), w_minus_denom, w_plus_denom, &tempPoly, PolyElement(&w_denominator, 2));
            if (l != 0)  *PolyElement(&w_denominator, 3) = -b_lmj(l, m, j);
            else  *PolyElement(&w_denominator, 3) = 0;
            AddPolyProduct(DenomPolyCoef(Wlm, l, m, j), w_minus_denom, w_plus_denom, &tempPoly, PolyElement(&w_denominator, 3));
            if (l != l_cutoff)      {
                *PolyElement(&KsqSubPoly, 1) = -a_lmj_sq(l, m, j);
                AddPolyProduct(DenomPolyCoef(Wlm, l, m, j), w_minus_num, w_plus_denom, &tempPoly, PolyElement(&KsqSubPoly, 1));     }
            *PolyElement(&KsqSubPoly, 1) = -a_lmj_sq(l+1, m, j);
            AddPolyProduct(DenomPolyCoef(Wlm, l, m, j), w_minus_denom, w_plus_num, &tempPoly, PolyElement(&KsqSubPoly, 1));
            *PolyElement(&KsqSubPoly, 1) = 1;
            
            SimplifyPoly(NumPolyCoef(Wlm, l, m, j));
            SimplifyPoly(DenomPolyCoef(Wlm, l, m, j));
        }
    }}
    
    deleteLinkedList(&w_denominator);
    deleteLinkedList(&KsqSubPoly);
    deleteLinkedList(&ZeroPoly);
    deleteLinkedList(&OnePoly);
    deleteLinkedList(&tempPoly);
    
    return passed;
}


// Adds multiplier * Poly1 * Poly2 to DestPoly

int AddPolyProduct(linkedlist *DestPoly, linkedlist *Poly1, linkedlist *Poly2, linkedlist *tempPoly, double *multiplier)
{
    int linkedLretrn;
    
    linkedLretrn = MultiplyPoly(Poly1, Poly2, tempPoly);
    if ((linkedLretrn == passed) && (multiplier != 0))  ScalarMultiplyPoly(tempPoly, multiplier);
    if (linkedLretrn == passed)  linkedLretrn = AddPoly(DestPoly, tempPoly);
    return linkedLretrn;
}


// Finds the roots of the w-/w+/W polynomials (numerator & denominator)

int GetWFactors(int argc, char **argv)
{
    int c1, c2, c3, NumRoots, NumberOfKs, PolyLength, EP_no_derivs[2] = { 0, 0 }, EP_deriv[2] = { 1, 0 };
    linkedlist *W_LL, *Factor_LL, Reduced_LL, RootsLL, WorkspaceLL;
    double *loopTerm, loopK, *DTopElement;
    cdouble EP_inputs[2] = { c_zero, c_zero };
    int FoundRoots;
    cdouble *TopElement;
    arg_info *ArgInfo = (arg_info *) *(argv+argc);
    int linkedLretrn;
    ccBool NotAllRootsFound = ccFalse;
    
    const int ArgTypes[] = { string_type, string_type, double_type, double_type };
    const int ArgIndices[] = { -1, -1, 1, 1 };
    
    if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "GetWFactors") != passed)  return 1;
    
    getArgs(argc, argv, fromArg(2), byValue(&K_step), byValue(&K_max));
    
    if (K_max < 0)  {  printf("GetWFactors error: K_max must be greater than 0\n");  return 3;  }
    NumberOfKs = (int) floor(0.5 + K_max/K_step);
    
    RootsLL.memory = 0;
    linkedLretrn = newLinkedList(&RootsLL, 0, sizeof(ComplexRootType), 255, ccFalse);
    if (linkedLretrn != passed)  return 1;
    
    WorkspaceLL.memory = 0;
    linkedLretrn = newLinkedList(&WorkspaceLL, 0, 1, 255, ccFalse);
    if (linkedLretrn != passed)  return 1;
    
    Reduced_LL.memory = 0;
    linkedLretrn = newLinkedList(&Reduced_LL, 0, 1, 255, ccFalse);
    if (linkedLretrn != passed)  return 1;
    
    for (c1 = 0; c1 < (ArgInfo+1)->argIndices; c1++)       {
        W_LL = ((linkedlist *) (*argv)) + c1;
        Factor_LL = ((linkedlist *) (*(argv+1))) + c1;
        
        
            // calculate number of roots based on max power in variable 1 (iK), and allocate storage
        
        NumRoots = 0;
        for (c2 = 1; c2 <= PolyTerms(W_LL); c2++)   {
            loopTerm = PolyElement(W_LL, c2);
            if (*PolyCoef(loopTerm, 1) > NumRoots)  NumRoots = *PolyCoef(loopTerm, 1);
        }
        
        if (Factor_LL->memory != 0)  deleteLinkedList(Factor_LL);
        linkedLretrn = newLinkedList(Factor_LL, NumberOfKs*(2*NumRoots+2)*sizeof(double), 1, 0, ccFalse);
        if (linkedLretrn != passed)  return 1;
        
        
            // call the root solver repeatedly for each K value (each round extracts all roots)
        
        for (c2 = 0; c2 < NumberOfKs; c2++)     {
            
            loopK = (1.*c2+0.5)*K_step;
            ReducePoly_iK(W_LL, &Reduced_LL, 2, loopK);
            
            PolyLength = NumRoots+1;
            
            deleteElements(&WorkspaceLL, 1, WorkspaceLL.elementNum);
            linkedLretrn = addElements(&WorkspaceLL, (PolyLength + NumRoots)*sizeof(cdouble), ccTrue);
            if (linkedLretrn == passed)  linkedLretrn = defragmentLinkedList(&WorkspaceLL);
            if (linkedLretrn == passed)  linkedLretrn = defragmentLinkedList(&RootsLL);
            if (linkedLretrn != passed)  return 1;
            TopElement = (cdouble *) element(&WorkspaceLL, 1);
            DTopElement = (double *) element(&WorkspaceLL, 1);
            
            
                // convert the polynomial into a list of complex coefficients, from p^N to p^0
            
            /*for (c3 = 1; c3 <= PolyTerms(&Reduced_LL); c3++)      {
                loopTerm = PolyElement(&Reduced_LL, c3);
                if (PolyVariables(&Reduced_LL) == 2)        {
                    if (*PolyCoef(loopTerm, 2) == 0)  (TopElement+NumRoots-(*PolyCoef(loopTerm, 1)))->real += *loopTerm;
                    else  (TopElement+NumRoots-(*PolyCoef(loopTerm, 1)))->imag += *loopTerm;        }
                else  (TopElement+NumRoots-(*PolyCoef(loopTerm, 1)))->real += *loopTerm;            }
            
            
            
                // call the root solver (in PolyRoots.cpp)
            
            *(cdouble *) element(Factor_LL, c2*(2*NumRoots+2)*sizeof(double)+1) = *TopElement;
            if (NumRoots >= 1)      {
                
                Gradient_ComplexPolyRoots(TopElement, &PolyLength, (cdouble *) element(Factor_LL, (c2*(2*NumRoots+2)+2)*sizeof(double)+1), 1.e-15);
                element(Factor_LL, (c2+1)*(2*NumRoots+1)*sizeof(double));      // what was this supposed to be doing?
                if (PolyLength > 1)  NotAllRootsFound = ccTrue;
                
                {
                    double *tst_Dtop = (double *) TopElement;
                    double *tst = (double *) element(Factor_LL, (c2*(2*NumRoots+2)+2)*sizeof(double)+1);
                    int cX;
                    
                    printf("(%i, 1)\n", PolyLength);
                    for (cX = (NumRoots+1)*2; cX < (2*NumRoots+1)*2; cX++)  printf("%g ", *(tst_Dtop+cX));
                    printf("\n");
                    for (cX = 0; cX < NumRoots*2; cX++)  printf("%g ", *(tst+cX));
                    printf("\n");
                }
            }*/
            //FClearElements(&WorkspaceLL, 1, WorkspaceLL.elementNum);
            //PolyLength = NumRoots+1;
            
            
                // convert the polynomial into a list of complex coefficients, from p^N to p^0
            
            for (c3 = 1; c3 <= PolyTerms(&Reduced_LL); c3++)        {
                loopTerm = PolyElement(&Reduced_LL, c3);
                if (PolyVariables(&Reduced_LL) == 2)        {
                    if (*PolyCoef(loopTerm, 2) == 0)  *(DTopElement+NumRoots-(*PolyCoef(loopTerm, 1))) += *loopTerm;
                    else  *(DTopElement+PolyLength+NumRoots-(*PolyCoef(loopTerm, 1))) += *loopTerm;     }
                else  *(DTopElement+NumRoots-(*PolyCoef(loopTerm, 1))) += *loopTerm;            }
            
            
                // call the root solver (in cpoly.cpp)
            
            *(double *) element(Factor_LL, c2*(2*NumRoots+2)*sizeof(double)+1) = *DTopElement;         // the 'multiplier' of the polynomial
            *(double *) element(Factor_LL, (c2*(2*NumRoots+2)+1)*sizeof(double)+1) = *(DTopElement + PolyLength);
            if (NumRoots >= 1)      {
                
                double *FirstRootTemp = DTopElement + 2*PolyLength;
                double *FirstRoot = LL_Double(Factor_LL, (c2*(2*NumRoots+2)+2)*sizeof(double)+1);
                int cR;
                
                FoundRoots = cpoly(DTopElement, DTopElement+PolyLength, (int) PolyLength-1, FirstRootTemp, FirstRootTemp+NumRoots);
                if (FoundRoots != NumRoots)  NotAllRootsFound = ccTrue;
                
                for (cR = 0; cR < NumRoots; cR++)       {           // polish roots using Newton's method
                    
                    cdouble poly_f, poly_df, one_root, last_root, droot;
                    double droot_norm, last_droot_norm;
                    
                    one_root.real = *(FirstRootTemp + cR);
                    one_root.imag = *(FirstRootTemp + NumRoots + cR);
                    
                    droot_norm = -1.;
                    do  {
                        last_root = one_root;
                        last_droot_norm = droot_norm;
                        
                        EP_inputs[0] = one_root;
                        poly_f = EvaluatePoly(&Reduced_LL, EP_inputs, EP_no_derivs, 2);
                        poly_df = EvaluatePoly(&Reduced_LL, EP_inputs, EP_deriv, 2);
                        if ((poly_df.real == 0.) && (poly_df.imag == 0.))  break;
                        
                        droot = cDiv(poly_f, poly_df);
                        droot_norm = droot.real*droot.real + droot.imag*droot.imag;
                        if ((last_droot_norm >= 0.) && (droot_norm >= last_droot_norm))  break;
                        
                        one_root = cSub(one_root, droot);
                    }  while ((last_root.real != one_root.real) || (last_root.imag != one_root.imag));
                    
                    *(FirstRoot + cR*2) = one_root.real;
                    *(FirstRoot + cR*2 + 1) = one_root.imag;     }
            }
            
            /*if (NumRoots >= 1)      {
                int cdbg_root;
                int cdr, cdi;
                double *FirstRoot = LL_Double(Factor_LL, (c2*(2*NumRoots+2)+2)*sizeof(double)+1);
                cdouble poly_f, one_root;
                
                for (cdbg_root = 0; cdbg_root < NumRoots; cdbg_root++)  {
                    
                    printf("\nRoot:  %g %g\n", *(FirstRoot+2*cdbg_root+0), *(FirstRoot+2*cdbg_root+1));
                    
                    for (cdr = -1; cdr <= 1; cdr++)  {
                    for (cdi = -1; cdi <= 1; cdi++)  {
                        
                        int tmp_int;
                        
                        poly_f = c_zero;
                        one_root.real = *(FirstRoot+2*cdbg_root+0);
                        one_root.imag = *(FirstRoot+2*cdbg_root+1);
                        
                        tmp_int = *(int *) &(one_root.real);
                        tmp_int += 5*cdr;
                        one_root.real = *(double *) &tmp_int;
                        
                        tmp_int = *(int *) &(one_root.imag);
                        tmp_int += 5*cdi;
                        one_root.imag = *(double *) &tmp_int;
                        
                        EP_inputs[0] = one_root;
                        poly_f = EvaluatePoly(&Reduced_LL, EP_inputs, EP_no_derivs, 2);
                        printf("(%g %g) + (%i %i):  %g %g\n", one_root.real, one_root.imag, cdr, cdi, poly_f.real, poly_f.imag);
                    }}
                    */
/*                    printf("\nRoot:  %g %g\n", *(FirstRoot+2*cdbg_root+0), *(FirstRoot+2*cdbg_root+1));
                    
                    for (cdr = -1; cdr <= 1; cdr++)  {
                    for (cdi = -1; cdi <= 1; cdi++)  {
                        
                        int tmp_int;
                        
                        poly_f = c_zero;
                        one_root.real = *(FirstRoot+2*cdbg_root+0);
                        one_root.imag = *(FirstRoot+2*cdbg_root+1);
                        
                        tmp_int = *(int *) &(one_root.real);
                        tmp_int += 5*cdr;
                        one_root.real = *(double *) &tmp_int;
                        
                        tmp_int = *(int *) &(one_root.imag);
                        tmp_int += 5*cdi;
                        one_root.imag = *(double *) &tmp_int;
                        
                        EP_inputs[0] = one_root;
                        poly_f = EvaluatePoly(&Reduced_LL, EP_inputs, EP_no_derivs, 2);
                        printf("(%g %g) + (%i %i):  %g %g\n", one_root.real, one_root.imag, cdr, cdi, poly_f.real, poly_f.imag);
                    }}
            }   }*/
        }
    }
    
    deleteLinkedList(&RootsLL);
    deleteLinkedList(&WorkspaceLL);
    deleteLinkedList(&Reduced_LL);
    
    if (NotAllRootsFound == ccTrue)  return 4;
    
    return passed;
}


// Next two routines:  calculate the Green's function

int call_CalcG(int argc, char **argv)
{
    double N, eta, twist;
    cdouble *cresult, *G_list, *G_K_space;
    linkedlist workspace[3];
    int l_max, top_l, c1, PredictedIndices, theta_step_num, phi_step_num, cg_retrn, sum_R;
    ccBool sum_tangent, sum_twist;
    vector_3 R, u0, n0, uf, nf;
    arg_info *ArgInfo = (arg_info *) *(argv+argc);
    int retrn;
    
    const int ArgTypes[] = { int_type, double_type, double_type, double_type, double_type, double_type, double_type, int_type, int_type,
                                                double_type, double_type, int_type, bool_type, bool_type, double_type, string_type, double_type, double_type,
                                                double_type, double_type };
    const int ArgIndices[] = { 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, -1, -2, -3, 1, 1 };
    
    if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "CalcG") != passed)  return 1;
    
    getArgs(argc, argv, byValue(&top_l), byValue(&R), byValue(&u0), byValue(&n0), byValue(&uf), byValue(&nf), byValue(&N),
                byValue(&theta_step_num), byValue(&phi_step_num), byValue(&eta), byValue(&twist), byValue(&sum_R), byValue(&sum_tangent), byValue(&sum_twist),
                &cresult, &Wlm_roots, &G_list, &G_K_space, byValue(&K_step), byValue(&K_max));
    
    if (sum_R != 2)      {
        if ((K_max <= 0) || (K_step <= 0))
                {  printf("CalcG() error: K_max or K_step not initialized correctly\n");  return 4;  }
        
        l_max = 0;
        PredictedIndices = 6;
        while (PredictedIndices < ArgInfo[15].argIndices)        {
            l_max++;
            PredictedIndices += 3*(l_max + 1)*(l_max + 2);
            if (PredictedIndices > ArgInfo[15].argIndices)  {  printf("CalcG() error: wrong number of strings\n");  return 2;  }
        }
        w_plus_roots = Wlm_roots + PredictedIndices/3;
        w_minus_roots = w_plus_roots + PredictedIndices/3;
        
        if (ArgInfo[16].argIndices < 2*(K_max/K_step+0.5))  {  printf("CalcG() error: wrong number of elements in G_list\n");  return 5;  }
        if (ArgInfo[17].argIndices < 2*(K_max/K_step+0.5)*theta_step_num*phi_step_num)
            {  printf("CalcG() error: wrong number of elements in G_K_space\n");  return 6;  }
        if (top_l > l_max)  {  printf("CalcG() error: top_l must be <= l_max from Init()\n");  return 7;  }
        if ((theta_step_num == 0) || (phi_step_num == 0))  {  printf("CalcG() error: theta/phi_step_num must be greater than zero\n");  return 4;  }        }
    
    for (c1 = 0; c1 < 2; c1++)      {
        workspace[c1].memory = 0;
        retrn = newLinkedList(workspace+c1, 0, sizeof(cdouble), 100, ccFalse);
        if (retrn != passed)  return retrn;         }
    
    FactorialList = &(workspace[2]);
    FactorialList->memory = 0;
    retrn = newLinkedList(FactorialList, 1, sizeof(double), 255, ccFalse);        // first element is cleared to zero
    if (retrn != passed)  return retrn;
    *LL_Double(FactorialList, 1) = 1;
    
    cg_retrn = CalcG(top_l, R, u0, n0, uf, nf, N, theta_step_num, phi_step_num, eta, twist, sum_R, sum_tangent, sum_twist, workspace, cresult, G_list, G_K_space);
    
    for (c1 = 0; c1 < 3; c1++)                      // the third is the factorial list
        deleteLinkedList(&(workspace[c1]));
    
    return cg_retrn;
}


int CalcG(int top_l, vector_3 R, vector_3 u0, vector_3 n0, vector_3 uf, vector_3 nf, double N, int theta_step_num, int phi_step_num, double eta, double twist,
            unsigned char sum_R, ccBool sum_tangent, ccBool sum_twist, linkedlist workspace[2], cdouble *result, cdouble *G_list, cdouble *G_K_space)
{
    int l0, lf, min_l, l_cutoff, m, j, K_counter, NumberOfKs, loop_m, loop_j, lm_step, lj_step, eval_m, eval_j, eval_temp;
    int lmj_it, loop_mj_iterations, theta_counter, phi_counter;
    double K_norm, theta, phi, delta_phi_f0, delta_psi_f0, W_theta_0, W_theta_f, real_factor, theta_step, unit_K_dot_R = 0., R_norm = 0., phi_step;
    cdouble Gllm_result, D_sum, exp_factor, inv_Fourier_factor, *one_K_element;
    vector_3 K;
    int retrn;
    ccBool Cyclization;
    const double one_over_8_pi_3 = 1./(8.*pi*pi*pi), one_over_16_pi_4 = 1./(16.*pi*pi*pi*pi);
    const double pi_4 = 4.*pi, pi_2 = 2.*pi;
    
    if (Renorm_3v(&u0) != passed)  {  printf("CalcG() error:  u0, n0 must be nonzero and non-parallel\n");  return 1;  }
    GramSchmidt(&n0, &u0);
    if (Renorm_3v(&n0) != passed)  {  printf("CalcG() error:  u0, n0 must be nonzero and non-parallel\n");  return 1;  }
    
    if (Renorm_3v(&uf) != passed)  {  printf("CalcG() error:  uf, nf must be nonzero and non-parallel\n");  return 1;  }
    GramSchmidt(&nf, &uf);
    if (Renorm_3v(&nf) != passed)  {  printf("CalcG() error:  uf, nf must be nonzero and non-parallel\n");  return 1;  }
    
    if ((dot(&R, &R) == 0) && (u0.x == uf.x) && (u0.y == uf.y) && (u0.z == uf.z) && (n0.x == nf.x) &&
                    (n0.y == nf.y) && (n0.z == nf.z) && (sum_R != 2) && (sum_twist == ccFalse) && (sum_tangent == ccFalse))  Cyclization = ccTrue;
    else  Cyclization = ccFalse;
    
    if (sum_tangent == ccTrue)  sum_twist = ccTrue;
    *result = c_zero;
    
    NumberOfKs = (int) floor(0.5 + K_max/K_step);
    theta_step = pi/theta_step_num;
    phi_step = 2*pi/phi_step_num;
    
    one_K_element = G_K_space;
    if ((Cyclization == ccFalse) && (sum_R != 2))      {
    for (theta_counter = 0; theta_counter < theta_step_num; theta_counter++)        {
    for (phi_counter = 0; phi_counter < phi_step_num; phi_counter++)      {
    for (K_counter = 0; K_counter < NumberOfKs; K_counter++)        {
        *one_K_element = c_zero;
        one_K_element++;
    }}}}
    
    for (l0 = 0; l0 <= top_l; l0++) {       // 2 lines below:  if sum_tangent (implies sum_twist) then l0 must equal lf for the integral
    for (lf = 0; lf <= top_l; lf++) {       // of G(l0,m,j)*G(lf,m,j) over the relative angle to be nonzero
    if ((l0 == lf) || ((Cyclization == ccFalse) && (sum_R != 2)))      {
    if ((lf == 0) || (sum_tangent == ccFalse))        {
        if (l0 < lf)  min_l = l0;
        else  min_l = lf;
        
        for (m = 0; m <= min_l; m++)        {
        if ((sum_tangent == ccFalse) || (m == 0))      {
        for (j = 0; j <= m; j++)        {               // only up to 'm' because we only computed the inv-Laplace transforms for |j| <= |m|
        if ((sum_twist == ccFalse) || (j == 0))        {
            l_cutoff = m;
            if ((sum_R != 2) && (Cyclization == ccFalse))     {
                
                
                    // compute the Laplace-transformed script-G in Spakowitz
                    // don't precompute this because we want to be able to change N w/o solving for roots again
                    
                for (K_counter = 0; K_counter < NumberOfKs; K_counter++)    {
                    retrn = Gllmj(workspace, l0, lf, m, j, (1.*K_counter+0.5)*K_step, K_counter, N, G_list+K_counter);
                    if (retrn != passed)  return 10;        }
                
                
                    // now integrate inv-Fourier over the Wigner functions
                
                for (theta_counter = 0; theta_counter < theta_step_num; theta_counter++)        {
                    for (phi_counter = 0; phi_counter < phi_step_num; phi_counter++)        {
                        theta = (1.*theta_counter + 0.5)*theta_step;
                        phi = (1.*phi_counter + 0.5)*phi_step;
                        
                        fillVector(&K, 1, theta, phi);
                        
                        GetAngles(&K, &u0, &uf, &n0, &nf, &delta_phi_f0, &delta_psi_f0);    // delta_*_f0 is in coordinate system where K // z
                        W_theta_0 = acos(dot(&u0, &K));
                        W_theta_f = acos(dot(&uf, &K));
                        
                        if (m == j)  loop_mj_iterations = 1;
                        else  loop_mj_iterations = 2;
                        
                        eval_m = m;
                        eval_j = j;
                        for (lmj_it = 1; lmj_it <= loop_mj_iterations; lmj_it++)    {
                            if (eval_m == 0)  lm_step = 1;      else  lm_step = 2*eval_m;
                            if (eval_j == 0)  lj_step = 1;      else  lj_step = 2*eval_j;
                            for (loop_m = -eval_m; loop_m <= eval_m;  loop_m += lm_step)        {
                            for (loop_j = -eval_j; loop_j <= eval_j;  loop_j += lj_step)        {
                                real_factor = exp(-N*(l_cutoff*l_cutoff + eta*loop_j*loop_j)) * Wigner(l0, loop_m, loop_j, W_theta_0)
                                                * Wigner(lf, loop_m, loop_j, W_theta_f) * sin(theta);
                                exp_factor.real = cos(loop_m*delta_phi_f0 + loop_j*delta_psi_f0 - loop_j*N*twist);      // no need to C.C. -- using loop_j rather than j
                                exp_factor.imag = sin(loop_m*delta_phi_f0 + loop_j*delta_psi_f0 - loop_j*N*twist);      // gives it the correct sign
                                
                                one_K_element = G_K_space + (theta_counter*phi_step_num + phi_counter)*NumberOfKs;
                                if ((loop_m*loop_j > 0) || ((loop_m*loop_j == 0) && (loop_m+loop_j >= 0)))   {      // use either computed Gllmj's or their complex conj's
                                    exp_factor = crMult(exp_factor, real_factor);
                                    for (K_counter = 0; K_counter < NumberOfKs; K_counter++)    {
                                        //*one_K_element = cAdd(*one_K_element, cMult(exp_factor, *(G_list+K_counter)));  <-- too slow
                                        one_K_element->real += exp_factor.real * (G_list+K_counter)->real - exp_factor.imag * (G_list+K_counter)->imag;
                                        one_K_element->imag += exp_factor.imag * (G_list+K_counter)->real + exp_factor.real * (G_list+K_counter)->imag;
                                        one_K_element++;                
                                }   }
                                else        {
                                    if ((l0+lf) % 2 == 0)  exp_factor = crMult(exp_factor, real_factor);
                                    else  exp_factor = crMult(exp_factor, -real_factor);           // undo the complex-conj. of the iK multiplier between the w_plus, Wlms
                                    for (K_counter = 0; K_counter < NumberOfKs; K_counter++)    {
                                            //*one_K_element = cAdd(*one_K_element, cMult(exp_factor, cConj(*(G_list+K_counter))));
                                        one_K_element->real += exp_factor.real * (G_list+K_counter)->real + exp_factor.imag * (G_list+K_counter)->imag;
                                        one_K_element->imag += exp_factor.imag * (G_list+K_counter)->real - exp_factor.real * (G_list+K_counter)->imag;
                                        one_K_element++;
                            }}  }   }
                            eval_temp = eval_m;
                            eval_m = eval_j;
                            eval_j = eval_temp;
                }   }   }
            }
            
            else if (sum_R == 2)        {               // if we sum over all space
                fillVector(&K, 1, 0, 0);
                
                GetAngles(&K, &u0, &uf, &n0, &nf, &delta_phi_f0, &delta_psi_f0);
                W_theta_0 = acos(dot(&u0, &K));
                W_theta_f = acos(dot(&uf, &K));
                
                if (m == j)  loop_mj_iterations = 1;
                else  loop_mj_iterations = 2;
                
                eval_m = m;
                eval_j = j;
                for (lmj_it = 1; lmj_it <= loop_mj_iterations; lmj_it++)    {
                    if (eval_m == 0)  lm_step = 1;      else  lm_step = 2*eval_m;
                    if (eval_j == 0)  lj_step = 1;      else  lj_step = 2*eval_j;
                    for (loop_m = -eval_m; loop_m <= eval_m;  loop_m += lm_step)        {
                    for (loop_j = -eval_j; loop_j <= eval_j;  loop_j += lj_step)        {
                        real_factor = exp(-N*(l0*(l0+1) + eta*loop_j*loop_j))
                                * Wigner(l0, loop_m, loop_j, W_theta_0) * Wigner(lf, loop_m, loop_j, W_theta_f);
                        exp_factor.real = cos(loop_m*delta_phi_f0 + loop_j*delta_psi_f0 - loop_j*N*twist);
                        exp_factor.imag = sin(loop_m*delta_phi_f0 + loop_j*delta_psi_f0 - loop_j*N*twist);
                        exp_factor = crMult(exp_factor, real_factor);
                        
                        *result = cAdd(*result, exp_factor);
                    }}
                    eval_temp = eval_m;
                    eval_m = eval_j;
                    eval_j = eval_temp;
            }   }
            
            else if (Cyclization == ccTrue)        {
                cdouble Gllm_temp;
                
                K_counter = 0;
                for (K_counter = 0; K_counter < NumberOfKs; K_counter++)    {
                    K_norm = (1.*K_counter+0.5)*K_step;
                    retrn = Gllmj(workspace, l0, lf, m, j, K_norm, K_counter, N, &Gllm_result);
                    if (retrn != passed)  return 10;
                    
                    Gllm_temp = Gllm_result;
                    exp_factor.real = cos(-j*N*twist);
                    exp_factor.imag = sin(-j*N*twist);
                    real_factor = exp(N*(-l_cutoff*l_cutoff - eta*j*j))*K_norm*K_norm;
                    
                    if (m > 0)  {  Gllm_temp.real *= 2;  Gllm_temp.imag = 0;  }
                    D_sum = cMult(crMult(Gllm_temp, real_factor), exp_factor);
                    if (j > 0)  {  D_sum.real *= 2;  D_sum.imag = 0;  }
                    
                    *result = cAdd(*result, D_sum);
                    
                    if (m != j)     {
                        exp_factor.real = cos(-m*N*twist);
                        exp_factor.imag = sin(-m*N*twist);
                        real_factor = exp(N*(-l_cutoff*l_cutoff - eta*m*m))*K_norm*K_norm;
                        
                        if (m > 0)  {  Gllm_result.real *= 2;  Gllm_result.imag = 0;  }
                        D_sum = cMult(crMult(Gllm_result, real_factor), exp_factor);
                        if (j > 0)  {  D_sum.real *= 2;  D_sum.imag = 0;  }
                        
                        *result = cAdd(*result, D_sum);                                 }
            }   }
        }}}}
    }}}}
    
    if (sum_R == 1)  {
        inv_Fourier_factor.imag = 0;
        R_norm = sqrt(dot(&R, &R));     }
    
    if ((Cyclization == ccFalse) && (sum_R != 2))      {
        for (theta_counter = 0; theta_counter < theta_step_num; theta_counter++)        {
            theta = (1.*theta_counter + 0.5)*theta_step;
            for (phi_counter = 0; phi_counter < phi_step_num; phi_counter++)      {
                phi = (1.*phi_counter + 0.5)*phi_step;
                
                fillVector(&K, 1, theta, phi);
                if (sum_R == 0)  unit_K_dot_R = dot(&K, &R);
                
                one_K_element = G_K_space + (theta_counter*phi_step_num + phi_counter)*NumberOfKs;
                for (K_counter = 0; K_counter < NumberOfKs; K_counter++)        {
                    K_norm = (1.*K_counter+0.5)*K_step;
                    if (sum_R == 0)     {
                        inv_Fourier_factor.real = cos(unit_K_dot_R*K_norm);
                        inv_Fourier_factor.imag = -sin(unit_K_dot_R*K_norm);        }
                    else                {
                        inv_Fourier_factor.real = sin(R_norm*K_norm)/K_norm;        }
                    *one_K_element = crMult(cMult(*one_K_element, inv_Fourier_factor), K_norm*K_norm);
                    *result = cAdd(*result, *one_K_element);
                    one_K_element++;
        }   }   }
        *result = crMult(*result, K_step*theta_step*phi_step);
        if (sum_R == 1)  *result = crMult(*result, pi_4*R_norm);
    }
    
    
        // Cyclization:  by rotating K the Ds should only integrate over tangents, whereas integral(G* x G) = 1 only when twists are integrated over as well.
        // For this reason so we need an extra 1/2pi.
    
    if (Cyclization == ccTrue)  *result = crMult(*result, one_over_16_pi_4*K_step);
    else  {
        if (sum_R != 2)  *result = crMult(*result, one_over_8_pi_3);
        if (sum_tangent == ccTrue)  *result = crMult(*result, pi_4);
        if (sum_twist == ccTrue)  *result = crMult(*result, pi_2);         }
    
    return passed;
}


// Returns the inverse-Laplace transformed G(l0, lf, m, j; p) given the roots in the numerator/denominator of the continued fraction

int Gllmj(linkedlist WorkLLs[2], int l0, int lf, int m, int j, double K, int K_index, double N, cdouble *result)
{
    cdouble NumDenom_Multiplier, *firstNumRoot = NULL, *firstDenomRoot = NULL;
    double AK_Multiplier;
    linkedlist *NumRoots, *DenomRoots;
    int NumberOfKs, l, linkedLretrn, poly_retrn;
    
    NumRoots = &(WorkLLs[0]);
    DenomRoots = &(WorkLLs[1]);
    
    NumberOfKs = (int) floor(0.5 + K_max/K_step);
    
    deleteElements(NumRoots, 1, NumRoots->elementNum);
    deleteElements(DenomRoots, 1, DenomRoots->elementNum);
    
    AK_Multiplier = 1;
    NumDenom_Multiplier = c_one;
    
    
        // first stack the num/denom root lists
    
    linkedLretrn = AddToRootList(NumRoots, DenomRoots, NumPolyCoef(Wlm_roots, l0, m, j), DenomPolyCoef(Wlm_roots, l0, m, j), &NumDenom_Multiplier, NumberOfKs, K_index);
    if (linkedLretrn != passed)  return 1;
        
    if (l0 < lf)        {
    for (l = l0+1; l <= lf; l++)    {
        linkedLretrn = AddToRootList(NumRoots, DenomRoots, NumPolyCoef(w_plus_roots, l, m, j), DenomPolyCoef(w_plus_roots, l, m, j), &NumDenom_Multiplier, NumberOfKs, K_index);
        if (linkedLretrn != passed)  return 1;
        AK_Multiplier *= K*K*a_lmj_sq(l, m, j);
        NumDenom_Multiplier = cMult(NumDenom_Multiplier, c_i);
    }}
    else        {
    for (l = l0-1; l >= lf; l--)    {
        linkedLretrn = AddToRootList(NumRoots, DenomRoots, NumPolyCoef(w_minus_roots, l, m, j), DenomPolyCoef(w_minus_roots, l, m, j), &NumDenom_Multiplier, NumberOfKs, K_index);
        if (linkedLretrn != passed)  return 1;
        AK_Multiplier *= K*K*a_lmj_sq(l+1, m, j);
        NumDenom_Multiplier = cMult(NumDenom_Multiplier, c_i);
    }}
    
    linkedLretrn = defragmentLinkedList(NumRoots);
    if (linkedLretrn == passed)  linkedLretrn = defragmentLinkedList(DenomRoots);
    if (linkedLretrn != passed)  return 1;
    
    
        // now call the root-canceller & multiply by the alpha factors
    
    if (NumRoots->elementNum > 0)  firstNumRoot = (cdouble *) element(NumRoots, 1);
    if (DenomRoots->elementNum > 0)  firstDenomRoot = (cdouble *) element(DenomRoots, 1);
    poly_retrn = InverseLaplace(firstNumRoot, NumRoots->elementNum, firstDenomRoot, DenomRoots->elementNum, N, &NumDenom_Multiplier);
    if (poly_retrn != passed)  {            // above line -- don't use element() since Num/DenomRoots may have 0 elements
        printf("Warning: root solver gave error %i; disregarding roots\n", poly_retrn);
        *result = c_zero;       }
    else        {
        result->real = NumDenom_Multiplier.real*sqrt(AK_Multiplier);
        result->imag = NumDenom_Multiplier.imag*sqrt(AK_Multiplier);                }
    
    return passed;
}


// Adds the roots in OneNum/DenomW (byte-strings) to Num/DenomRoots (cdoubles)
// Updates Multiplier, which is a constant prefactor, with the ratio of mult_num / mult_denom

int AddToRootList(linkedlist *NumRoots, linkedlist *DenomRoots, linkedlist *OneNumW, linkedlist *OneDenomW, cdouble *Multiplier, int NumberOfKs, int K_index)
{
    int linkedLretrn, c1, OldNumTop = NumRoots->elementNum, OldDenomTop = DenomRoots->elementNum;
    int NumRootsNum = OneNumW->elementNum/sizeof(cdouble)/NumberOfKs-1, DenomRootsNum = OneDenomW->elementNum/sizeof(cdouble)/NumberOfKs-1;
    cdouble *loopNumRoot, *loopDenomRoot;
    
    linkedLretrn = addElements(NumRoots, NumRootsNum, ccFalse);
    if (linkedLretrn == passed)  linkedLretrn = addElements(DenomRoots, DenomRootsNum, ccFalse);
    if (linkedLretrn == passed)  linkedLretrn = defragmentLinkedList(OneNumW);
    if (linkedLretrn == passed)  linkedLretrn = defragmentLinkedList(OneDenomW);
    if (linkedLretrn != passed)  return linkedLretrn;
    
    loopNumRoot = (cdouble *) element(OneNumW, K_index*(NumRootsNum+1)*sizeof(cdouble)+1);
    loopDenomRoot = (cdouble *) element(OneDenomW, K_index*(DenomRootsNum+1)*sizeof(cdouble)+1);
    *Multiplier = cMult(*Multiplier, cDiv(*loopNumRoot, *loopDenomRoot));
    
    for (c1 = OldNumTop+1; c1 <= NumRoots->elementNum; c1++)        {
        loopNumRoot++;
        *(cdouble *) element(NumRoots, c1) = *(loopNumRoot);       }
    for (c1 = OldDenomTop+1; c1 <= DenomRoots->elementNum; c1++)        {
        loopDenomRoot++;
        *(cdouble *) element(DenomRoots, c1) = *(loopDenomRoot);       }
    
    return passed;
}


// Multiplies 'result' by the inverse-Laplace transform of a rational function (specified by the num/denom root lists), using the residue method
// N is the conjugate of the Laplace variable, which is the normalized chain length

int InverseLaplace(cdouble *NumCoefs, int NumCoefTop, cdouble *DenomCoefs, int DenomCoefTop, double N, cdouble *result)
{
    cdouble complex_exp, *loopRoot, one_weighted_residue, OneFactorTerm, ResidueSum;
    int loopRootCounter, loopRoot2Counter, NumSimilarRoots, DenomSimilarRoots;
    const double SimilarityThreshold = 1e-15;               // max difference in two roots for them to be considered the same root
    
    ResidueSum.real = 0;
    ResidueSum.imag = 0;
    for (loopRootCounter = 0; loopRootCounter < DenomCoefTop; loopRootCounter++)    {
        loopRoot = DenomCoefs+loopRootCounter;
        complex_exp.real = exp((loopRoot->real)*N);
        complex_exp.imag = complex_exp.real;
        complex_exp.real *= cos((loopRoot->imag)*N);
        complex_exp.imag *= sin((loopRoot->imag)*N);
        
        NumSimilarRoots = 0;
        DenomSimilarRoots = 0;
        one_weighted_residue.real = 1;
        one_weighted_residue.imag = 0;
        
        
            // factor in the num/denom roots
            
        for (loopRoot2Counter = 0; loopRoot2Counter < NumCoefTop; loopRoot2Counter++)       {
            OneFactorTerm = cSub(*loopRoot, *(NumCoefs+loopRoot2Counter));
            if (OneFactorTerm.real*OneFactorTerm.real + OneFactorTerm.imag*OneFactorTerm.imag > SimilarityThreshold)  {
                one_weighted_residue = cMult(one_weighted_residue, OneFactorTerm);          }
            else  NumSimilarRoots += 1;                     }
        
        for (loopRoot2Counter = 0; loopRoot2Counter < DenomCoefTop; loopRoot2Counter++)     {
        if (loopRoot2Counter != loopRootCounter)        {
            OneFactorTerm = cSub(*loopRoot, *(DenomCoefs+loopRoot2Counter));
            if (OneFactorTerm.real*OneFactorTerm.real + OneFactorTerm.imag*OneFactorTerm.imag > SimilarityThreshold)  {
                one_weighted_residue = cDiv(one_weighted_residue, OneFactorTerm);           }
            else  DenomSimilarRoots += 1;                       }}
        
        
            // update the residue if it is a simple pole
        
        if (NumSimilarRoots == DenomSimilarRoots)       {
            one_weighted_residue = crMult(one_weighted_residue, 1./(NumSimilarRoots+1));
            one_weighted_residue = cMult(one_weighted_residue, complex_exp);
            ResidueSum = cAdd(ResidueSum, one_weighted_residue);            }
    }
    
    *result = cMult(ResidueSum, *result);
    
    return passed;
}


// Calculates the change in the 'phi' and 'psi' Euler angles between the two ends of the chain

void GetAngles(vector_3 *K, vector_3 *u0, vector_3 *uf, vector_3 *n0, vector_3 *nf, double *delta_phi_f0, double *delta_psi_f0)
{
    vector_3 u0_perp;
    vector_3 uf_perp;
    double psi_0, psi_f;
    
    GetOneAngle(u0, n0, K, &u0_perp, &psi_0);
    GetOneAngle(uf, nf, K, &uf_perp, &psi_f);
    *delta_phi_f0 = AngleBetween(u0_perp, uf_perp, *K);
    *delta_psi_f0 = psi_f - psi_0;    
    
/*    {             // Don't delete (used for checking)
        vector_3 temp_x;
        
        fillVector(&temp_x, 1, pi/2, 0);
        printf("%g, %g, %g    ----    ", acos(dot(u0, K)), AngleBetween(temp_x, u0_perp, *K), psi_0);
        printf("%g, %g, %g\n", acos(dot(uf, K)), AngleBetween(temp_x, uf_perp, *K), psi_f);
    } */
}


// The Spakowitz method defines 'z' as the direction of the K vector, so the 'n' lies in the K-u plane if the third Euler angle 'psi' is zero.
// GetOneAngle() calculates psi as between the in-plane 'n0' and the given n, rotated in the sense of u.
// Also returns the component of u perp. to K.

void GetOneAngle(vector_3 *u, vector_3 *n, vector_3 *K, vector_3 *u_perp, double *psi)
{
    vector_3 n_rot, neg_u_par;
    double u_dot_K = dot(u, K);
    
    if (1. - u_dot_K < 1.e-15)  {
        *u_perp = *n;
        *psi = 0;                }             // magnitude doesn't matter for GetAngles()
    else if (1.+u_dot_K < 1.e-15)  {
        v3ScalarMult(n, -1., u_perp);
        *psi = 0;                    }
    else        {
        v3ScalarMult(K, -u_dot_K, &neg_u_par);
        v3Add(u, &neg_u_par, u_perp);
        
        RotateVector(K, u, u_perp, &n_rot);
        *psi = AngleBetween(n_rot, *n, *u);       }
}


// Returns the angle needed to rotate v0 into vf.
// The out-of-plane component of v_clockwise determines the positive sense of rotation.

double AngleBetween(vector_3 v0, vector_3 vf, vector_3 v_clockwise)
{
    double answer, norm, norm_dot;
    vector_3 cross_0f;
    
    norm = sqrt(dot(&v0, &v0)*dot(&vf, &vf));
    if (norm == 0)  return 0;
    
    norm_dot = dot(&v0, &vf)/norm;
    if (norm_dot < -1)  norm_dot = -1;          // catch numerical problems close to +-1
    else if (norm_dot > 1)  norm_dot = 1;
    answer = acos(norm_dot);        // [0, pi] -- assumes v1 X v2 > 0
    
    cross_0f = cross(&v0, &vf);
    if (dot(&cross_0f, &v_clockwise) > 0)  return answer;
    else  return -1.0*answer;
}


// Writes a vector having the given polar coordinates

inline void fillVector(vector_3 *V, double norm, double theta, double phi)
{
    V->x = norm * sin(theta) * cos(phi);
    V->y = norm * sin(theta) * sin(phi);
    V->z = norm * cos(theta);
}


// Applies the same rotation n --> n_rot as would rotate K onto u
// Assumes both u and K are unit vectors, and that they are not parallel

void RotateVector(vector_3 *K, vector_3 *u, vector_3 *n, vector_3 *n_rot)
{
    vector_3 omega, n_par, n_perp_old, n_perp, n_cross, n_perp_rot;
    double omega_sq;
    
    omega = cross(K, u);
    omega_sq = dot(&omega, &omega);
    if (omega_sq < 1.e-30)      {
        *n_rot = *n;
        return;             }
    else if (omega_sq > 1)  omega_sq = 1;
    
    v3ScalarMult(&omega, dot(n, &omega)/omega_sq, &n_par);
    v3Sub(n, &n_par, &n_perp_old);
    if (dot(K, u) < 0)
        v3ScalarMult(&n_perp_old, -sqrt(1-omega_sq), &n_perp);
    else
        v3ScalarMult(&n_perp_old, sqrt(1-omega_sq), &n_perp);
    n_cross = cross(&omega, &n_perp_old);
    
    v3Add(&n_perp, &n_cross, &n_perp_rot);
    v3Add(&n_perp_rot, &n_par, n_rot);
}


// Normalizes vector 'v' to 1

int Renorm_3v(vector_3 *v)
{
    double norm = sqrt(dot(v, v));
    
    if (norm == 0)  return 1;
    v->x /= norm;
    v->y /= norm;
    v->z /= norm;
    
    return passed;
}


// Orthogonalizes v1 with respect to v2
// (assumes v1 is a unit vector)

void GramSchmidt(vector_3 *v1, vector_3 *v2)
{
    vector_3 parallel_comp;
    
    vScalarMult((double *) v2, 3, dot(v1, v2), (double *) &parallel_comp);
    vSub((double *) v1, (double *) &parallel_comp, 3, (double *) v1);
}



// Next 2 routines:  return the (Yamakawa-defined) Wigner function, excluding the complex exponentials

int call_Wigner(int argc, char **argv)
{
    linkedlist workspace;
    int retrn, l, m, j;
    double theta, *ans;
    arg_info *ArgInfo = (arg_info *) argv[argc];
    
    const int ArgTypes[] = { int_type, int_type, int_type, double_type, double_type };
    const int ArgIndices[] = { 1, 1, 1, 1, 1 };
    
    if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "Wigner") != passed)  return 1;
    
    getArgs(argc, argv, byValue(&l), byValue(&m), byValue(&j), byValue(&theta), &ans);
    
    if ((m*m > l*l) || (j*j > l*l)) {
        printf("Wigner():  |l| must be >= |m|, |j|\n");
        return 2;           }
    
    FactorialList = &workspace;
    FactorialList->memory = 0;
    retrn = newLinkedList(FactorialList, 1, sizeof(double), 255, ccFalse);        // first element is cleared to zero
    if (retrn != passed)  return retrn;
    *LL_Double(FactorialList, 1) = 1;
    
    *ans = Wigner(l, m, j, theta);
    
    deleteLinkedList(FactorialList);
    
    return passed;
}

double Wigner(int l, int m, int j, double theta)
{
    int i, n, hold_j, alpha, beta;
    double x, n_plus_i_pow, answer, FinalMultiplier = 1;
    const double one_over_8pipi = 1./(8*pi*pi);
    
    if (abs(m) > abs(j))    {
        if ((l+m) % 2 != 0)  FinalMultiplier *= -1;
        hold_j = j;
        j = -m;
        m = hold_j;
        theta = pi-theta;       }
    
    if (j < 0)      {
        if ((l+j) % 2 != 0)  FinalMultiplier *= -1;
        j = -j;
        theta = pi+theta;       }
    
    n = l-j;
    alpha = j-m;
    beta = j+m;
    x = cos(theta);
    
    answer = 0;
    n_plus_i_pow = 1;
    
    for (i = 0; i <= n; i++)        {
        answer += (Factorial(alpha+n) / Factorial(alpha+n-i) / Factorial(i)) * (Factorial(beta+n) / Factorial(beta+i) / Factorial(n-i)) * IntegerPow(1-x, n-i) * n_plus_i_pow;
        n_plus_i_pow *= -(1+x);         }           // -1 since each term is a higher derivative of 1-x
    
        // NO divide by n! -- the reason is that, in taking derivatives, the different orderings contribute, so answer should be multiplied by n! which cancels the 1/n!
    return answer * FinalMultiplier * IntegerPow(cos(theta/2), j+m) * IntegerPow(sin(theta/2), j-m) / IntegerPow(-2, n)
            * sqrt( (2*l + 1) * one_over_8pipi * Factorial(l+j) * Factorial(l-j) / Factorial(l+m) / Factorial(l-m) );
}



// Copies one polynomial into another

int CopyPoly(linkedlist *SourcePoly, linkedlist *DestPoly)
{
    int linkedLretrn;
    
    linkedLretrn = ResizePoly(DestPoly, PolyVariables(SourcePoly), PolyTerms(SourcePoly));
    if (linkedLretrn != passed)  return linkedLretrn;
    
    copyElements(SourcePoly, 1, DestPoly, 1, SourcePoly->elementNum);
    return passed;
}



// Next 3 routines:  add, multiply and scalar-multiply polynomials
// Note that the scalar-multiply actually multiplies by a 'term':  scalar * var1 * var2 ...

int AddPoly(linkedlist *DestPoly, linkedlist *PolyToAdd)
{
    int linkedLretrn, FirstCopiedTerm_ID;
    
    FirstCopiedTerm_ID = PolyTerms(DestPoly)+1;
    linkedLretrn = ResizePoly(DestPoly, PolyVariables(DestPoly), PolyTerms(PolyToAdd)+PolyTerms(DestPoly));
    if (linkedLretrn != passed)  return linkedLretrn;
    copyElements(PolyToAdd, PolyElementIndex(PolyToAdd, 1), DestPoly, PolyElementIndex(DestPoly, FirstCopiedTerm_ID),
//              PolyTerms(PolyToAdd)*sizeof(struct{ double Coef;  int Pows[PolyVariables(DestPoly)]; }));
                PolyTerms(PolyToAdd)*(sizeof(double) + PolyVariables(DestPoly)*sizeof(int)));
    
    SimplifyPoly(DestPoly);
    return passed;
}

int MultiplyPoly(linkedlist *Poly1, linkedlist *Poly2, linkedlist *DestPoly)
{
    int Term1Counter, Term2Counter, VariableCounter, DestTop, DestCounter;
    int MaxVariable = PolyVariables(Poly1), MaxTerm1 = PolyTerms(Poly1), MaxTerm2 = PolyTerms(Poly2);
    double *Term1, *Term2, *DestTerm;
    int linkedLretrn;
    
    ResizePoly(DestPoly, PolyVariables(Poly1), 0);      // get it to clear the memory
    linkedLretrn = ResizePoly(DestPoly, PolyVariables(Poly1), PolyTerms(Poly1)*PolyTerms(Poly2));
    if (linkedLretrn != passed)  return linkedLretrn;
    DestTop = 0;
    DestTerm = 0;           // to get rid of compile-time warning
    
    for (Term1Counter = 1; Term1Counter <= MaxTerm1; Term1Counter++)        {
    for (Term2Counter = 1; Term2Counter <= MaxTerm2; Term2Counter++)        {
        Term1 = PolyElement(Poly1, Term1Counter);
        Term2 = PolyElement(Poly2, Term2Counter);
        for (DestCounter = 1; DestCounter <= DestTop; DestCounter++)        {
            DestTerm = PolyElement(DestPoly, DestCounter);
            for (VariableCounter = 1; VariableCounter <= MaxVariable; VariableCounter++)        {
            if (*PolyCoef(DestTerm, VariableCounter) != *PolyCoef(Term1, VariableCounter) + *PolyCoef(Term2, VariableCounter))      {
                break;
            }}
            if (VariableCounter == MaxVariable+1)  break;
        }
        
        if (DestCounter == DestTop+1)       {
            DestTerm = PolyElement(DestPoly, DestCounter);
            DestTop = DestCounter;
            for (VariableCounter = 1; VariableCounter <= MaxVariable; VariableCounter++)        {
                *PolyCoef(DestTerm, VariableCounter) = *PolyCoef(Term1, VariableCounter) + *PolyCoef(Term2, VariableCounter);       }
        }
        *DestTerm += (*Term1) * (*Term2);
    }}
    
    ResizePoly(DestPoly, PolyVariables(Poly1), DestTop);        // get it to clear the memory
    return passed;
}

void ScalarMultiplyPoly(linkedlist *thePoly, double *ConstantMultiplier)
{
    int c1, c2, TopVar = PolyVariables(thePoly), TopTerm = PolyTerms(thePoly);
    double *loopTerm;
    
    for (c1 = 1; c1 <= TopTerm; c1++)       {
        loopTerm = PolyElement(thePoly, c1);        // multiply the coefficients
        *loopTerm *= *ConstantMultiplier;
        for (c2 = 1; c2 <= TopVar; c2++)        {
            *PolyCoef(loopTerm, c2) += *PolyCoef(ConstantMultiplier, c2);   // multiply the powers of each variable
    }   }
}



// Evaluates one of the variables in a polynomial (i.e. p or K), but leaves the power of that variable in place.
// This is because we reduce the iK variable by evaluating the Ks and reinterpret the power as that of the 'i's alone.

int ReducePoly_iK(linkedlist *FullPoly, linkedlist *ReducedPoly, int CollapsedVariable, double CV_Value)
{
    int linkedLretrn, c1, c2, PolyTop = PolyTerms(FullPoly), OldPolyVars = PolyVariables(FullPoly);
    double *loopFullTerm, *loopReducedTerm;
    
    linkedLretrn = ResizePoly(ReducedPoly, OldPolyVars, PolyTop);
    if (linkedLretrn != passed)  return linkedLretrn;
    
    for (c1 = 1; c1 <= PolyTop; c1++)       {
        loopFullTerm = PolyElement(FullPoly, c1);
        loopReducedTerm = PolyElement(ReducedPoly, c1);
        
        *loopReducedTerm = *loopFullTerm;
        for (c2 = 1; c2 <= *PolyCoef(loopFullTerm, CollapsedVariable); c2++)        {
            *loopReducedTerm *= CV_Value;           }
        
        for (c2 = 1; c2 <= OldPolyVars; c2++)       {
            *PolyCoef(loopReducedTerm, c2) = *PolyCoef(loopFullTerm, c2);       }
    }
    
    iReducePoly(ReducedPoly, CollapsedVariable);
    
    return passed;
}



// Sets each power of i in a polynomial to 0 or 1, negating the coefficient if necessary; then simplifies

void iReducePoly(linkedlist *thePoly, int i_term)
{
    int TermCounter, PolyTop = PolyTerms(thePoly);
    double *loopTerm;
    
    for (TermCounter = 1; TermCounter <= PolyTop; TermCounter++)        {
        loopTerm = PolyElement(thePoly, TermCounter);
        if ((*PolyCoef(loopTerm, i_term) % 4) >= 2)     {
            *loopTerm *= -1;            }
        *PolyCoef(loopTerm, i_term) = *PolyCoef(loopTerm, i_term) % 2;
    }
    
    SimplifyPoly(thePoly);
}



// Complex-conjugates a polynomial

void ComplexConjugatePoly(linkedlist *thePoly, int i_term)
{
    int TermCounter, PolyTop = PolyTerms(thePoly);
    double *loopTerm;
    
    for (TermCounter = 1; TermCounter <= PolyTop; TermCounter++)        {
        loopTerm = PolyElement(thePoly, TermCounter);
        if (*PolyCoef(loopTerm, i_term) % 2 == 1)  *loopTerm *= -1;     }
}



// Evaluates the Nth derivative (a vector in each variable besides i) / prod(N!)
// of a polynomial for given values of the variables (also a vector).

cdouble EvaluatePoly(linkedlist *thePoly, cdouble *CV_Values, int *Derivatives, int i_term)
{
    cdouble Answer, SubAnswer;
    double *OneTerm, temp;
    int c1, c2, c3, loopCoef, loopDeriv, VarTop = PolyVariables(thePoly), TermTop = PolyTerms(thePoly);
    
    Answer.real = 0;
    Answer.imag = 0;
    
    for (c1 = 1; c1 <= TermTop; c1++)       {
        OneTerm = PolyElement(thePoly, c1);
        SubAnswer.real = *OneTerm;
        SubAnswer.imag = 0;
        for (c2 = 1; c2 <= VarTop; c2++)        {
            loopCoef = *PolyCoef(OneTerm, c2);
            if (c2 == i_term)       {
            if (loopCoef == 1)      {
                temp = SubAnswer.imag;
                SubAnswer.imag = SubAnswer.real;
                SubAnswer.real = -temp;
            }}
            else        {
                loopDeriv = *(Derivatives+c2-1);
                if (loopDeriv <= loopCoef)  {
                    for (c3 = 1; c3 <= loopCoef-loopDeriv; c3++)        {
                        SubAnswer = cMult(SubAnswer, *(CV_Values+c2-1));            }
                    for (c3 = 0; c3 < loopDeriv; c3++)      {
                        SubAnswer.real *= (loopCoef-c3);
                        SubAnswer.imag *= (loopCoef-c3);            }
                    for (c3 = 2; c3 <= loopDeriv; c3++)     {
                        SubAnswer.real /= c3;
                        SubAnswer.imag /= c3;           }       }
                else        {
                    SubAnswer.real = 0;
                    SubAnswer.imag = 0;     }
        }   }
        Answer.real += SubAnswer.real;
        Answer.imag += SubAnswer.imag;
    }
    
    return Answer;
}



// Sums terms with the same powers of the variables into a single term, and eliminates terms with coefficients of zero

void SimplifyPoly(linkedlist *thePoly)
{
    int Term1Counter, Term2Counter, VariableCounter, MaxVariable = PolyVariables(thePoly), MaxTerm = PolyTerms(thePoly);
    double *Term1, *Term2;
    
    for (Term1Counter = 1; Term1Counter <= MaxTerm; Term1Counter++) {
    for (Term2Counter = Term1Counter+1; Term2Counter <= MaxTerm; Term2Counter++)    {
        
        Term1 = PolyElement(thePoly, Term1Counter);
        Term2 = PolyElement(thePoly, Term2Counter);
        
        for (VariableCounter = 1; VariableCounter <= MaxVariable; VariableCounter++)    {
        if (*PolyCoef(Term1, VariableCounter) != *PolyCoef(Term2, VariableCounter)) {
            break;
        }}
        
        if (VariableCounter == MaxVariable+1)   {
            *Term1 += *Term2;
            *Term2 = 0;             }
    }}
    
    for (Term1Counter = 1; Term1Counter <= MaxTerm; Term1Counter++)     {
    if (*PolyElement(thePoly, Term1Counter) == 0)       {
        ((PolyType *) element(thePoly, 1))->terms_num--;
        deleteElements(thePoly, PolyElementIndex(thePoly, Term1Counter), PolyElementIndex(thePoly, Term1Counter+1)-1);
        Term1Counter--;
        MaxTerm--;
    }}
}



// Sizes the polynomial to store a given number of variables and terms

int ResizePoly(linkedlist *thePoly, int VariablesNum, int TermsNum)
{
    PolyType *PolyHeader;
    int RequiredStorage, linkedLretrn;
    
    RequiredStorage = 2*sizeof(int) + TermsNum*sizeof(double) + TermsNum*VariablesNum*sizeof(int);
    if (RequiredStorage > thePoly->elementNum)      {
        linkedLretrn = addElements(thePoly, (RequiredStorage-thePoly->elementNum), ccTrue);
        if (linkedLretrn != passed)  return linkedLretrn;           }
    else if (RequiredStorage < thePoly->elementNum)     {
        deleteElements(thePoly, RequiredStorage+1, thePoly->elementNum);       }
    
    linkedLretrn = defragmentLinkedList(thePoly);
    if (linkedLretrn != passed)  return linkedLretrn;
    
    PolyHeader = (PolyType *) element(thePoly, 1);
    PolyHeader->variables_num = VariablesNum;
    PolyHeader->terms_num = TermsNum;
    
    return passed;
}



// Sizes thePoly with the correct number of variables & terms, and initializes it with 0s or, optionally, InitialValues

int InitPoly(linkedlist *thePoly, int VariablesNum, int TermsNum, char *InitialValues)
{
    int linkedLretrn;
    
    if (InitialValues == 0)     {               // to force it to clear the new memory
        linkedLretrn = ResizePoly(thePoly, VariablesNum, 0);
        if (linkedLretrn != passed)  return linkedLretrn;           }
    
    linkedLretrn = ResizePoly(thePoly, VariablesNum, TermsNum);
    if (linkedLretrn != passed)  return linkedLretrn;
    
    if (InitialValues != 0)     {
        setElements(thePoly, PolyElementIndex(thePoly, 1), thePoly->elementNum, (void *) InitialValues);
    }
    
    return passed;
}



// Next 2 routines:  return the given element/its pointer in a polynomial given its byte-array

inline int PolyElementIndex(linkedlist *thePoly, int theElement)
{
    int PolyVars = PolyVariables(thePoly);
    if (PolyVars % 2 == 1)  PolyVars += 1;
    return 2*sizeof(int) + (theElement-1)*sizeof(double) + (theElement-1)*PolyVars*sizeof(int) + 1;
}

double *PolyElement(linkedlist *thePoly, int theElement)
{
    return (double *) element(thePoly, PolyElementIndex(thePoly, theElement));
}



// Next 2 routines:  return the number of variables/terms of the given polynomial

inline int PolyVariables(linkedlist *thePoly)
{
    return ((PolyType *) (((char *) thePoly->memory)+sublistHeaderSize))->variables_num;
}

int PolyTerms(linkedlist *thePoly)
{
    return ((PolyType *) (((char *) thePoly->memory)+sublistHeaderSize))->terms_num;
}



// Returns the pointer to the power of the given variable in a given term of a polynomial

int *PolyCoef(double *theTerm, int VariableNumber)
{
    return ((int *) (theTerm+1)) + VariableNumber - 1;
}



// Next 2 routines:  return the polynomial in the numerator/denominator of the { l, i, j, m } fraction in one of the w arrays

inline linkedlist *NumPolyCoef(linkedlist *theTables, int l, int m, int j)
{
    linkedlist *theList = theTables;
    int loop_l;
    
    for (loop_l = 0; loop_l < l; loop_l++)
        theList += 2*TriangleSum(1, loop_l+1);
    
    return theList + 2*(TriangleSum(1, m) + j);
}

inline linkedlist *DenomPolyCoef(linkedlist *theTables, int l, int m, int j)
{
    linkedlist *theList = theTables;
    int loop_l;
    
    for (loop_l = 0; loop_l < l; loop_l++)
        theList += 2*TriangleSum(1, loop_l+1);
    
    return theList + 2*(TriangleSum(1, m) + j) + 1;
}



// Alternative to pow() for integers -- faster?

double IntegerPow(double x, int pow)
{
    int c1;
    double answer = 1;
    
    for (c1 = 1; c1 <= pow; c1++)
        answer *= x;
    
    return answer;
}



// Next 2 routines:  return [log] x!

double LogFactorial(int x)
{
    int c1, OldTop;
    
    if (x+1 > LogFactorialList->elementNum)     {
        OldTop = LogFactorialList->elementNum;
        addElements(LogFactorialList, x + 1 - OldTop, ccFalse);      // ignore return value
        for (c1 = OldTop+1; c1 <= x+1; c1++)    {
            *LL_Double(LogFactorialList, c1) = *LL_Double(LogFactorialList, c1-1) + log(c1-1);      }
    }
    return *LL_Double(LogFactorialList, x+1);
}

double Factorial(int x)
{
    int c1, OldTop;
    
    if (x+1 > FactorialList->elementNum)        {
        OldTop = FactorialList->elementNum;
        addElements(FactorialList, x + 1 - OldTop, ccFalse);     // ignore return value
        defragmentLinkedList(FactorialList);
        for (c1 = OldTop+1; c1 <= x+1; c1++)    {
            *LL_Double(FactorialList, c1) = *LL_Double(FactorialList, c1-1)*(c1-1);     }
    }
    
    return *LL_Double(FactorialList, x+1);
}
