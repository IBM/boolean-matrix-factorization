// Copyright (C) 2019, International Business Machines
// Corporation.  All Rights Reserved.

// This program is distributed under the terms of the
// Eclipse Public License - v 2.0

#include "Vec.hpp"
#include "Alg1.hpp"
#include "BuildMatA.hpp"
#include "Alg7.hpp"
#include "Alg8.hpp"
#include "Alg6.hpp"

#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#ifndef WIN32
#include <sys/times.h>
#endif

using namespace std;


void readInputData(string fname, int& n, int& m, Vec<Vec<int> >& x);
int evaluateSolution(Vec<Vec<int> > &x, Vec<Vec<int> > &B, Vec<Vec<int> > &Stranp);
Vec<Vec<int> > product(Vec<Vec<int> > &Stranp, Vec<Vec<int> > &B, int rowsS, int colsB);
int optimize(Vec<Vec<int> > &X, int &k, double &trsh, int Alg, Vec<Vec<int> >& B, Vec<Vec<int> >& Stranp);
void perturbB(Vec<Vec<int> > &B);
void doubleLoop(Vec<Vec<int> > &X, int &k, double &tresholdStep,
                double &perturbStep, double &maxPerturb,
                int Alg);
int tresholdLoop(Vec<Vec<int> > &X, int &k, double &tresholdStep,
                 double &optTresh, int &Alg, Vec<Vec<int> >& currB, Vec<Vec<int> >& currStranp);
void copySol(Vec<Vec<int> >& inB, Vec<Vec<int> >& inStranp, Vec<Vec<int> >& outB, Vec<Vec<int> >& outStranp);
void writeBandS(Vec<Vec<int> >& B, Vec<Vec<int> >& Stranp);

double toPerturb=-1.;

int main(int argc, const char * argv[]) {
    // insert code here...
    // start time measurement
    double t0;
    struct tms timearr; //clock_t tres;
    //tres = times(&timearr);
    times(&timearr);
    t0 = timearr.tms_utime;
    
    // read input data file name
    assert(argc==7);
    string inputFileName=argv[1];
    int k=atoi(argv[2]);
    double tresholdStep=atof(argv[3]);
    double perturbStep=atof(argv[4]);
    double maxPerturb=atof(argv[5]);
    int Alg=atoi(argv[6]);
    cout<<"inputFileName="<<inputFileName<<", k="<<k<<", tresholStep="<<tresholdStep<<
    ", perturbStep="<<perturbStep<<
    ", maxPerturbation="<<maxPerturb<<
    ", Alg "<< Alg<<endl;
    

    
    int n, m;
    Vec<Vec<int> > X;
    readInputData(inputFileName,n, m, X);
    printf("n= %d, m= %d\n", n, m);
    
    doubleLoop(X, k, tresholdStep, perturbStep, maxPerturb, Alg);
    
    
    
    // end time measurement
    //tres = times(&timearr);
    times(&timearr);
    double t = (timearr.tms_utime-t0)/100.;
    printf("Total Time: %f secs\n", t);
    
    return 0;
}
void doubleLoop(Vec<Vec<int> > &X, int &k, double &tresholdStep,
                double &perturbStep, double &maxPerturb,
                int Alg){
    Vec<Vec<int> > bestB;
    Vec<Vec<int> > bestStranp;
    Vec<Vec<int> > currB;
    Vec<Vec<int> > currStranp;
    int nrows = (int) X.size();
    int ncols = (int) X[0].size();
    cout<<"DEBUG: nrows = "<<nrows<<", ncols = "<<ncols<<endl;
    for(int i=0; i<k; i++) {
      Vec<int> row1(ncols);
      bestB.push_back(row1);
      Vec<int> row2(ncols);
      currB.push_back(row2);
      Vec<int> col1(nrows);
      bestStranp.push_back(col1);
      Vec<int> col2(nrows);
      currStranp.push_back(col2);
    }
    int d5 = Alg %10;
    double optTresh=-1.;
    if ( d5==1 ){ // loop for perturbation
        toPerturb=0.0;
        int minValue=9999999;
        double minPert=-1.;
        double minTresh=-1.;
        while (toPerturb < maxPerturb){
            printf("\nPerturbation %g --------------------------\n",toPerturb);
            int value=tresholdLoop(X, k, tresholdStep, optTresh, Alg, currB, currStranp);
            if ( value < minValue){
                minValue=value;
                minPert=toPerturb;
                minTresh=optTresh;
		copySol(currB, currStranp, bestB, bestStranp);
            }
            toPerturb+=perturbStep;
        }
        printf("\nBest value %d treshold %g perturbation %g\n",
               minValue, minTresh, minPert);
    }
    else { // no perturbation
        int value=tresholdLoop(X, k, tresholdStep, optTresh, Alg, currB, currStranp);
	copySol(currB, currStranp, bestB, bestStranp);
        printf("\nBest value %d treshold %g \n", value, optTresh);
    }
    writeBandS(bestB,bestStranp);
}
int tresholdLoop(Vec<Vec<int> > &X, int &k, double &tresholdStep,
                 double &optTresh, int &Alg, Vec<Vec<int> >& currB, Vec<Vec<int> >& currStranp){
    int minValue=9999999;
    double minTrsh=-1.;
    double trsh=tresholdStep;
    while ( trsh < 1. ){
        printf("\nTreshold %g\n",trsh);
	Vec<Vec<int> > B;
	Vec<Vec<int> > Stranp;
        int value=optimize(X, k, trsh, Alg, B, Stranp);
        if ( value < minValue){
            minValue=value;
            minTrsh=trsh;
            copySol(B, Stranp, currB, currStranp);
        }
        trsh+=tresholdStep;
    }
    //printf("\nBest value %d treshold %g\n", minValue, minTrsh);
    optTresh=minTrsh;
    return minValue;
}

int optimize(Vec<Vec<int> > &X, int &k, double &trsh, int Alg, Vec<Vec<int> >& B, Vec<Vec<int> >& Stranp)
{
    // get digits from Alg
    int d5 = Alg %10;
    Alg /= 10;
    int d4 = Alg %10;
    Alg /= 10;
    int d3 = Alg %10;
    Alg /= 10;
    int d2 = Alg %10;
    Alg /= 10;
    int d1 = Alg %10;
    printf("Algorithms %d %d %d %d %d\n", d1, d2, d3, d4, d5);

    BuildMatA buildA(X, k, trsh);
    Vec<Vec<int> > &A=buildA.A;
    //    Vec<Vec<int> > B;
    //    Vec<Vec<int> > Stranp;
    if (d1==1){
        Alg1 alg1(X, k, A);
        B=alg1.B;
        Stranp=alg1.Stranp;
    }
    else
    {
        Alg6 alg6(X, k, A);
        B=alg6.B;
        Stranp=alg6.Stranp;
    }
    
    int result=evaluateSolution(X, B, Stranp);
    if (d2==7){
        Alg7 alg7(X, k, B);
        Stranp=alg7.Stranp;
        printf("Alg7 ");
        result=evaluateSolution(X, B, Stranp);
        if (d5==1){
            perturbB(B);
            printf("perturbed ");
            result=evaluateSolution(X, B, Stranp);
        }
    }
    if (d2==8){
        Alg8 alg8(X, k, Stranp, B);
        printf("Alg8 ");
        result=evaluateSolution(X, B, Stranp);
        if (d5==1){
            perturbB(B);
            printf("perturbed ");
            result=evaluateSolution(X, B, Stranp);
        }
    }
    if (d3==8){
        Alg8 alg8(X, k, Stranp, B);
        printf("Alg8 ");
        result=evaluateSolution(X, B, Stranp);
    }
    if (d4==7){
        Alg7 alg7(X, k, B);
        Stranp=alg7.Stranp;
        printf("Alg7 ");
        result=evaluateSolution(X, B, Stranp);
    }
    return result;
}

int evaluateSolution(Vec<Vec<int> > &x, Vec<Vec<int> > &B, Vec<Vec<int> > &Stranp)
{
    const int nrows=(int)x.size();
    const int ncols=(int)x[0].size();
    Vec<Vec<int> > SB=product(Stranp, B, nrows, ncols);
    int count=0;
    for (int i=0; i<nrows; ++i)
        for (int j=0; j<ncols; ++j)
            if ( x[i][j] != SB[i][j] ) ++count;
    printf("discrepancy %d \n", count);
    return count;
}
Vec<Vec<int> > product(Vec<Vec<int> > &Stranp, Vec<Vec<int> > &B, int rowsS, int colsB)
{
    Vec<Vec<int> > SB;
    const int rowsB=(int)B.size();
    if ( rowsB==0 ){
        Vec<int> row(colsB,0);
        for ( int irow=0; irow < rowsS; ++irow)
            SB.push_back(row);
    }
    else {
        Vec<int> row(colsB);
        for (int irow=0; irow<rowsS; ++irow){
            for (int icol=0; icol<colsB; ++icol){
                int sum=0;
                for (int k=0; k<rowsB; ++k)
                    sum+=Stranp[k][irow]*B[k][icol];
                row[icol]= (sum > 0) ? 1 : 0;
            }
            SB.push_back(row);
        }
    }
    return SB;
}
void perturbB(Vec<Vec<int> > &B)
{
    const int nrows=(int)B.size();
    const int ncols=(int)B[0].size();
    for (int i=0; i<nrows; ++i){
        for (int j=0; j < ncols; ++j) {
            if (drand48() < 1. - toPerturb)continue;
            B[i][j]= 1 - B[i][j];
        }
    }
}

void copySol(Vec<Vec<int> >& inB, Vec<Vec<int> >& inStranp, Vec<Vec<int> >& outB, Vec<Vec<int> >& outStranp)
{
  int k = (int)inB.size();
  int colsB = (int)inB[0].size();
  int rowsS = (int)inStranp[0].size();
  for(int i=0; i<k; i++) {
    for(int j=0; j<colsB; j++)
      outB[i][j] = inB[i][j];
    for(int j=0; j<rowsS; j++)
      outStranp[i][j] = inStranp[i][j];
  }
}

void writeBandS(Vec<Vec<int> >& B, Vec<Vec<int> >& Stranp)
{
    int k = (int)B.size();
    int colsB = (int)B[0].size();
    int rowsS = (int)Stranp[0].size();
    FILE *sfile=fopen("Smatrix.txt","w");
    for (int i=0; i<rowsS; ++i){
        for (int j=0; j<k; ++j) fprintf(sfile,"%d ",Stranp[j][i]);
        fprintf(sfile,"\n");
    }
    fclose(sfile);
    FILE *bfile=fopen("Bmatrix.txt","w");
    for (int i=0; i<k; ++i){
        for (int j=0; j<colsB; ++j) fprintf(bfile,"%d ",B[i][j]);
        fprintf(bfile,"\n");
    }
    fclose(bfile);

}
