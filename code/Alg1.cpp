// Copyright (C) 2019, International Business Machines
// Corporation.  All Rights Reserved.

// This program is distributed under the terms of the
// Eclipse Public License - v 2.0

#include "Vec.hpp"
#include "Alg1.hpp"

Alg1::Alg1(Vec<Vec<int> > &C, int k, Vec<Vec<int> > &A){
    wp=wm=1.0;
    for (int l=0; l < k; ++l){
        int i=-1;
        Vec<int> s;
        find_max(C, A, i, s);
        B.push_back(A[i]);
        Stranp.push_back(s);
    }
}
void Alg1::find_max(Vec<Vec<int> > &C, Vec<Vec<int> > &A, int &imax, Vec<int> &smax)
{
    const int m=(int)A.size();
    const int rowsS=(int)C.size();
    const int colsB=(int)C[0].size();
    double vmax=-1.e31;
    imax=-1;
    Vec<Vec<int> > SB = product(rowsS, colsB);
    for (int k=0; k<m; ++k){ // rows of A
        Vec<int> s(rowsS,-1);
        for (int is=0; is<rowsS; ++is){
            int zeroOne=0;
            if ( cover(C, SB, is, 1, A, k) > cover(C, SB, is, 0, A, k) ) zeroOne=1;
            s[is] = zeroOne;
        }
        double value=evaluate(C, SB, A, k, s);
        if (value > vmax){
            vmax=value;
            imax=k;
            smax=s;
        }
    }
}
Vec<Vec<int> > Alg1::product(int rowsS, int colsB){
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
double Alg1::cover(Vec<Vec<int> > &C, Vec<Vec<int> > &SB, const int &is,
             int s, Vec<Vec<int> > &A, int &k) {
    const int ncols=(int)C[0].size();
    int countp=0, countm=0;
    for ( int j=0; j<ncols; ++j){
        int sob = SB[is][j] + s * A[k][j];
        if ( C[is][j] == 1 && sob > 0 ) ++countp;
        if ( C[is][j] == 0 && sob > 0 ) ++countm;
    }
    return wp*countp - wm*countm;
}
double Alg1::evaluate(Vec<Vec<int> > &C, Vec<Vec<int> > &SB, Vec<Vec<int> > &A,
                const int & k, Vec<int> &s){
    const int nrows=(int)C.size();
    const int ncols=(int)C[0].size();
    int countp=0, countm=0;
    for (int i=0; i<nrows; ++i){
        for (int j=0; j<ncols; ++j){
            int sob=SB[i][j]+s[i]*A[k][j];
            if ( C[i][j] == 1 && sob > 0 ) ++countp;
            if ( C[i][j] == 0 && sob > 0 ) ++countm;
        }
    }
    return wp*countp - wm*countm;
}

