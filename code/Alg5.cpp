// Copyright (C) 2019, International Business Machines
// Corporation.  All Rights Reserved.

// This program is distributed under the terms of the
// Eclipse Public License - v 2.0

#include "Alg5.hpp"
#include <cassert>

Alg5::Alg5(Vec<Vec<int> > &C, int k, Vec<Vec<int> > &Stranp, Vec<Vec<int> > &B){
    assert(Stranp.size()==k);
    assert(B.size()==k);
    wp=wm=1.0;
    const int n=(int)C.size();
    Vec<Vec<int> > Stranp_new;
    Stranp_new.resize(k);
    for(int l=0; l<k; l++) {
        for(int i=0; i<(int)Stranp[l].size(); i++) {
            Stranp_new[l].push_back(Stranp[l][i]);
        }
    }
    
    int minError = overall_error(C, Stranp, B);
    int error;
    do {
        for (int l=0; l < k; ++l){
            Vec<int> s;
            find_max_alg5(B, Stranp_new, C, l, s);
            for(int i=0; i<n; i++)
                Stranp_new[l][i]=s[i];
        }
        error = overall_error(C, Stranp_new, B);
        //std::cout<<"error="<<error<<", minError="<<minError<<std::endl;
        if(error < minError) {
            minError=error;
            for (int l=0; l < k; ++l) {
                for(int i=0; i<n; i++) {
                    Stranp[l][i] = Stranp_new[l][i];
                }
            }
        }
        else
            break;
    } while(true);
}
int Alg5::overall_error(Vec<Vec<int> > &C, Vec<Vec<int> > &Stranp, Vec<Vec<int> > &B){
    int rowsS=(int)Stranp[0].size();
    int rowsB=(int)B.size();
    int colsB=(int)B[0].size();
    int error=0;
    for (int irow=0; irow<rowsS; ++irow){
        for (int icol=0; icol<colsB; ++icol){
            int sum=0;
            for (int l=0; l<rowsB; ++l) {
                sum+=Stranp[l][irow]*B[l][icol];
            }
            if((sum>0 && C[irow][icol]==0) || (sum==0 && C[irow][icol]==1))
                error++;
        }
    }
    return error;
}
void Alg5::find_max_alg5(Vec<Vec<int> > &B, Vec<Vec<int> > &Stranp,
                   Vec<Vec<int> > &C, int l, Vec<int> &smax)
{
    const int rowsS=(int)C.size();
    assert(Stranp.size()>0);
    assert(B.size()>0);
    Vec<Vec<int> > SB = product_alg5(Stranp, B, l);
    assert(SB.size()>0);
    Vec<int> s(rowsS,-1);
    for (int is=0; is<rowsS; ++is){
        int zeroOne=0;
        if ( cover_alg5(C, SB, is, 1, B, l) > cover_alg5(C, SB, is, 0, B, l) ) zeroOne=1;
        s[is] = zeroOne;
    }
    smax=s;
}
Vec<Vec<int> > Alg5::product_alg5(Vec<Vec<int> > &Stranp, Vec<Vec<int> > &B, int colS){
    Vec<Vec<int> > SB;
    int rowsS=(int)Stranp[0].size();
    int rowsB=(int)B.size();
    int colsB=(int)B[0].size();
    Vec<int> row(colsB);
    for (int irow=0; irow<rowsS; ++irow){
        for (int icol=0; icol<colsB; ++icol){
            int sum=0;
            for (int k=0; k<rowsB; ++k) {
                if(k != colS)
                    sum+=Stranp[k][irow]*B[k][icol];
            }
            row[icol]= (sum > 0) ? 1 : 0;
        }
        SB.push_back(row);
    }
    return SB;
}
double Alg5::cover_alg5(Vec<Vec<int> > &C, Vec<Vec<int> > &SB, const int &is,
                  int s, Vec<Vec<int> > &B, int l) {
    const int ncols=(int)C[0].size();
    int countp=0, countm=0;
    for ( int j=0; j<ncols; ++j){
        int sob = SB[is][j] + s * B[l][j];
        if ( C[is][j] == 1 && sob > 0 ) ++countp;
        if ( C[is][j] == 0 && sob > 0 ) ++countm;
    }
    return wp*countp - wm*countm;
}
