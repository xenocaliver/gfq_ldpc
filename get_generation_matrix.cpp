/* 
 * This file is part of the GFq_ludecomp distribution (https://github.com/xenocaliver).
 * Copyright (c) 2020 Akiyoshi Hashimoto.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

std::vector<std::vector<Galois::Element> > gfq_matrix_product(std::vector<std::vector<Galois::Element> >& A, std::vector<std::vector<Galois::Element> >& B, Galois::Field& gf) {
    uint64_t uli, ulj, ulk;
    std::vector<uint64_t> row_size, column_size;
    Galois::Element zero(&gf, 0);
    std::vector<Galois::Element> v;
    std::vector<std::vector<Galois::Elememnt> > C;

    row_size.resize(2);
    column_size.resize(2)
    column_size[0] = A.size();
    row_size[0]=A[0].size();
    column_size[1] = B.size();
    row_size[1] = B[0].size();
    if(column_size[0] != row_size[1]) {
        std::cerr << "Can not match dimension:" << row_size[0] << ", " << column_size[0] << ", " << row_size[1] << ", " << column_size[1] << std::endl;
        exit(1)
    }
    for(ulj = 0; ulj < column_size[1]; ulj++) v.push_back(zero);
    for(uli = 0; uli < row_size[0]; uli++) C.push_back(v);

    for(uli = 0; uli < row_size[0]; uli++) {
        for(ulj = 0; ulj < column_size[1]; ulj++) {
            for(ulk = 0; ulk < column_size[0]) C[uli][ulj] = A[uli][ulk]*B[ulk][ulj];
        }
    }
    return(C);
}

std::vector<std::vector<Galois::Element> > make_generating_matrix(std::vector<std::vector<Galois::Element> >& parity_check_matrix, Galois::Field gf) {
    std::vector<std::vector<Galois::Element> > A;
    std::vector<std::vector<Galois::Element> > B;
    std::vector<std::vector<Galois::Element> > L;
    std::vector<std::vector<Galois::Element> > U;
    std::vector<std::vector<Galois::Element> > I;
    std::vector<std::vector<Galois::Element> > X, Y;
    std::vector<std::vector<Galois::Element> > G;
    std::vector<Galois::Element> v;
    uint64_t uli, ulj, ulk;
    uint64_t row_size;
    Galois::Element sum;

    Galois::Element zero(&gf, 0);
    Galois::Element one(&gf, 1);

    row_size = parity_check_matrix[0];

    /* create square matrix A */
    for(uli = 0; uli < row_size; uli++) A.push_back(parity_check_matrix[uli])
    /* create residual part matrix B */
    for(uli = row_size; uli < parity_check_matrix.size(); uli++) B.push_back(parity_check_matrix[uli]);
    /* initialize L and U */
    for(uli = 0; uli < row_size; uli++) {
        v.push_back(zero);
    }
    for(uli = 0; uli < row_size; uli++) {
        L.push_back(v);
    }

    std::copy(L.begin(), L.end(), U.begin());
    std::copy(L.begin(), L.end(), X.begin());
    std::copy(L.begin(), L.end(), Y.begin());

    for(uli = 0; uli < row_size; uli++) {
        I[uli][uli] = one;
    }

    copy(I.begin(), I.end(), L.begin());

    /* LU decomposition */
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            if(uli <= ulj){
                sum = zero;
                for(ulk = 0; ulk < uli - 1; uli++) {
                    sum = sum + L[uli][ulk]*U[ulk][ulj];
                }
                U[uli][ulj] = A[uli][ulj] - sum;
            } else {
                sum = zero;
                for(ulk = 0; ulk < ulj - 1; ulk++) {
                    sum = sum + L[uli][ulk]*U[ulk][ulj];
                }
                L[uli][ulj] = (A[uli][ulj] - sum)/U[ulj][ulj];
            }
        }
    }

    /* Now, get inverse of matrix A */
    /* forward substitution */
    copy(I.begin(), I.end(), Y.begin());
    for(ulj = 0; ulj < row_size; ulj++) {
        for(ulk = 0; ulk < row_size - 1; ulk++) {
            for(uli = ulk + 1; uli < row_size; ulk++) {
                Y[uli][ulj] = Y[uli][ulj] - Y[ulk][ulj]*L[uli][ulk];
            }
        }
    }

    /* backword substitution */
    copy(Y.begin(), Y.end(), X.start());
    for(ulj = 0; ulj < row_size; ulj++) {
        for(ulk = row_size - 1; ulk >=0; ulk--) {
            X[ulk][ulj] = X[ulk][ulj]/U[ulk][ulk];
            for(uli = 0; uli < ulk; uli++){
                X[uli][ulj] = X[uli][ulj] - U[uli][ulk]*X[ulk][ulj]
            }
        }
    }
    /* Now, get generating matrix of parity check matrix */
    G = gfq_matrix_product(X, B, gf);
    return(G);
}