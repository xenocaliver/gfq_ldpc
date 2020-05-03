/* 
 * This file is part of the gfq_ldpc distribution (https://github.com/xenocaliver).
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
#include <list>
#include <cstdint>
#include <iomanip>
#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

std::vector<std::vector<Galois::Element> > gfq_matrix_product(std::vector<std::vector<Galois::Element> >& A, std::vector<std::vector<Galois::Element> >& B, const Galois::Field* gf) {
    uint64_t uli, ulj, ulk;
    std::vector<uint64_t> row_size, column_size;
    Galois::Element zero(gf, 0);
    Galois::Element sum = zero;
    std::vector<Galois::Element> v;
    std::vector<std::vector<Galois::Element> > C;

    row_size.resize(2);
    column_size.resize(2);
    column_size[0] = A[0].size();
    row_size[0]=A.size();
    column_size[1] = B[0].size();
    row_size[1] = B.size();
    if(column_size[0] != row_size[1]) {
        std::cerr << "Can not match dimension:" << row_size[0] << ", " << column_size[0] << ", " << row_size[1] << ", " << column_size[1] << std::endl;
        exit(1);
    }
    v.clear();
    for(ulj = 0; ulj < column_size[1]; ulj++) v.push_back(zero);
    for(uli = 0; uli < row_size[0]; uli++) C.push_back(v);

    for(uli = 0; uli < row_size[0]; uli++) {
        for(ulj = 0; ulj < column_size[1]; ulj++) {
            sum = zero;
            for(ulk = 0; ulk < column_size[0]; ulk++) sum += A[uli][ulk]*B[ulk][ulj];
            C[uli][ulj] = sum;
        }
    }
    return(C);
}

void clear_matrix(std::vector<std::vector<Galois::Element> >& M, const Galois::Field* gf) {
    uint64_t uli, ulj;
    Galois::Element zero(gf, 0);
    uint64_t row_size = M.size();
    uint64_t column_size = M[0].size();
    
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < column_size; ulj++) {
            M[uli][ulj] = zero;
        }
    }
}

bool make_matrix_full_rank(std::vector<std::vector<Galois::Element> >& parity_check_matrix, std::list<std::vector<std::vector<Galois::Element> > > list_of_Q, const Galois::Field* gf) {
    uint64_t uli, ulj, ulk;
    uint64_t row_size = parity_check_matrix.size();
    uint64_t column_size = parity_check_matrix[0].size();
    Galois::Element zero(gf, 0);
    Galois::Element one(gf, 1);
    std::list<std::vector<std::vector<Galois::Element> > > permutation_matrices;
    std::vector<Galois::Element> v, tmp;
    std::vector<std::vector<Galois::Element> > Q;
    bool found_flag = false;

    /* preparation */
    for(uli = 0; uli < row_size; uli++) tmp.push_back(zero);
    for(ulj = 0; ulj < column_size; ulj++) v.push_back(zero);
    for(ulj = 0; ulj < column_size; ulj++) Q.push_back(v);

    for(uli = 0; uli < row_size; uli++) {
        if(parity_check_matrix[uli][uli] != zero) continue;
        found_flag = false;
        for(ulj = uli + 1; ulj < column_size; ulj++) {
            if(parity_check_matrix[uli][ulj] != zero) { /* if uli-th row ulj-column element does not equal to zero */
                /* column swap */
                for(ulk = 0; ulk < row_size; ulk++) { 
                    tmp[ulk] = parity_check_matrix[ulk][uli];
                }
                for(ulk = 0; ulk < row_size; ulk++) { 
                    parity_check_matrix[ulk][uli] = parity_check_matrix[ulk][ulj];
                }
                for(ulk = 0; ulk < row_size; ulk++) { 
                    parity_check_matrix[ulk][ulj] = tmp[ulk];
                }
                clear_matrix(Q, gf);
                for(ulk = 0; ulk < uli; ulk++) Q[ulk][ulk] = one;
                Q[uli][ulj] = one;
                for(ulk = uli + 1; ulk < ulj; ulk++) Q[ulk][ulk] = one;
                Q[ulj][uli] = one;
                for(ulk = ulj + 1; ulk < column_size; ulk++) Q[ulk][ulk] = one;
                permutation_matrices.push_back(Q);
                found_flag = true;
                break;
            }
        }
        if(found_flag == false) {
            std::cerr << "This matrix may be singular." << std::endl;
            return(found_flag);
        }
    }
    list_of_Q = permutation_matrices;
    return(true);
}
std::vector<std::vector<Galois::Element> > make_generating_matrix(std::vector<std::vector<Galois::Element> >& parity_check_matrix, const Galois::Field* gf) {
    std::vector<std::vector<Galois::Element> > A;
    std::vector<std::vector<Galois::Element> > B;
    std::vector<std::vector<Galois::Element> > I;
    std::vector<std::vector<Galois::Element> > X, Y, Z, P;
    std::vector<std::vector<Galois::Element> > G, Gdash;
    std::vector<std::vector<Galois::Element> > full_rank_matrix;
    std::list<std::vector<std::vector<Galois::Element> > > list_of_Q;
    std::vector<Galois::Element> v;
    Galois::Element sum(gf, 0);
    uint64_t uli, ulj, ulk;
    int64_t i;
    uint64_t row_size;
    uint64_t column_size;
    bool full_rank;

    Galois::Element zero(gf, 0);
    Galois::Element one(gf, 1);

    row_size = parity_check_matrix.size();
    column_size = parity_check_matrix[0].size();
    std::cout << "*** parity check matrix ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < column_size; ulj++) {
            std::cout << std::setw(2) << parity_check_matrix[uli][ulj].value() << " ";
        }
        std::cout << std::endl;
    }

    v.clear();
    for(uli = 0; uli < parity_check_matrix[0].size(); uli++) v.push_back(zero);
    for(uli = 0; uli < parity_check_matrix.size(); uli++) full_rank_matrix.push_back(v);
    copy(parity_check_matrix.begin(), parity_check_matrix.end(), full_rank_matrix.begin());
    full_rank = make_matrix_full_rank(full_rank_matrix, list_of_Q, gf);
    if(full_rank == false) {
        std::cerr << "This matrix may be singular." << std::endl;
        exit(-1);
    }
    std::cout << "*** full rank parity check matrix ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < column_size; ulj++) {
            std::cout << std::setw(2) << full_rank_matrix[uli][ulj].value() << " ";
        }
        std::cout << std::endl;
    }
    v.clear();
    for(uli = 0; uli < row_size; uli++) v.push_back(zero);
    for(uli = 0; uli < row_size; uli++) A.push_back(v);
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            A[uli][ulj] = full_rank_matrix[uli][ulj];
        }
    }
    v.clear();
    /* create residual part matrix B */
    for(uli = 0; uli < column_size - row_size; uli++) v.push_back(zero);
    for(uli = 0; uli < row_size; uli++) B.push_back(v);
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < column_size - row_size; ulj++) {
            B[uli][ulj] = full_rank_matrix[uli][ulj + row_size];
        }
    }
    v.clear();
    /* initiaize L and U */
    for(uli = 0; uli < row_size; uli++) {
        v.push_back(zero);
    }
    for(uli = 0; uli < row_size; uli++) {
        X.push_back(v);
        I.push_back(v);
    }
    
    for(uli = 0; uli < row_size; uli++) {
        I[uli][uli] = one;
    }

    std::cout << "*** A ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << A[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }

    /* LU decomposition */
    for(uli = 0; uli < row_size; ++uli){
        // calculating L (i <= j)
        for(ulj = 0; ulj <= uli; ++ulj){
            sum = A[uli][ulj];
            for(ulk = 0; ulk < ulj; ++ulk){
                sum -= A[uli][ulk]*A[ulk][ulj];    // l_ik * u_kj
            }
            A[uli][ulj] = sum;
        }
 
        // calculating U (i < j)
        for(ulj = uli + 1; ulj < row_size; ++ulj){
            sum = A[uli][ulj];
            for(ulk = 0; ulk < uli; ++ulk){
                sum -= A[uli][ulk]*A[ulk][ulj];    // l_ik * u_kj
            }
            A[uli][ulj] = sum/A[uli][uli];
        }
    }
    std::cout << "LU decomposition completed." << std::endl;
    std::cout << "*** LU ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << A[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
    /* forward substitution */
    for(ulj = 0; ulj < row_size; ulj++) {
        for(uli = 0; uli < row_size; uli++) {
            sum = I[uli][ulj];
            for(ulk = 0; ulk < uli; ulk++) sum -= A[uli][ulk]*X[ulk][ulj];
            X[uli][ulj] = sum/A[uli][uli];
        }
    }
    std::cout << "*** Y ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << X[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }

    /* backward substitution */
    for(ulj = 0; ulj < row_size; ulj++) {
        for(i = row_size - 1; i >= 0; i--) {
            sum = X[i][ulj];
            for(ulk = i + 1; ulk < row_size; ulk++) {
                sum -= A[i][ulk]*X[ulk][ulj];
            }
            X[i][ulj] = sum;
        }
    }
    std::cout << "*** X ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << X[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }

    v.clear();
    A.clear();
    for(uli = 0; uli < row_size; uli++) v.push_back(zero);
    for(uli = 0; uli < row_size; uli++) A.push_back(v);
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            A[uli][ulj] = full_rank_matrix[uli][ulj];
        }
    }

    Y = gfq_matrix_product(X, A, gf);
    std::cout << "*** Inverse matrix check result ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << Y[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Inverse matrix calculation completed." << std::endl;
    G = gfq_matrix_product(X, B, gf);

    /* create complete generating matrix */
    v.clear();
    Z.clear();
    std::cout << "row_size = " << row_size << " column_size = " << column_size << std::endl;
    for(uli = 0; uli < row_size; uli++) v.push_back(zero);
    for(uli = 0; uli < column_size; uli++) Z.push_back(v);

    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            Z[uli][ulj] = G[uli][ulj];
        }
    }
    for(uli = 0; uli < row_size; uli++) {
        Z[uli + row_size][uli] = one;
    } 

    std::cout << "*** transeposed generating matrix ***" << std::endl;
    for(uli = 0; uli < column_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << Z[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }

    X = gfq_matrix_product(full_rank_matrix, Z, gf);
    std::cout << "*** check result ***" << std::endl;
    for(uli = 0; uli < X.size(); uli++) {
        for(ulj = 0; ulj < X[0].size(); ulj++) {
            std::cout << std::setw(2) << X[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }

    return(G);
}