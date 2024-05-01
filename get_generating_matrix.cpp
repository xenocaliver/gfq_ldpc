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

uint64_t search_upper_triangle_vector(uint64_t n, std::vector<std::vector<Galois::Element> >& parity_check_matrix, Galois::Field* gf) {
    uint64_t uli, ulj;
    Galois::Element zero(gf, 0);
    uint64_t flag = 0;

    uint64_t row_size = parity_check_matrix.size();
    uint64_t column_size = parity_check_matrix[0].size();
    for(ulj = n; ulj < column_size; ulj++) {
        flag = 0;
        for(uli = n; uli < row_size; uli++) {
            if(parity_check_matrix[uli][n] != zero) {
                break;
            } else {
                flag++;
            }
        }
        if(flag == row_size - n) {
            return(ulj);
        } 
    }
}

uint64_t make_matrix_full_rank(std::vector<std::vector<Galois::Element> >& parity_check_matrix, std::list<std::vector<std::vector<Galois::Element> > >& list_of_Q, const Galois::Field* gf) {
    uint64_t uli, ulj, ulk;
    uint64_t row_size = parity_check_matrix.size();
    uint64_t column_size = parity_check_matrix[0].size();
    Galois::Element zero(gf, 0);
    Galois::Element one(gf, 1);
    std::list<std::vector<std::vector<Galois::Element> > > permutation_matrices;
    std::vector<Galois::Element> v, tmp;
    std::vector<std::vector<Galois::Element> > Q;
    uint64_t found_flag = 0;

    /* preparation */
    for(uli = 0; uli < row_size; uli++) tmp.push_back(zero);
    for(ulj = 0; ulj < column_size; ulj++) v.push_back(zero);
    for(ulj = 0; ulj < column_size; ulj++) Q.push_back(v);

    for(uli = 0; uli < row_size; uli++) {
        if(parity_check_matrix[uli][uli] != zero) {
            found_flag++;
            continue;
        }
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
                found_flag++;
                break;
            }
        }
    }
    list_of_Q = permutation_matrices;
    return(found_flag);
}

std::vector<std::vector<Galois::Element> > make_upper_triangluar_matrix(std::vector<std::vector<Galois::Element> >& full_rank_matrix, const Galois::Field* gf) {
    uint64_t uli, ulj;
    Galois::Element zero(gf, 0);
    std::vector<std::vector<Galois::Element> > rtnv;
    std::vector<std::vector<Galois::Element> > P;
    std::vector<Galois::Element> v;
    uint64_t row_size, column_size;

    // preparation
    row_size = full_rank_matrix.size();
    column_size = row_size;
    for(uli = 0; uli < column_size; uli++) {
        v.push_back(zero);
    }
    for(uli = 0; uli < row_size; uli++) {
        P.push_back(v);
        rtnv.push_back(v);
    }

    // search lower non zero element
    for(ulj = 0; ulj < column_size; ulj++) {
        for(uli = ulj + 1; uli < row_size; uli++) {
            if(full_rank_matrix[uli][ulj] != zero) {}
        }
    }
}

std::vector<std::vector<Galois::Element> > make_generating_matrix(std::vector<std::vector<Galois::Element> >& parity_check_matrix, const Galois::Field* gf) {
    std::vector<std::vector<Galois::Element> > T, A, B, C, D, E;
    std::vector<std::vector<Galois::Element> > I, X, Y, Z;
    std::vector<std::vector<Galois::Element> > G, Gdash, Q;
    std::vector<std::vector<Galois::Element> > full_rank_matrix;
    std::list<std::vector<std::vector<Galois::Element> > > list_of_Q;
    std::vector<Galois::Element> v;
    Galois::Element sum(gf, 0);
    uint64_t uli, ulj, ulk;
    int64_t i;
    uint64_t row_size;
    uint64_t column_size;
    uint64_t parity_check_size;
    uint64_t gap;
    uint64_t T_size;

    Galois::Element zero(gf, 0);
    Galois::Element one(gf, 1);

    row_size = parity_check_matrix.size();
    column_size = parity_check_matrix[0].size();
    parity_check_size = column_size - row_size;
#ifdef DEBUG
    std::cout << "*** parity check matrix ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < column_size; ulj++) {
            std::cout << std::setw(2) << parity_check_matrix[uli][ulj].value() << " ";
        }
        std::cout << std::endl;
    }
#endif

    v.clear();
    for(uli = 0; uli < parity_check_matrix[0].size(); uli++) v.push_back(zero);
    for(uli = 0; uli < parity_check_matrix.size(); uli++) full_rank_matrix.push_back(v);
    copy(parity_check_matrix.begin(), parity_check_matrix.end(), full_rank_matrix.begin());
    T_size = make_matrix_full_rank(full_rank_matrix, list_of_Q, gf);
    gap = parity_check_matrix.size() - T_size;
    std::cout << "gap = " << gap << std::endl;
#ifdef DEBUG
    std::cout << "*** full rank parity check matrix ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < column_size; ulj++) {
            std::cout << std::setw(2) << full_rank_matrix[uli][ulj].value() << " ";
        }
        std::cout << std::endl;
    }
#endif
    v.clear();
    for(uli = 0; uli < T_size; uli++) v.push_back(zero);
    for(uli = 0; uli < T_size; uli++) T.push_back(v);
    for(uli = 0; uli < T_size; uli++) {
        for(ulj = 0; ulj < T_size; ulj++) {
            T[uli][ulj] = full_rank_matrix[uli][ulj];
        }
    }
    v.clear();
    /* initiaize L and U */
    for(uli = 0; uli < T_size; uli++) {
        v.push_back(zero);
    }
    for(uli = 0; uli < T_size; uli++) {
        X.push_back(v);
        I.push_back(v);
    }
    
    for(uli = 0; uli < T_size; uli++) {
        I[uli][uli] = one;
    }

#ifdef DEBUG
    std::cout << "*** T ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << T[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
#endif
    /* LU decomposition */
    for(uli = 0; uli < T_size; ++uli){
        // calculating L (i <= j)
        for(ulj = 0; ulj <= uli; ++ulj){
            sum = T[uli][ulj];
            for(ulk = 0; ulk < ulj; ++ulk){
                sum -= T[uli][ulk]*T[ulk][ulj];    // l_ik * u_kj
            }
            T[uli][ulj] = sum;
        }
 
        // calculating U (i < j)
        for(ulj = uli + 1; ulj < T_size; ++ulj){
            sum = T[uli][ulj];
            for(ulk = 0; ulk < uli; ++ulk){
                sum -= T[uli][ulk]*T[ulk][ulj];    // l_ik * u_kj
            }
            T[uli][ulj] = sum/T[uli][uli];
        }
    }
#ifdef DEBUG
    std::cout << "LU decomposition completed." << std::endl;
    std::cout << "*** LU ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << T[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
#endif
    /* forward substitution */
    for(ulj = 0; ulj < T_size; ulj++) {
        for(uli = 0; uli < T_size; uli++) {
            sum = I[uli][ulj];
            for(ulk = 0; ulk < uli; ulk++) sum -= T[uli][ulk]*X[ulk][ulj];
            X[uli][ulj] = sum/T[uli][uli];
        }
    }
#ifdef DEBUG
    std::cout << "*** X ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << X[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
#endif

    /* backward substitution */
    for(ulj = 0; ulj < T_size; ulj++) {
        for(i = T_size - 1; i >= 0; i--) {
            sum = X[i][ulj];
            for(ulk = i + 1; ulk < row_size; ulk++) {
                sum -= T[i][ulk]*X[ulk][ulj];
            }
            X[i][ulj] = sum;
        }
    }
#ifdef DEBUG
    std::cout << "*** X ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << X[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
#endif

    v.clear();
    A.clear();
    for(uli = 0; uli < T_size; uli++) v.push_back(zero);
    for(uli = 0; uli < T_size; uli++) A.push_back(v);
    for(uli = 0; uli < T_size; uli++) {
        for(ulj = 0; ulj < T_size; ulj++) {
            A[uli][ulj] = full_rank_matrix[uli][ulj];
        }
    }

    Y = gfq_matrix_product(X, A, gf);
    for(uli = 0; uli < T_size; uli++) {
        for(ulj = 0; ulj < T_size; ulj++) {
            if(Y[uli][ulj] != I[uli][ulj]) {
                std::cerr << "Inverse matrix is not correct." << std::endl;
                exit(-1);
            }
        }
    }
#ifdef DEBUG
    std::cout << "*** Inverse matrix check result ***" << std::endl;
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << Y[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Inverse matrix calculation completed." << std::endl;
#endif

    v.clear();
#ifdef DEBUG
    std::cout << "parity_check_size = " << parity_check_size << std::endl;
#endif
    for(uli = 0; uli < parity_check_size; uli++) v.push_back(zero);
    for(uli = 0; uli < row_size; uli++) B.push_back(v);
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < parity_check_size; ulj++) {
            B[uli][ulj] = full_rank_matrix[uli][T_size + ulj];
        }
    }
    G = gfq_matrix_product(X, B, gf);

    /* create complete generating matrix */
    v.clear();
    Z.clear();
#ifdef DEBUG
    std::cout << "row_size = " << row_size << " column_size = " << column_size << std::endl;
#endif
    for(uli = 0; uli < row_size; uli++) v.push_back(zero);
    for(uli = 0; uli < column_size; uli++) Z.push_back(v);

    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            Z[uli][ulj] = -G[uli][ulj];
        }
    }
    for(uli = 0; uli < row_size; uli++) {
        Z[uli + row_size][uli] = one;
    } 

#ifdef DEBUG
    std::cout << "*** transeposed generating matrix ***" << std::endl;
    for(uli = 0; uli < column_size; uli++) {
        for(ulj = 0; ulj < row_size; ulj++) {
            std::cout << std::setw(2) << Z[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
#endif

    X = gfq_matrix_product(full_rank_matrix, Z, gf);
#ifdef DEBUG
    std::cout << "*** check result ***" << std::endl;
    for(uli = 0; uli < X.size(); uli++) {
        for(ulj = 0; ulj < X[0].size(); ulj++) {
            std::cout << std::setw(2) << X[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
#endif

    /* column permutation matrix */
    v.clear();
    I.clear();
    for(uli = 0; uli < column_size; uli++) v.push_back(zero);
    for(uli = 0; uli < column_size; uli++) I.push_back(v);
    for(uli = 0; uli < column_size; uli++) I[uli][uli] = one;
    std::list<std::vector<std::vector<Galois::Element> > >::iterator qit;
    Q = I;
    for(qit = list_of_Q.begin(); qit != list_of_Q.end(); ++qit) {
#ifdef DEBUG
        X = *qit;
        std::cout << "*** permutation matrix ***" << std::endl;
#endif
        Q = gfq_matrix_product(Q, *qit, gf);
#ifdef DEBUG
        for(uli = 0; uli < X.size(); uli++) {
            for(ulj = 0; ulj < X[0].size(); ulj++) {
                std::cout << std::setw(2) << X[uli][ulj] << " ";
            }
            std::cout << std::endl;
        }
#endif
    }

    G = gfq_matrix_product(Q, Z, gf);

    /* final check */
    X = gfq_matrix_product(parity_check_matrix, G, gf);
#ifdef DEBUG
    std::cout << "*** final check result ***" << std::endl;
    for(uli = 0; uli < X.size(); uli++) {
        for(ulj = 0; ulj < X[0].size(); ulj++) {
            std::cout << std::setw(2) << X[uli][ulj] << " ";
        }
        std::cout << std::endl;
    }
#endif
    for(uli = 0; uli < X.size(); uli++) {
        for(ulj = 0; ulj < X[0].size(); ulj++) {
            if(X[uli][ulj] != zero) {
                std::cerr << "Incorrect generating matrix." << std::endl;
                exit(-1);
            }
        }
    }

    return(G);
}