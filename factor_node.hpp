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
#ifndef FACTORNODE_HPP_
#define FACTORNODE_HPP_
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <iomanip>

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

#include <fftw3.h>

#include "edge.hpp"
#include "constants.hpp"


class factor_node {                                    /*                     */
public:                                                /*                     */
    uint64_t ID;                                       /* identifier          */
    std::vector<edge*> edges;                          /* edges               */
    std::vector<std::pair<uint64_t, uint64_t> > parity_check_matrix_elements;
    factor_node(void) : ID(0) {}                       /* default constructor */

    std::vector<double> fourier_transform(std::vector<double> q) {
        std::vector<double> Q(q.size(), 0.0);

        fftw_make_planner_thread_safe();
        fftw_plan plan = fftw_plan_r2r_1d(q.size(), q.data(), Q.data(), FFTW_REDFT11, FFTW_MEASURE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        return(Q);
    }

    void update_messages(const Galois::Field* gf) {
        uint64_t uli, ulj;
        uint64_t g;
        std::vector<std::vector<double> > Q(this->edges.size(), std::vector<double>((uint64_t)(gf->q), 0.0));
        std::vector<double> R((uint64_t)(gf->q), 0.0);
        std::vector<double> product((uint64_t)(gf->q), 1.0);

        for(uli = 0; uli < this->edges.size(); uli++) {
            Q[uli] = fourier_transform(this->edges[uli]->variable_to_factor_message);
        }

        for(uli = 0; uli < this->edges.size(); uli++) {
            for(g = 0; g < (uint64_t)(gf->q); g++) {
                for(ulj = uli + 1; ulj < this->edges.size(); ulj++) {
                    product[g] *= Q[ulj][g];
                }
                for(ulj = 0; ulj < uli; ulj++) {
                    product[g] *= Q[ulj][g];
                }
            }
            R = fourier_transform(product);
            for(g = 0; g < (uint64_t)(gf->q); g++) R[g] = R[g]/(double)(2*(uint64_t)(gf->q)); /* renormalize for inverse fourier transformation */
            this->edges[uli]->factor_to_variable_message = R;
            for(g = 0; g < (uint64_t)(gf->q); g++) product[g] = 1.0;
        }
    }
};
#endif /* SRC_FACTORNODE_HPP_ */
