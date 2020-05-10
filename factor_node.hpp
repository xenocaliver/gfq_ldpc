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

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

#include "edge.hpp"


class factor_node {                                    /*                     */
public:                                                /*                     */
    uint64_t ID;                                       /* identifier          */
    std::vector<edge*> edges;                          /* edges               */
    std::vector<std::pair<uint64_t, uint64_t> > parity_check_matrix_elements;
    std::vector<std::vector<uint64_t> > fullfill_table; /*                    */
    factor_node(void) : ID(0) {}                       /* default constructor */

    std::vector<std::vector<uint64_t> > search_for_other_elements(uint64_t uli, uint64_t g) {
        std::vector<std::vector<uint64_t> > rtnv;
        std::vector<uint64_t> u, v;
        std::vector<std::vector<uint64_t> >::iterator vit;

        for(vit = this->fullfill_table.begin(); vit != this->fullfill_table.end(); ++vit) {
            u = *vit;
            if(u[uli] == g) {
                v = *vit;
                rtnv.push_back(v);
            }
        }
        return(rtnv);
    }

    void update_messages(const Galois::Field* gf) {
        uint64_t uli, ulj, ull;
        uint64_t g;
        double sum = 1.0;
        double prod;
        std::vector<std::vector<uint64_t> > nonzero_combination;

        for(uli = 0; uli < edges.size(); uli++) {
            for(g = 0; g < (uint64_t)(gf->q); g++) {
                sum = 0.0;
                nonzero_combination = this->search_for_other_elements(uli, g);
                for(ull = 0; ull < nonzero_combination.size(); ull++) {
                    prod = 1.0;
                    for(ulj = uli + 1; ulj < edges.size(); ulj++) {
                        prod *= this->edges[ulj]->variable_to_factor_message[nonzero_combination[ull][ulj]];
                    }
                    for(ulj = 0; ulj < uli; uli++) {
                        prod *= this->edges[ulj]->variable_to_factor_message[nonzero_combination[ull][ulj]];
                    }
                    sum += prod;
                }
                this->edges[uli]->factor_to_variable_message[g] = sum;
            }
        }
    }
};
#endif /* SRC_FACTORNODE_HPP_ */
