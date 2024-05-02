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
#ifndef VARIABLE_NODE_HPP_
#define VARIABLE_NODE_HPP_
#include <cstdint>
#include <vector>

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

#include "edge.hpp"
#include "constants.hpp"

class variable_node {                                  /*                     */
public:                                                /*                     */
	uint64_t ID;                                       /* identifier          */
	uint64_t value;                           /* galois field element's value */
	std::vector<edge*> edges;                          /* edges               */
    std::vector<double> channel_output_probability;    /*                     */
	variable_node(void) : ID(0), value(0){}            /* default constructor */
	variable_node(uint64_t id) : ID(id), value(0) {}   /* standard constructor*/
	~variable_node(void) {}                            /* destructor          */
                                                       /*                     */
    void update_messages(const Galois::Field* gf) {
        uint64_t uli, ulj;
        uint64_t characteristic = (uint64_t)(gf->q);
        uint64_t g;                             /* galois field element value */
        double product;                                /*                     */
        double sum;
        std::vector<double> new_message(characteristic, 1.0);

        for(uli = 0; uli < this->edges.size(); uli++) {
            for(g = 0; g < characteristic; g++) {
                product = this->channel_output_probability[g];
                for(ulj = uli + 1; ulj < this->edges.size(); ulj++) {
                    product *= this->edges[ulj]->factor_to_variable_message[g];
                }
                for(ulj = 0; ulj < uli; ulj++) {
                    product *= this->edges[ulj]->factor_to_variable_message[g];
                }
                new_message[g] = product;
            }
            sum = 0.0;
            for(g = 0; g < characteristic; g++) {
                sum += new_message[g];
            }
            for(g = 0; g < characteristic; g++) {
                new_message[g] /= sum;
            }
            this->edges[uli]->variable_to_factor_message = new_message;
        }
    }
    uint64_t speculate_temporal_symbol(const Galois::Field* gf) {
        uint64_t ulj;
        uint64_t characteristic = (uint64_t)(gf->q);
        uint64_t g;                             /* galois field element value */
        double product;                                /*                     */
        double sum = 0.0;
        double max = -1000.0;
        uint64_t rtnv = 0;
        std::vector<double> new_message(characteristic, 1.0);

        for(g = 0; g < characteristic; g++) {
            product = this->channel_output_probability[g];
            for(ulj = 0; ulj < this->edges.size(); ulj++) {
                if(this->edges[ulj]->factor_to_variable_message[g] > lower_limit) {
                    product *= this->edges[ulj]->factor_to_variable_message[g];
                } else {
                    product *= lower_limit;
                }
            }
            new_message[g] = product;
        }

        /* normalize messages */
        for(g = 0; g < characteristic; g++) {
            sum += new_message[g];
        }
        
        for(g = 0; g < characteristic; g++) {
            new_message[g] = new_message[g]/sum;
        }

        /* find argmax */
        for(g = 0; g < characteristic; g++) {
            if(new_message[g] > max) {
                max = new_message[g];
                rtnv = g;
            }
        }
        return(rtnv);
    }
};                                                     /*                     */
#endif