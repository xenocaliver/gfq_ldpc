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

#include <cstdlib>
#include <cstdint>
#include <vector>

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

#include "variable_node.hpp"
#include "factor_node.hpp"
#include "edge.hpp"
#include "gfq_alist.hpp"

std::vector<Galois::Element> split_bitpattern(uint64_t bit_pattern, uint64_t extend_degree, uint64_t number_of_variables, const Galois::Field* gf) {
    uint64_t full_one = (1 << extend_degree) - 1;
    Galois::Element zero(gf, 0);
    Galois::Element gfe = zero;
    std::vector<Galois::Element> rtnv;
    uint64_t x;

    for(uint64_t uli = 0; uli < number_of_variables; uli++) {
        x = bit_pattern & full_one;
        gfe.setValue(x);
        rtnv.push_back(gfe);
        bit_pattern = bit_pattern >> extend_degree;
    }
    return(rtnv);
}

std::vector<std::vector<uint64_t> > search_for_fullfill_set(std::vector<std::pair<uint64_t, uint64_t> >& row_elements, const Galois::Field* gf) {
    Galois::Element zero(gf, 0);
    Galois::Element sum = zero;
    std::vector<Galois::Element> parity_check_elements;
    std::vector<Galois::Element> random_variables;
    uint64_t extend_degree = (uint64_t)(log2(gf->q));
    uint64_t extend_max = 1 << (extend_degree*row_elements.size());
    uint64_t x;
    std::vector<uint64_t> v;
    std::vector<std::vector<uint64_t> > rtnv;

    for(uint64_t uli = 0; uli < row_elements.size(); uli++) {
        parity_check_elements.push_back(zero);
        random_variables.push_back(zero);
    }
    for(uint64_t uli = 0; uli < parity_check_elements.size(); uli++) {
        parity_check_elements[uli].setValue(row_elements[uli].second);
    }

    for(x = 1; x < extend_max; x++) {
        random_variables = split_bitpattern(x, extend_degree, row_elements.size(), gf);
        sum = zero;
        for(uint64_t uli = 0; uli < random_variables.size(); uli++) {
            sum += parity_check_elements[uli]*random_variables[uli];
        }
        if(sum == zero) {
            v.clear();
            for(uint64_t uli = 0; uli < random_variables.size(); uli++) {
                v.push_back(random_variables[uli].value());
            }
            rtnv.push_back(v);
        }
    }
    return(rtnv);
}

void construct_factor_graph(std::vector<variable_node>& variable_nodes, std::vector<factor_node>& factor_nodes, std::vector<edge>& edges, gfq_alist alist, Galois::Field* gf) {
    uint64_t number_of_columns = alist.number_of_columns;
    uint64_t number_of_rows = alist.number_of_rows;
    std::vector<uint64_t> row_weights = alist.mwlist;
    uint64_t uli, ulj;
    uint64_t total_number_of_edges = 0;
    uint64_t edge_id = 0;
    std::vector<std::vector<uint64_t> > fullfill_table;

    for(uli = 0; uli < row_weights.size(); uli++) {
        total_number_of_edges += row_weights[uli];
    }
    variable_nodes.resize(number_of_columns);
    factor_nodes.resize(number_of_rows);
    edges.resize(total_number_of_edges);
    for(uli = 0; uli < variable_nodes.size(); uli++) variable_nodes[uli].ID = uli;
    for(uli = 0; uli < factor_nodes.size(); uli++) factor_nodes[uli].ID = uli;
    for(uli = 0; uli < total_number_of_edges; uli++) edges[uli].ID = uli;
    for(uli = 0; uli < number_of_rows; uli++) factor_nodes[uli].parity_check_matrix_elements = alist.mlist[uli];
    
    for(uli = 0; uli < number_of_rows; uli++) {
        for(ulj = 0; ulj < alist.mlist[uli].size(); ulj++) {
            factor_nodes[uli].edges.push_back(&(edges[edge_id]));
            variable_nodes[alist.mlist[uli][ulj].first - 1].edges.push_back(&(edges[edge_id]));
            edge_id++;
        }
        factor_nodes[uli].parity_check_matrix_elements = alist.mlist[uli];
        fullfill_table = search_for_fullfill_set(alist.mlist[uli], gf);
        factor_nodes[uli].fullfill_table = fullfill_table;
    }
}