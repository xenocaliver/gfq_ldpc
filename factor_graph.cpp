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
#include <algorithm>

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

#include "variable_node.hpp"
#include "factor_node.hpp"
#include "edge.hpp"
#include "gfq_alist.hpp"


void construct_factor_graph(std::vector<variable_node>& variable_nodes, std::vector<factor_node>& factor_nodes, std::vector<edge>& edges, gfq_alist alist, Galois::Field* gf) {
    uint64_t number_of_columns = alist.number_of_columns;
    uint64_t number_of_rows = alist.number_of_rows;
    std::vector<uint64_t> row_weights = alist.mwlist;
    uint64_t uli, ulj;
    uint64_t total_number_of_edges = 0;
    uint64_t edge_id = 0;
    std::vector<std::vector<uint64_t> > u;
    std::vector<uint64_t> v;
    Galois::Element sum(gf, 0);
    Galois::Element zero(gf, 0);
    Galois::Element gfe(gf, 0);
    Galois::Element parity(gf, 0);

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
    }
}