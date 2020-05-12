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


std::vector<std::vector<uint64_t> > search_for_fullfill_set(std::vector<std::pair<uint64_t, uint64_t> >& row_elements, const Galois::Field* gf) {
    uint64_t uli, ulj;
    Galois::Element zero(gf, 0);
    Galois::Element sum = zero;
    std::vector<Galois::Element> parity_check_elements;
    std::vector<Galois::Element> random_variables;
    std::vector<uint64_t> v;
    std::vector<std::vector<uint64_t> > rtnv;
    uint64_t depth = row_elements.size();
    std::vector<uint64_t> slots(depth, 0);
    uint64_t index;
    uint64_t max = (uint64_t)(gf->q);

    v.resize(depth);
    for(uint64_t uli = 0; uli < row_elements.size(); uli++) {
        parity_check_elements.push_back(zero);
        random_variables.push_back(zero);
    }
    for(uint64_t uli = 0; uli < parity_check_elements.size(); uli++) {
        parity_check_elements[uli].setValue(row_elements[uli].second);
    }

    index = 0;
    while(true) {
        sum = zero;
        for(ulj = 0; ulj < depth; ulj++) {
            random_variables[ulj].setValue(slots[ulj]);
            sum += parity_check_elements[ulj]*random_variables[ulj];
        }
        if(sum == zero) {
            for(uli = 0; uli < depth; uli++) v[uli] = random_variables[uli].value();
            rtnv.push_back(v);
        }
        slots[0]++;
        while(slots[index] == max) {
            if(index == depth - 1) {
                rtnv.erase(rtnv.begin());    // 0, 0, 0,... case is invalid
                return(rtnv);
            }
            slots[index++] = 0;
            slots[index]++;
        }
        index = 0;
    }

}

void construct_factor_graph(std::vector<variable_node>& variable_nodes, std::vector<factor_node>& factor_nodes, std::vector<edge>& edges, gfq_alist alist, Galois::Field* gf) {
    uint64_t number_of_columns = alist.number_of_columns;
    uint64_t number_of_rows = alist.number_of_rows;
    std::vector<uint64_t> row_weights = alist.mwlist;
    uint64_t uli, ulj, ulk;
    uint64_t total_number_of_edges = 0;
    uint64_t edge_id = 0;
    std::vector<std::vector<uint64_t> > fullfill_table;
    std::vector<std::vector<uint64_t> >::iterator vit;
    std::vector<std::vector<std::vector<uint64_t> > > edgewise_fullfill_table;
    std::vector<std::vector<uint64_t> > u;
    std::vector<uint64_t> v;
    Galois::Element sum(gf, 0);
    Galois::Element zero(gf, 0);
    Galois::Element gfe(gf, 0);
    Galois::Element parity(gf, 0);
    uint64_t g;

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
        for(ulk = 0; ulk < factor_nodes[uli].edges.size(); ulk++) {
            u.clear();
            for(g = 0; g < (uint64_t)(gf->q); g++) {
                for(vit = fullfill_table.begin(); vit < fullfill_table.end(); ++vit) {
                    v = *vit;
                    if(v[ulk] == g) {
                        u.push_back(v);
                    }
                }
                edgewise_fullfill_table.push_back(u);
            }
            factor_nodes[uli].edges[ulk]->edgewise_fullfill_table = edgewise_fullfill_table;
            edgewise_fullfill_table.clear();
        }
    }
}