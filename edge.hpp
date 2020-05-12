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
#ifndef EDGE_HPP_
#define EDGE_HPP_
#include <cstdlib>
#include <cstdint>
#include <vector>

class edge {
public:
    uint64_t ID;
    std::vector<double> factor_to_variable_message;
    std::vector<double> variable_to_factor_message;
    std::vector<std::vector<std::vector<uint64_t> > > edgewise_fullfill_table;

    /* default constructor */
    edge(void) : ID(0) {}
    edge(uint64_t id) : ID(id) {}
};

#endif