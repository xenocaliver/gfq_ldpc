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

#ifndef GFQ_ALIST_H_
#define GFQ_ALIST_H_
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include <string>
#include <utility>
#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>
#include <boost/algorithm/string.hpp>  // for split function

#include "split.hpp"

class gfq_alist {
public:
    uint64_t number_of_columns;
    uint64_t number_of_rows;
    uint64_t characteristic;
    uint64_t max_column_weight;
    uint64_t max_row_weight;
    std::vector<uint64_t> nwlist;
    std::vector<uint64_t> mwlist;
    std::vector<std::vector<std::pair<uint64_t, uint64_t> > > nlist;
    std::vector<std::vector<std::pair<uint64_t, uint64_t> > > mlist;

    /* constructors */
    gfq_alist(void) : number_of_columns(0), number_of_rows(0), characteristic(0), max_column_weight(0), max_row_weight(0) {}
    gfq_alist(char *parity_check_file_name) {
        std::ifstream ifs;
        int rtnv = 0;

        ifs.open(parity_check_file_name);
        if(!ifs) {
            std::cerr << "Can not open file: " << parity_check_file_name << std::endl;
            exit(-1);
        }
        rtnv = gfq_parse_alist(ifs);
        if(rtnv != 1) {
            std::cerr << "Error in gfq_parse_alist" << std::endl;
            exit(-1);
        }
        ifs.close();
    }

    int gfq_parse_alist(std::ifstream& ifs) {
        std::string line;
        std::vector<std::string> v, w;
        std::vector<std::string>::iterator vit;
        uint64_t uli, ulk;
        std::vector<std::pair<uint64_t, uint64_t> > ev;
        std::string zero("0");

        getline(ifs, line);
        boost::trim(line);
        v = split(line, ' ');           /* split by empty string */
        if(v.size() != 3) {
            std::cerr << "Invalid ALIST FILE" << std::endl;
            return(-1);
        }
        try {
            this->number_of_columns = std::stoull(v[0], nullptr, 10);
            this->number_of_rows = std::stoull(v[1], nullptr, 10);
            this->characteristic = std::stoull(v[2], nullptr, 10);
        } catch(std::invalid_argument e) {
            std::cerr << e.what() << std::endl;
            return(-1);
        }
        std::pair<uint64_t, uint64_t> pos(0, 0);
        getline(ifs, line);
        boost::trim(line);
        v = split(line, ' ');
        if(v.size() != 2) {
            std::cerr << "Invalid ALIST FILE" << std::endl;
            return(-1);
        }
        try {
            this->max_column_weight = std::stoull(v[0], nullptr, 10);
            this->max_row_weight = std::stoull(v[1], nullptr, 10);
        } catch(std::invalid_argument e) {
            std::cerr << e.what() << std::endl;
            return(-1);
        }

        getline(ifs, line);
        boost::trim(line);
        v = split(line, ' ');
        for(vit = v.begin(); vit != v.end(); ++vit) {
            this->nwlist.push_back(std::stoull(*vit, nullptr, 10));
        }

        getline(ifs, line);
        boost::trim(line);
        v = split(line, ' ');
        for(vit = v.begin(); vit != v.end(); ++vit) {
            this->mwlist.push_back(std::stoull(*vit, nullptr, 10));
        }

        this->nlist.resize(this->number_of_columns);
        for(uli = 0; uli < this->number_of_columns; uli++) {
            getline(ifs, line);
            boost::trim(line);
            v = split(line, ' ');
            for(ulk = 0; ulk < v.size(); ulk += 2) {
                if(v[ulk + 1] == zero) continue;
                pos.first = std::stoull(v[ulk], nullptr, 10);
                pos.second = std::stoull(v[ulk + 1], nullptr, 10);
                this->nlist[uli].push_back(pos);
            }
        }

        this->mlist.resize(this->number_of_rows);
        for(uli = 0; uli < this->number_of_rows; uli++) {
            getline(ifs, line);
            boost::trim(line);
            v = split(line, ' ');
            for(ulk = 0; ulk < v.size(); ulk += 2) {
                if(v[ulk + 1] == zero) continue;
                pos.first = std::stoull(v[ulk], nullptr, 10);
                pos.second = std::stoull(v[ulk + 1], nullptr, 10);
                this->mlist[uli].push_back(pos);
            }
        }
        return(1);
    } 

    std::vector<std::vector<uint64_t> > make_dense(void) {
        uint64_t uli;
        std::vector<std::vector<uint64_t> > H(this->number_of_rows, std::vector<uint64_t>(this->number_of_columns, 0));
        std::vector<uint64_t> v;
        std::vector<std::pair<uint64_t, uint64_t> >::iterator vit;

        Galois::Field gf(this->characteristic);
        for(uli = 0; uli < this->number_of_rows; uli++) {
            for(vit = this->mlist[uli].begin(); vit != this->mlist[uli].end(); ++vit) {
                H[uli][vit->first - 1] = vit->second;
            }
        }
        return(H);
    }
};
#endif /* GFQ_ALIST_H_ */