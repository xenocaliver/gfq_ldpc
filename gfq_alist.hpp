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
#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>
#include <boost/algorithm/string.hpp>  // for split function

using nonzero_element = std::pair<uint64_t, Galois::Element>;

class gfq_alist {
public:
    uint64_t number_of_columns;
    uint64_t number_of_rows;
    uint64_t characteristic;
    uint64_t max_column_weight;
    uint64_t max_row_weight;
    std::vector<uint64_t> nwlist;
    std::vector<uint64_t> mwlist;
    std::vector<std::vector<nonzero_element> > nlist;
    std::vector<std::vector<nonzero_element> > mlist

    /* constructors */
    gfq_alist(void) : number_of_columns(0), number_of_rows(0), characteristic(0), max_column_weight(0), max_row_weight(0) {}
    gfq_alist(char *parity_check_file_name) {
        std::ifstream ifs;
        int rtnv = 0;

        ifs.open(parity_check_file_name);
        if(!ifs) {
            std::cerr << "Can not open file: " << parity_check_matrix_file_name << std::endl;
            exit(-1);
        }
        rtnv = gfq_parse_alist(ifs);
        if(rtnv != 1) {
            std::cerr << "Error in gfq_parse_alist" << std::endl;
            exit(-1);
        }
    }

    int gfq_parse_alist(std::ifstream& ifs) {
        std::string line;
        std::vector<std::string> v;
        std::vector<std::string>::iterator vit;
        uint64_t uli, ulj;
        nonzero_element e;
        std::vector<nonzero_element> ev;
        Galois::Element gfe;
        Galois::Field gf;

        getline(ifs, line);
        boost::algorithm::split(v, line, boost::is_space());           /* split by empty string */
        if(v.size() != 3) {
            std::cerr << "Invalid ALIST FILE" << std::endl;
            return(-1);
        }
        this->number_of_columns = strtoul(v[0], nullptr, 10);
        this->number_of_rows = strtoul(v[1], nullptr, 10);
        this->characteristic = strtoul(v[2], nullptr, 10);
        gf.resize(this->characteristic);
        gfe(&gf, 0);
        getline(ifs, line);
        boost::algorithm::split(v, line, boost::is_space());
        if(v.size() != 2) {
            std::cerr << "Invalid ALIST FILE" << std::endl;
            return(-1);
        }
        this->max_column_weight = strtoul(v[0], nullptr, 10);
        this->max_row_weight = strtoul(v[1], nullptr, 10);

        getline(ifs, line);
        boost::algorithm::split(v, line, boost::is_space());
        for(vit = v.begin(), vit != v.end(), ++vit) {
            this->nwlist.push_back(strtoul(*vit, nullptr, 10));
        }

        getline(ifs, line);
        boost::algorithm::split(v, line, boost::is_space());
        for(vit = v.begin(), vit != v.end(), ++vit) {
            this->mwlist.push_back(strtoul(*vit, nullptr, 10));
        }

        this->nlist.resize(this->number_of_columns);
        for(uli = 0; uli < this->number_of_columns; uli++) {
            getline(ifs, line);
            boost::algorithm::split(v, line, is_any_of("\t "), boost::token_compress_on);
            for(ulj = 0; ulj < v.size(); ulj += 2){
                if(v[ulj + 1] == '0') continue; 
                e.first = strtoul(v[ulj], nullptr, 10);
                gfe.setValue(strtorul(v[ulj + 1], nullptr, 10));
                e.second = gfe
                this->nlist[uli].push_back(e);
            }
        }

        this->mlist.resize(this->number_of_rows);
       for(uli = 0; uli < this->number_of_rows; uli++) {
            getline(ifs, line);
            boost::algorithm::split(v, line, is_any_of("\t "), boost::token_compress_on);
            for(ulj = 0; ulj < v.size(); ulj += 2){
                if(v[ulj + 1] == '0') continue; 
                e.first = strtoul(v[ulj], nullptr, 10);
                gfe.setValue(strtorul(v[ulj + 1], nullptr, 10));
                e.second = gfe
                this->mlist[uli].push_back(e);
            }
        }
        return(1);
    }
};
#endif /* GFQ_ALIST_H_ */