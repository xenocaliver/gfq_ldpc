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

class gfq_alist {
public:
    uint64_t number_of_columns;
    uint64_t number_of_rows;
    uint64_t characteristic;
    uint64_t max_column_weight;
    uint64_t max_row_weight;
    std::vector<uint64_t> nwlist;
    std::vector<uint64_t> mwlist;
    std::vector<std::vector<std::pair<uint64_t, Galois::Element> > > nlist;
    std::vector<std::vector<std::pair<uint64_t, Galois::Element> > > mlist;
    std::string zero_string = "0";

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
        std::vector<std::string> v;
        std::vector<std::string>::iterator vit;
        uint64_t uli, ulj;
        std::vector<std::pair<uint64_t, Galois::Element> > ev;

        getline(ifs, line);
        boost::algorithm::split(v, line, boost::is_space());           /* split by empty string */
        if(v.size() != 3) {
            std::cerr << "Invalid ALIST FILE" << std::endl;
            return(-1);
        }
        this->number_of_columns = std::stoi(v[0]);
        this->number_of_rows = std::stoi(v[1]);
        this->characteristic = std::stoi(v[2]);
        Galois::Field gf(this->characteristic);
        Galois::Element gfe = Galois::Element(&gf, 0);
        std::pair<uint64_t, Galois::Element> e(0, gfe);
        getline(ifs, line);
        boost::algorithm::split(v, line, boost::is_space());
        if(v.size() != 2) {
            std::cerr << "Invalid ALIST FILE" << std::endl;
            return(-1);
        }
        this->max_column_weight = std::stoi(v[0]);
        this->max_row_weight = std::stoi(v[1]);

        getline(ifs, line);
        boost::algorithm::split(v, line, boost::is_space());
        for(vit = v.begin(); vit != v.end(); ++vit) {
            this->nwlist.push_back(std::stoi(*vit));
        }

        getline(ifs, line);
        boost::algorithm::split(v, line, boost::is_space());
        for(vit = v.begin(); vit != v.end(); ++vit) {
            this->mwlist.push_back(std::stoi(*vit));
        }

        this->nlist.resize(this->number_of_columns);
        for(uli = 0; uli < this->number_of_columns; uli++) {
            getline(ifs, line);
            boost::algorithm::split(v, line, boost::is_any_of("\t "), boost::token_compress_on);
            for(ulj = 0; ulj < v.size(); ulj += 2){
                if(v[ulj + 1] == zero_string) continue; 
                e.first = std::stoi(v[ulj]);
                gfe.setValue(std::stoi(v[ulj + 1]));
                e.second = gfe;
                this->nlist[uli].push_back(e);
            }
        }

        this->mlist.resize(this->number_of_rows);
       for(uli = 0; uli < this->number_of_rows; uli++) {
            getline(ifs, line);
            boost::algorithm::split(v, line, boost::is_any_of("\t "), boost::token_compress_on);
            for(ulj = 0; ulj < v.size(); ulj += 2){
                if(v[ulj + 1] == zero_string) continue; 
                e.first = std::stoi(v[ulj]);
                gfe.setValue(std::stoi(v[ulj + 1]));
                e.second = gfe;
                this->mlist[uli].push_back(e);
            }
        }
        return(1);
    } 

    std::vector<std::vector<Galois::Element> > make_dense(void) {
        Galois::Field gf(this->characteristic);
        Galois::Element zero(&gf, 0);
        uint64_t uli, ulj;
        std::vector<std::vector<Galois::Element> > H;
        std::vector<Galois::Element> v;
        std::vector<std::pair<uint64_t, Galois::Element> >::iterator vit;

        for(ulj = 0; ulj < this->number_of_columns; ulj++) {
            v.push_back(zero);
        }
        for(uli = 0; uli < this->number_of_rows; uli++) {
            H.push_back(v);
        }

        for(uli = 0; uli < this->number_of_rows; uli++) {
            for(vit = this->mlist[uli].begin(); vit != this->mlist[uli].end(); ++vit) {
                H[uli][vit->first - 1] = vit->second;
            }
        }
        return(H);
    }
};
#endif /* GFQ_ALIST_H_ */