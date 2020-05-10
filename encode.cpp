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
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include <iomanip>
#include <filesystem>

#include <msgpack.hpp>

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

#include "generating_matrix.hpp"

extern std::vector<std::vector<Galois::Element> > gfq_matrix_product(std::vector<std::vector<Galois::Element> >&, std::vector<std::vector<Galois::Element> >&, const Galois::Field*);

std::vector<std::vector<Galois::Element> > load_generating_matrix(std::string file_name, const Galois::Field* gf) {
    std::ifstream ifs;
    Galois::Element gfe(gf, 0);
    generating_matrix G;
    uint64_t uli, ulj;
    uint64_t row_size, column_size;
    std::vector<Galois::Element> v;
    std::vector<std::vector<Galois::Element> > rtnv;
    uint64_t file_size;
    char* buf;

    try {
        file_size = std::filesystem::file_size(file_name);
    } catch(std::filesystem::filesystem_error e) {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }

    try {
        buf = new char(file_size);
    } catch(std::bad_alloc e) {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }

    ifs.open(file_name.c_str(), std::ios::in | std::ios::binary);
    if(!ifs) {
        std::cerr << "Can not open file: " << file_name << std::endl;
        exit(-1);
    }
    if(!ifs.read(buf, file_size)) {
        std::cerr << "Loading generating matrix file failed." << std::endl;
        exit(-1);
    }
    ifs.close();

    /* using msgpack, unpack generating matrix */
    try {
        msgpack::object_handle hd = msgpack::unpack(buf, file_size);
        const msgpack::object& obj = hd.get();
        obj.convert(G);
    } catch(std::bad_cast& e) {
        std::cerr << e.what() << std::endl;
    }

    if(G.characteristic != gf->q) {
        std::cerr << "characteristic do not match. " << G.characteristic << ", " << gf->q << std::endl;
        exit(-1);
    }

    row_size = G.contents.size();
    column_size = G.contents[0].size();
    for(ulj = 0; ulj < column_size; ulj++) v.push_back(gfe);
    for(uli = 0; uli < row_size; uli++) rtnv.push_back(v);
    for(uli = 0; uli < row_size; uli++) {
        for(ulj = 0; ulj < column_size; ulj++) {
            rtnv[uli][ulj].setValue(G.contents[uli][ulj]);
        }
    }

    delete[] buf;

    return(rtnv);
}

std::vector<Galois::Element> encode(std::vector<Galois::Element>& input, std::vector<std::vector<Galois::Element> >& G, const Galois::Field* gf) {
    uint64_t codeword_size = G.size();
    Galois::Element zero(gf, 0);
    std::vector<Galois::Element> v;
    uint64_t uli;
    std::vector<std::vector<Galois::Element> > B, C;

    for(uli = 0; uli < codeword_size; uli++) {
        v.push_back(zero);
    }
    B.push_back(input);
    C = gfq_matrix_product(G, B, gf);
    return(C[0]);
}