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
#include <string>
#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>
#include <msgpack.hpp>
#include "gfq_alist.hpp"
#include "generating_matrix.hpp"

extern std::vector<std::vector<Galois::Element> > make_generating_matrix(std::vector<std::vector<Galois::Element> >&, const Galois::Field*);

int main(int argc, char* argv[]) {
    std::ifstream ifs;
    std::ofstream ofs;
    std::vector<std::vector<Galois::Element> > G, H;
    uint64_t characteristic;
    generating_matrix g;
    const Galois::Field* gf;
    uint64_t number_of_rows, number_of_columns;
    uint64_t uli, ulj;
    msgpack::sbuffer buf;
    
    if(argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <GF(q) alist file name> <generating matrix file name>" << std::endl;
        return(EXIT_FAILURE);
    }
    gfq_alist alist = gfq_alist(argv[1]);
    characteristic = alist.characteristic;
    H = alist.make_dense();
    gf = H[0][0].field();
    G = make_generating_matrix(H, gf);
    number_of_rows = G.size();
    number_of_columns = G[0].size();
    std::vector<std::vector<uint64_t> > g_contents(number_of_rows, std::vector<uint64_t>(number_of_columns, 0));

    for(uli = 0; uli < number_of_rows; uli++) {
        for(ulj = 0; ulj < number_of_columns; ulj++) {
            g_contents[uli][ulj] = G[uli][ulj].value();
        }
    }

    g.characteristic = characteristic;
    g.contents = g_contents;
    msgpack::pack(buf, g);

    ofs.open(argv[2], std::ios::binary | std::ios::out);
    if(!ofs){
        std::cerr << "Can not open file: " << argv[2] << std::endl;
        return(EXIT_FAILURE);
    }
    ofs.write(buf.data(), buf.size());
    ofs.close();
    return(EXIT_SUCCESS);
}