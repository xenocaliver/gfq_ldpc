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
#include <cstdlib>
#include <vector>
#include <cstdint>
#include <iomanip>
#include <random>

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

extern double BPSK(double);

std::vector<std::vector<double> > transmit(std::vector<Galois::Element>& input, const Galois::Field* gf, double sigma, std::mt19937_64* mt) {
    uint64_t q = gf->q;
    uint64_t m = 0;                            /* multiplicative order */
    uint64_t uli, ulj;
    uint64_t value;
    uint64_t b;

    std::normal_distribution<double> gaussian(0, sigma);
    /* get multiplicative order */
    while(q != 1) {
        q = q >> 1;
        m++;
    }    
    std::vector<std::vector<double> > rtnv(input.size(), std::vector<double>(m, 0));
    for(uli = 0; uli < input.size(); uli++) {
        value = input[uli].value();
        for(ulj = 0; ulj < m; ulj++) {
            b = value & 1;
            rtnv[uli][ulj] = BPSK((double)b) + gaussian(*mt);
            value = value >> 1;
        }
    }
    return(rtnv);
}