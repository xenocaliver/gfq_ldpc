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
#include <cmath>

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

double BPSK(double x) {
    return(1.0 - 2.0*x);
}

std::vector<std::vector<double> > get_a_priori_probability(std::vector<std::vector<double> > channel_output, double sigma, const Galois::Field* gf) {
    uint64_t m = 0;
    uint64_t uli, ulj, ulk;
    uint64_t value;
    int64_t b;
    double prior_probability;
    
    uint64_t q = (uint64_t)(gf->q);
    /* get multiplicative degree */
    m = channel_output[0].size();
    std::vector<std::vector<double> > f(channel_output.size(), std::vector<double>(q, 0));

    for(uli = 0; uli < channel_output.size(); uli++) {
        for(ulj = 0; ulj < q; ulj++) {
            value = ulj;
            prior_probability = 1.0;
            for(ulk = 0; ulk < m; ulk++) {
                b = (int64_t)(value & 1);
                prior_probability *= (1 + exp(2.0*BPSK((double)b)*channel_output[uli][ulk]/(sigma*sigma)));
                value = value >> 1;
            }
            f[uli][ulj] = 1.0/prior_probability;
        }
    }
    return(f);
}

bool parity_check(std::vector<uint64_t> speculate_temporal_symbols, std::vector<std::vector<std::pair<uint64_t, uint64_t> > > mlist, const Galois::Field* gf) {
    uint64_t uli, ulj;
    Galois::Element temporal_symbol(gf, 0);
    Galois::Element parity_check_element(gf, 0);
    Galois::Element zero(gf, 0);
    Galois::Element sum(gf, 0);

    std::cout << "speculate code word" << std::endl;
    for(uli = 0; uli < speculate_temporal_symbols.size(); uli++) {
        std::cout << std::setw(2) << speculate_temporal_symbols[uli] << " ";
    }
    std::cout << std::endl;

    for(uli = 0; uli < mlist.size(); uli++) {
        sum = zero;
        for(ulj = 0; ulj < mlist[uli].size(); ulj++) {
            parity_check_element.setValue(mlist[uli][ulj].second);
            temporal_symbol.setValue(speculate_temporal_symbols[mlist[uli][ulj].first]);
            sum += parity_check_element*temporal_symbol;
            std::cout << std::setw(2) << temporal_symbol.value() << "*" << parity_check_element.value() << " ";
        }
        std::cout << sum.value() << std::endl;
        if(sum != zero) {
            return(false);
        }
    }
    return(true);
}