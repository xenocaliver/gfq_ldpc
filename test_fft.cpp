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
#include <cmath>
#include <vector>
#include <map>
#include <iomanip>
#include <random>

#include <fftw3.h>

int main(int argc, char* argv[]) {
    uint64_t size = 10;
    uint64_t uli;
    std::vector<double> in(size, 0.0);
    std::vector<double> out(size, 0.0);
    std::vector<double> final(size, 0.0);
    std::mt19937_64 mt(1234);
    std::normal_distribution<double> gaussian(0, 0.1);

    for(uli = 0; uli < size; uli++) in[uli] = gaussian(mt);

    std::cout << "/**************** input *****************/" << std::endl;
    for(uli = 0; uli < size; uli++) {
        std::cout << std::fixed << std::setprecision(5) << in[uli] << " ";
    }
    std::cout << std::endl;
    fftw_plan plan = fftw_plan_r2r_1d(size, in.data(), out.data(), FFTW_REDFT11, FFTW_MEASURE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    std::cout << "/**************** output *****************/" << std::endl;
    for(uli = 0; uli < size; uli++) {
        std::cout << std::fixed << std::setprecision(5) << out[uli] << " ";
    }
    std::cout << std::endl;

    plan = fftw_plan_r2r_1d(size, out.data(), final.data(), FFTW_REDFT11, FFTW_MEASURE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    std::cout << "/**************** final *****************/" << std::endl;
    for(uli = 0; uli < size; uli++) {
        std::cout << std::fixed << std::setprecision(5) << final[uli]/(double)(2*size) << " ";
    }
    std::cout << std::endl;
    return(0);
}