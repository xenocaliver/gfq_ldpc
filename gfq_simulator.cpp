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
#include <random>
#include <thread>
#include <chrono>
#include <format>

#include <galois++/array2d.h>
#include <galois++/element.h>
#include <galois++/field.h>
#include <galois++/fwd.h>
#include <galois++/primes.h>

#include "edge.hpp"
#include "gfq_alist.hpp"
#include "variable_node.hpp"
#include "factor_node.hpp"
#include "generating_matrix.hpp"

extern std::vector<std::vector<double> > get_a_priori_probability(std::vector<std::vector<double> >, double, const Galois::Field*);
extern void construct_factor_graph(std::vector<variable_node>&, std::vector<factor_node>&, std::vector<edge>&, gfq_alist, Galois::Field*);
extern std::vector<std::vector<double> > transmit(std::vector<Galois::Element>&, const Galois::Field*, double, std::mt19937_64*);
extern std::vector<std::vector<Galois::Element> > load_generating_matrix(std::string, const Galois::Field*);
extern bool parity_check(std::vector<uint64_t>, std::vector<std::vector<std::pair<uint64_t, uint64_t> > >, const Galois::Field*); 
extern std::vector<Galois::Element> encode(std::vector<Galois::Element>&, std::vector<std::vector<Galois::Element> >&, const Galois::Field*);
extern uint64_t rdtsc(void);

/* for multi-thread processing */
std::vector<std::vector<uint64_t> > get_modulo_ids(uint64_t number_of_nodes, uint64_t concurrency) {
    std::vector<uint64_t> v;
    std::vector<std::vector<uint64_t> > rtnv;
    uint64_t uli, ulj;
    uint64_t q, r;

    q = number_of_nodes/concurrency;
    r = number_of_nodes%concurrency;

    for(uli = 0; uli < concurrency; uli++) {
        for(ulj = 0; ulj < q; ulj++) {
            v.push_back(concurrency*uli + ulj);
        }
        rtnv.push_back(v);
        v.clear();
    }
    for(uli = q*concurrency; uli < q*concurrency + r; uli++) v.push_back(uli);
    rtnv.push_back(v);
    return(rtnv);
}

void update_factor_nodes(std::vector<factor_node>& factor_nodes, std::vector<uint64_t>& node_ids, Galois::Field* gf) {
    std::vector<uint64_t>::iterator vit;

    for(vit = node_ids.begin(); vit != node_ids.end(); ++vit) {
        factor_nodes[*vit].update_messages(gf);
    }
}

void update_variable_nodes(std::vector<variable_node>& variable_nodes, std::vector<uint64_t>& node_ids, Galois::Field* gf) {
    std::vector<uint64_t>::iterator vit;

    for(vit = node_ids.begin(); vit != node_ids.end(); ++vit) {
        variable_nodes[*vit].update_messages(gf);
    }
}

int main(int argc, char* argv[]) {
    double sigma;
    std::vector<uint64_t> input;
    uint64_t dimension;
    uint64_t trial = 0;
    std::vector<std::vector<Galois::Element> > generating_matrix;
    std::vector<factor_node> factor_nodes;
    std::vector<variable_node> variable_nodes;
    std::vector<edge> edges;
    std::vector<std::vector<double> > received_signals;
    std::vector<std::vector<double> > a_priori_probability;
    uint64_t uli, ulj;
    uint64_t g;
    std::vector<Galois::Element> input_vector;
    std::vector<Galois::Element> codeword;
    std::vector<Galois::Element> decoded_word;
    uint64_t iteration_limit;
    uint64_t iteration;
    bool parity_check_result;
    uint64_t error_count = 0;
    double frame_error_rate;
    uint64_t hardware_concurrency;
    std::vector<std::vector<uint64_t> > factor_modulo, variable_modulo;
    std::vector<std::thread> factor_threads, variable_threads;
    std::vector<std::thread>::iterator tit;
    uint64_t elapsed_count_in_microseconds;
    uint64_t sum_elapsed_time = 0;
    double average_execution_time = 0.0;
    uint64_t progress_interval = 100;

    if(argc != 7) {
        std::cerr << "Usage: " << argv[0] << " <alist file> <generating matrix file> <number of trial> <sigma> <iteration limit> <print progress interval>" << std::endl;
        return(EXIT_FAILURE);
    }

    /* preparation */
    gfq_alist alist(argv[1]);
    dimension = alist.number_of_columns - alist.number_of_rows;
    std::vector<uint64_t> speculated_codeword(alist.number_of_columns, 0);
    iteration_limit = std::strtoull(argv[5], nullptr, 10);
    Galois::Field gf(alist.characteristic);
    Galois::Element gfe(&gf, 0);
    sigma = std::stod(argv[4]);
    std::mt19937_64 mt(1234);
    progress_interval = std::strtoull(argv[6], nullptr, 10);
    std::uniform_int_distribution<uint64_t> input_symbol(0, alist.characteristic - 1);
    uint64_t number_of_trial = std::strtoull(argv[3], nullptr, 10);
    generating_matrix = load_generating_matrix(std::string(argv[2]), &gf);
    construct_factor_graph(variable_nodes, factor_nodes, edges, alist, &gf);
    for(uli = 0; uli < dimension; uli++) input_vector.push_back(gfe);
    for(uli = 0; uli < variable_nodes.size(); uli++) {
        for(ulj = 0; ulj < variable_nodes[uli].edges.size(); ulj++) {
            variable_nodes[uli].edges[ulj]->variable_to_factor_message.resize((uint64_t)(gf.q));
            variable_nodes[uli].edges[ulj]->factor_to_variable_message.resize((uint64_t)(gf.q));
        }
    }

    /* for multi-threading */
    hardware_concurrency = std::thread::hardware_concurrency();                         /* get system's max number of threads */
    if(hardware_concurrency == 0) {
        std::cout << "Hardware concurrency = 0. You must modify this code for single thread." << std::endl;
        return(-1);
    }

    std::cout << "hardware concurrency = " << hardware_concurrency << std::endl;
    factor_modulo = get_modulo_ids(alist.number_of_rows, hardware_concurrency);
    variable_modulo = get_modulo_ids(alist.number_of_columns, hardware_concurrency);

    try {
        std::cout << std::format("{:^12s}{:s}","count", "  average execution time[us]") << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
    } catch(const std::format_error& e) {
        std::cerr << e.what() << std::endl;
        return(-1);
    }
    /* main loop */
    for(trial = 0; trial < number_of_trial; trial++) {
        auto start = std::chrono::high_resolution_clock::now();
        /* generate input words */
        for(uli = 0; uli < dimension; uli++) {
            g = input_symbol(mt);
            gfe.setValue(g);
            input_vector[uli] = gfe;
        }
        /* encode */
        codeword = encode(input_vector, generating_matrix, &gf);
#ifdef DEBUG
        std::cout << "*** encoded word ***" << std::endl;
        for(uli = 0; uli < codeword.size(); uli++) {
            std::cout << std::setw(2) << codeword[uli].value() << " ";
        }
        std::cout << std::endl;
#endif
        /* add awgn noise */
        received_signals = transmit(codeword, &gf, sigma, &mt);
        /* calculate a priori probability */
        a_priori_probability = get_a_priori_probability(received_signals, sigma, &gf);
        for(uli = 0; uli < variable_nodes.size(); uli++) variable_nodes[uli].channel_output_probability = a_priori_probability[uli];
        /* initialize messages */
        for(uli = 0; uli < variable_nodes.size(); uli++) {
            for(ulj = 0; ulj < variable_nodes[uli].edges.size(); ulj++) {
                for(g = 0; g < (uint64_t)(gf.q); g++) {
                    variable_nodes[uli].edges[ulj]->variable_to_factor_message[g] = a_priori_probability[uli][g];
                }
            }
        }
        for(uli = 0; uli < factor_nodes.size(); uli++) {
            for(ulj = 0; ulj < factor_nodes[uli].edges.size(); ulj++) {
                for(g = 0; g < (uint64_t)(gf.q); g++) {
                    factor_nodes[uli].edges[ulj]->factor_to_variable_message[g] = 1.0;
                }
            }
        }
        parity_check_result = false;
        /* do decoding process */
        for(iteration = 0; iteration < iteration_limit; iteration++) {
            /* update messages */
            for(uli = 0; uli < hardware_concurrency + 1; uli++) {
                std::thread worker(update_factor_nodes, std::ref(factor_nodes), std::ref(factor_modulo[uli]), &gf);
                factor_threads.emplace_back(std::move(worker));
            }
            for(tit = factor_threads.begin(); tit != factor_threads.end(); ++tit) tit->join();
            factor_threads.clear();

            for(uli = 0; uli < hardware_concurrency + 1; uli++) {
                std::thread worker(update_variable_nodes, std::ref(variable_nodes), std::ref(variable_modulo[uli]), &gf);
                variable_threads.emplace_back(std::move(worker));
            }
            for(tit = variable_threads.begin(); tit != variable_threads.end(); ++tit) tit->join();
            variable_threads.clear();

            /* speculate code word */
            for(uli = 0; uli < variable_nodes.size(); uli++) speculated_codeword[uli] = variable_nodes[uli].speculate_temporal_symbol(&gf);
            parity_check_result = parity_check(speculated_codeword, alist.mlist, &gf);
            if(parity_check_result == true) break;
        }
        if(parity_check_result == false) {
            error_count++;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start;
        std::chrono::microseconds d = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        elapsed_count_in_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(d).count();
        sum_elapsed_time += elapsed_count_in_microseconds;
        if((trial + 1)%progress_interval == 0) {
            average_execution_time = (double)sum_elapsed_time/(double)(trial + 1);
            try{
                std::cout << std::format("{:>12d}{:>28f}", trial + 1, average_execution_time) << std::endl;
            } catch(const std::format_error& e) {
                std::cerr << e.what() << std::endl;
                return(-1);
            }
        }
    }
    frame_error_rate = (double)error_count/(double)number_of_trial;
    std::cout << "*** simulation result ***" << std::endl;
    std::cout << frame_error_rate << std::endl;
    return(EXIT_SUCCESS);
}