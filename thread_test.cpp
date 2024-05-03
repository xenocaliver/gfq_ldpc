#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include <iomanip>
#include <cmath>
#include <random>
#include <thread>

void go_thread(uint64_t a, uint64_t b) {
    uint64_t c;

    c = a + b;
    std::cout << "result = " << c << std::endl;
}

int main(int argc, char* argv[]) {
    uint64_t a = 10;
    uint64_t b = 21;
    uint64_t hardware_concurrency;
    uint64_t uli;
    std::vector<std::thread> threads;
    std::vector<std::thread>::iterator tit;

    hardware_concurrency = std::thread::hardware_concurrency();                         /* get system's max number of threads */
    if(hardware_concurrency == 0) {
        std::cout << "Hardware concurrency = 0. You must modify this code for single thread." << std::endl;
        return(-1);
    }
    std::cout << "hardware_concurrency = " << hardware_concurrency << std::endl;

    for(uli = 0; uli < hardware_concurrency; uli++) {
        std::thread worker(go_thread, std::ref(a), std::ref(b));
        threads.emplace_back(std::move(worker));
    }
    for(tit = threads.begin(); tit != threads.end(); ++tit) tit->join();
    return(0);
}