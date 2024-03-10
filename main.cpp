#include <iostream>
#include <fstream>
#include <omp.h>
#include "util.h"

#define N 10000000

bool is_prime(int64_t n) {
    bool output = ((n & 1) || (n == 2)) && (n > 1);
    
    if ((n & 1) && (n > 3)) {
        for (int64_t k = 3; k <= int64_t(ceil(sqrt(n))); k += 2) {
            if (n % k == 0) {
                output = false;
                break;
            }
        }
    }

    return output;
}

int main(int argc, char* argv[]) {

    double t0 = omp_get_wtime();
    int primelen = 0;
    int64_t* primearray = bitmap_to_array( run_sequential_sieve(N), N, primelen);
    std::cout << "Found all primes below " << N << std::endl;
    std::cout << "Array length: " << primelen << std::endl;

    /*
    //Testing the prime sieve, if needed
    for (int i = 0; i < 20; i++) {
        std::cout << primearray[i] << ", ";
    }
    std::cout << '\n';
    */

    std::ofstream f;
    f.open("PerfectNumbers.txt");
    
    std::cout << "Mersenne prime \t Perfect number \n";

    for (int i = 0; i < 15; i++) {
        int64_t p = pow(2, primearray[i]) - 1;
        if (is_prime(p)) {
            std::cout << p << '\t' << uint64_t(p * pow(2, primearray[i] - 1));
            std::cout << std::endl;

            //f << std::fixed << p << '\t' << int64_t(p * pow(2, primearray[i] - 1));
            //f << std::endl;
        }
    }

    f.close();
    
    std::cout << "\n";
    std::cout << "Time: " << omp_get_wtime() - t0 << std::endl;

    


    return 0;
}