#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <cstdio>
#include "util.h"
#include <omp.h>
#include <gmp.h>

#define N 1000000

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

bool is_prime_gmp(mpz_t n) {
    // Assume n is an odd number larger than 3
    //initialize sqrt_n and set it to the square root of n
    mpf_t sqrt_n, float_n;
    mpf_init(sqrt_n);
    mpf_init(float_n);
    mpf_set_z(float_n, n);
    mpf_sqrt(sqrt_n, float_n);
    
    //convert to an integer and round up
    mpz_t sqrt_n_z;
    mpz_init(sqrt_n_z);
    mpz_set_f(sqrt_n_z, sqrt_n);
    mpz_add_ui(sqrt_n_z, sqrt_n_z, 1);

    //Collect garbage
    mpf_clears(sqrt_n, float_n);

    std::cout << "n: ";
    mpz_out_str(NULL, 10, n);
    std::cout << "\nsqrt(n): ";
    mpz_out_str(NULL, 10, sqrt_n_z);
    std::cout << '\n';

    //Check divisibility
    bool output = true;
    mpz_t k;
    mpz_init_set_ui(k, 3);
    while ( mpz_cmp(k, sqrt_n_z) <= 0 ) {

        if ( mpz_divisible_p(n, k) ) {
            output = false;
            break;
        }

        mpz_add_ui(k, k, 2);
    }

    return output;
}

int main(int argc, char* argv[]) {

    double t0 = omp_get_wtime();
    int i = 0;
    int primelen = 0;
    uint64_t* primearray = bitmap_to_array( run_sequential_sieve(N), N, primelen);
    std::cout << "Found all primes below " << N << std::endl;
    std::cout << "Array length: " << primelen << std::endl;

    /*
    //Testing the prime sieve, if needed
    for (i = 0; i < 20; i++) {
        std::cout << primearray[i] << ", ";
    }
    std::cout << '\n';
    */

    //file output
    std::ofstream f_pn;
    f_pn.open("PerfectNumbers.txt");
    

    // Find the first few perfect numbers, up to the limit of 64-bit integers
    std::cout << "Mersenne prime \t Perfect number \n";

    for (i = 0; i <= 15; i++) {
        int64_t p = pow(2, primearray[i]) - 1;

        if (is_prime(p)) {
            std::cout << p << '\t' << uint64_t(p * pow(2, primearray[i] - 1));
            std::cout << std::endl;

            //f_pn << std::fixed << p << '\t' << int64_t(p * pow(2, primearray[i] - 1));
            //f_pn << std::endl;
        }
    }

    printarray(primearray, 30);

    int num_threads = 1;

    #pragma omp parallel 
    {
        num_threads = omp_get_num_threads();
    }

    std::cout << num_threads << " threads available.\n";

    std::list<std::string> perfnumlist;

    // check if each prime number p, results in a Mersenne prime
    for (i = 16; i < 650; i++) {

        uint64_t p = primearray[i];

        std::cout << "----------------------------\n";
        std::cout << "i = " << i << '\t' << "p = " << p << '\n';
        
        //Initialize prime_factor, and set it to 2**p - 1
        mpz_t prime_factor;
        mpz_init(prime_factor);
        mpz_ui_pow_ui(prime_factor, 2, p);
        mpz_sub_ui(prime_factor, prime_factor, 1);

        uint64_t pfui = mpz_get_ui(prime_factor);
        std::cout << "2^p - 1 = " << pfui << '\n';

        //Check if prime_factor is actually prime

        if ( mpz_probab_prime_p(prime_factor, 20) >= 1) {
            // Initialize and set power_of_2 to 2**(p - 1)
            mpz_t power_of_2;
            mpz_init(power_of_2);
            mpz_ui_pow_ui(power_of_2, 2, p - 1);

            //multiply to find the perfect number
            mpz_t perfnum;
            mpz_init(perfnum);
            mpz_mul(perfnum, prime_factor, power_of_2);

            mpz_out_str(NULL, 10, perfnum);
            std::cout << " is a perfect number.\n";
        }

        else {
            std::cout << "not a perfect number\n";
        }

        mpz_clear(prime_factor);
    }

    f_pn.close();
    
    std::cout << "\n";
    std::cout << "Time: " << omp_get_wtime() - t0 << std::endl;



    return 0;
}