#ifndef utility_fns_h
#define utility_fns_h

#include <iostream>
#include <cstdint>
#include <cmath>

// Prints an array to console, for debugging
void printarray(int64_t* a, int64_t len) {
	int_fast64_t i;
	for (i = 0; i < len; i++) {
		std::cout << a[i] << '\t';
	}
	std::cout << std::endl;
}

// This function allocates memory
bool* run_sequential_sieve(const int64_t N) {
	//bitmap for primality
	bool* s = new bool[N+1];

	// 0 and 1 are not prime
	s[0] = 0;
	s[1] = 0;
	// assume everything else is
	int64_t i;
	for (i = 2; i <= N; i++) { s[i] = 1; }


	//Remove multiples of primes
	for (int64_t x = 2; x <= ceil(sqrt(N)); x++) {
		if (s[x]) {
			for (int64_t multiple = x * x; multiple < N + 1; multiple += x) {
				s[multiple] = 0;
			}
		}
	}

	return s;
}

int64_t primecount(bool* s, int64_t len) {
    int64_t output = 0;
	for (int64_t i = 2; i < len; i++) { output += s[i]; }
    return output;
}

// This function allocates memory
int64_t* bitmap_to_array(bool* s, int64_t len, int& olen) {
	int64_t pcount = primecount(s, len);
    int64_t* output = new int64_t[pcount];
	int64_t i, j;
	j = 0;

	for (i = 2; i < len; i++) {
		if (s[i]) {
			output[j] = i;
			j++;
		}
	}

	olen = pcount;

	return output;
}

// This function checks if n is prime by first generating a list of prime numbers
bool standalone_is_prime(int64_t n) {

    // Generate an array of all prime numbers up to sqrt(n)
    const int32_t root_n = int32_t(ceil(sqrt(n)));
    bool* s = run_sequential_sieve(root_n);
    const int64_t L1 = primecount(s, root_n + 1);
	int x;
    int64_t* prime_array = bitmap_to_array(s, root_n + 1, x);
    
    // check for prime divisors
    bool output = (n > 1) && ( (n & 1) || (n == 2) );
    int32_t i = 1;
    while (output && (i < L1)) {
        output = (n % prime_array[i]);
        i++;
    }

    delete[] s;
	delete[] prime_array;

    return output;
}

int64_t pi(int64_t n) {
    bool* s = run_sequential_sieve(n);
    int64_t x = primecount(s, n);
    delete[] s;
    return x;
}


#endif
