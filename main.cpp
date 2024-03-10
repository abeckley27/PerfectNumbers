#include <iostream>
#include <cmath>
#include "util.h"

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

    for (int i = 1; i < 10; i++) {
        std::cout << i << '\t' << is_prime(i) << std::endl;
    }


    return 0;
}