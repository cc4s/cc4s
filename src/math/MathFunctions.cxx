#include <math/MathFunctions.hpp>
#include <cmath>

// TODO: Create a unit test specialized logger
int cc4s::permutationSign(std::string original, std::string permuted) {
    // Number of permutations
    int permutations = 0;
    int length = original.length();
    // Helping variables
    char prev;
    char next;
    std::string cycle = "";
    std::string visited = "";
    bool end_of_cycle = false;
    // At most we will find length cycles
    for (int c = 0; c < length; ++c) {
        cycle = "";
        end_of_cycle = false;
        prev = original[c];
        next = permuted[c];
        std::cout << "Prev " << prev << " next " << next<< std::endl;
        if (visited.find(prev) != std::string::npos)
            // If prev is already visited, skip it
            continue;
        if (prev == next) {
            // Do not bother with fix points
            visited += prev;
            continue;
        }
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j) {
                if (original[j] == next) {
                    cycle += prev;
                    next = permuted[j];
                    prev = original[j];
                    if (next == cycle[0]) {
                        cycle += prev;
                        end_of_cycle = true;
                        continue;
                    }
                }
            }
            if (end_of_cycle) {
                std::cout << "\tCycle = " << cycle << std::endl;
                break;
            }
        }
        permutations += cycle.length() - 1;
        visited += cycle;
        std::cout << "visited  " << visited << std::endl;
    }
    std::cout << "P ( " << original << " -> " << permuted << " ) "
              << "  =>  (-1)^" << permutations <<  std::endl;
    std::cout <<  pow(-1, permutations) <<  std::endl;
    return pow(-1, permutations);
}
