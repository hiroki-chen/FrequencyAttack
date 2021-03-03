#include "util.hpp"

#include <random>
#include <algorithm>
#include <iostream>
// Debug

int main() {
    /*std::vector<std::string> c_space = {"sduahdsu","adsdsa","assssskk293","sass_pOa","sass223ssas","sadas-2033", "OOIS==-S"};
    std::vector<std::string> c = {"sduahdsu","adsdsa","adsdsa","assssskk293","assssskk293","assssskk293","OOIS==-S","sass_pOa","sass_pOa","sass_pOa","sass_pOa","sass223ssas","sass223ssas","sass223ssas","sass223ssas","sass223ssas","sadas-2033","sadas-2033","sadas-2033","sadas-2033","sadas-2033","sadas-2033"};
    std::vector<std::string> z_space = {"null","null","U","N","T", "K", "I"};
    std::vector<std::string> z = {"N","T", "N","K","K","T", "T","U","U","U","U", "I"};

    auto alpha = FH_binominal_attack<std::string>(c, z, c_space, z_space, 5, 0.95);

    for (auto item : alpha) {
        std::cout << item.first << " -> " << item.second.first << ", " << item.second.second << std::endl;
    }*/

    std::string key = "123456";
    output_information(get_accuracy_information(key, 100, 100000, 100.0, 50.0, 2, 16));
}
/*
 * l_p
adsdsa -> null
OOIS==-S -> null
adsdsa -> I
assssskk293 -> K
sadas-2033 -> U
sass223ssas -> T
sass_pOa -> N
sduahdsu -> null
 */

/*
 * frequency_ana
OOIS==-S -> null
adsdsa -> I
assssskk293 -> K
sadas-2033 -> U
sass223ssas -> T
sass_pOa -> N
sduahdsu -> null
 */

/*
 * non-crossing
OOIS==-S -> N
sadas-2033 -> null
sduahdsu -> null
 */

/*
 * FH_binominal_attack
I -> -50.74048, 1.83333
K -> 1.83333, 8.37024
N -> 8.37024, 10.3148
T -> 10.3148, 16.3068
U -> 16.3068, 33.481
 */