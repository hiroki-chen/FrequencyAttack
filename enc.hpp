#ifndef _ENC_HPP_
#define _ENC_HPP_

#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <random>

char chars[] = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};

std::string gen_random(int length) {
    static auto& chrs = "0123456789"
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    thread_local static std::mt19937 rg{std::random_device{}()};
    thread_local static std::uniform_int_distribution<std::string::size_type> pick(0, sizeof(chrs) - 2);

    std::string s;

    s.reserve(length);

    while(length--)
        s += chrs[pick(rg)];

    return s;
}

std::vector<std::string> get_random(const long& size) {
    std::vector<std::string> ret(size);
    for (long i = 0; i < size; i++) {
        ret[i] = gen_random(12);
    }
    return ret;
}


std::map<std::string, std::string> getEncryption(const std::string& key,
                                                 const std::vector<std::string>& plaintext) {
    /* mapping between plaintext and its ciphertext. */
    std::map<std::string, std::string> ret;
    int s = key.size();
    for (int i = 0; i < plaintext.size(); i++) {
        std::string tmp = "";
        int sum = 0;
        for (int j = 0; j < plaintext[i].size(); j++) {
            sum += int(plaintext[i][j]);
        }

        for (int j = 0; j < plaintext[i].size(); j++) {
            tmp.push_back(chars[(long(key[j % s]) * long(plaintext[i][j]) + sum % long(plaintext[i][j])) % 26]);
        }

        ret[plaintext[i]] = tmp;
    }

    return ret;
}

#endif