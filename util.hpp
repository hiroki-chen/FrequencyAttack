 #include "AttackAlgorithms.hpp"
 #include "frequency_attack.hpp"
 #include "enc.hpp"

 #include <iostream>

template <typename DataType>
double KL_divergence(std::map<int, int> p,
                     std::map<int, int> q,
                     const std::vector<DataType>& p_data,
                     const std::vector<DataType>& q_data,
                     const int& p_size,
                     const int& q_size
                     ) {
    double ret = 0;

    for (auto item : p) {
        if (find(q_data.begin(), q_data.end(), p_data[item.first]) != q_data.end()) {
            double px = item.second * 1.0 / p_size;
            double qx = q[item.first] * 1.0 / q_size;
            ret += px * log(px  / qx);
        }
    }

    return ret;
}

std::vector<std::pair<double, double>> get_accuracy_information(const std::string& key,
                                                                long long z_space_size = 50,
                                                                long long z_size = 100000,
                                                                double mean = 20.0,
                                                                double variance = 10.0,
                                                                int interval = 10,
                                                                int p = 2,
                                                                bool distort = true) {
    /* z_space Obtained by adversary. */
    std::vector<std::string> z_space = get_random(z_space_size);

    /* Plaintext to ciphertext map */
    std::map<std::string, std::string> enc = getEncryption(key, z_space);
    std::sort(z_space.begin(), z_space.end());

    /* Generate the distribution of z. */
    std::vector<std::string> z(z_size);
    int i = 0;
    int index = 0;
    std::default_random_engine e;
    std::normal_distribution<double> normal(mean, variance);
    while(i < z_size) {
        int time = abs(std::lround(normal(e)));
        for (int j = i; j < z.size() && j < i + time; j++) {
            z[j] = z_space[index % z_space_size];
        }
        index ++;
        i += time;
    }
    auto hist_z = hist<std::string>(z, z_space);

    /* Distort the z to real distribution by either adding to it, or deleting from the it. */
    std::vector<std::string> real_z;
    for (int i = 0; i < z.size(); i++) {
        std::string s = z[i];
        real_z.push_back(s);
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(real_z.begin(), real_z.end(), g);

    std::normal_distribution<double> dist2(interval, interval * 1.0 / 10);
    
    if (true == distort) {
        int time = abs(std::lround(dist2(e)));
        int i = 0;
        for (std::string& s : real_z) {
            if (i++ % time == 0) {
                /* How to add ? How many times will an item be repeated? */
                real_z.push_back(s);
                /*auto pos_in_z_space = find(z_space.begin(), z_space.end(), s) - z_space.begin() + 1;
                int count = hist_z[pos_in_z_space];

                for (int j = 0; j < std::lround(abs(z_space_size - count)); j++) {
                    real_z.push_back(s);
                }*/
                time = abs(std::lround(dist2(e)));
            }
        }
    }

    /* Generate c and c_space. */
    auto map = getEncryption(key, real_z);
    std::vector<std::string> c_space;
    std::vector<std::string> c;
    for (auto item : map) {
        c_space.push_back(item.second);
    }
    for (int i = 0; i < real_z.size(); i++) {
        c.push_back(map[real_z[i]]);
    }
    std::vector<std::pair<double, double>> ret;

    auto hist_real_z = hist<std::string>(real_z, z_space);

    auto alpha = l_p_optimization(c, z, c_space, z_space, p);
    auto beta = frequency_ana(c, z, c_space, z_space);
    int correct = 0;
    
    for (auto item : alpha) {
        if (enc[item.second] == item.first) {
            auto pos_in_z_space = find(z_space.begin(), z_space.end(), item.second) - z_space.begin() + 1;
            correct += hist_real_z[pos_in_z_space];
        }
    }
    double accuracy = correct * 1.0 / real_z.size();
    double kld = KL_divergence<std::string>(hist_real_z, hist_z, z_space, z_space, real_z.size(), z.size());
    ret.push_back(std::make_pair(accuracy, kld));

    correct = 0;
    for (auto item : beta) {
        if (enc[item.second] == item.first) {
            auto pos_in_z_space = find(z_space.begin(), z_space.end(), item.second) - z_space.begin() + 1;
            correct += hist_real_z[pos_in_z_space];
        }
    }
    accuracy = correct * 1.0 / real_z.size();
    ret.push_back(std::make_pair(accuracy, kld));

    return ret;
}

void output_information(const std::vector<std::pair<double, double>>& info) {
    std::cout << "KL_divergence: " << info[0].second << std::endl;
    std::cout << "l_p_optimization:" << std::endl;
    std::cout << "accuracy: " << info[0].first << std::endl;
    std::cout << std::endl;
    std::cout << "frequency_analysis:" << std::endl;
    std::cout << "accuracy: " << info[1].first << std::endl;
}