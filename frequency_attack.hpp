/*
 * DataType: The type of the data. You can denote it as string as well.
 * !! You need to ensure that m_space = c_space by using padding !
 */

#include <vector>
#include <cfloat>
#include <limits>
#include <algorithm>
#include <utility>
#include <set>
#include <map>
#include <cmath>
#include <iostream>

#include "AttackAlgorithms.hpp"

/*
 * Input: 
 *      c: ciphertext;
 *      z: auxiliary dataset;
 *      c_space: the domain of c；
 *      m_space: the domain of z, i.e., the message. (Because z is over M_k)
 * Output: 
 *      alpha: mapping relationship between ciphertext space C and message space M.
 *             That is, outputs the corresponding plain text of each ciphertext.
 * Step 1: Get the histograms of c and z;
 * Step 2: Sort these two histograms;
 * Step 3: Match the frequency between these two;
*/
template<typename DataType>
std::map<DataType, DataType> frequency_ana(const std::vector<DataType>& c,
                                           const std::vector<DataType>& z,
                                           const std::vector<DataType>& c_space,
                                           const std::vector<DataType>& m_space) {
    //Step 1 & 2: Get the histograms of c and z; sort;
    std::map<int, int> hist_c_unsorted = hist<DataType>(c, c_space);
    std::map<int, int> hist_z_unsorted = hist<DataType>(z, m_space);

    // Sort the histogram.
    std::vector<std::pair<int, int>> hist_c(hist_c_unsorted.begin(), hist_c_unsorted.end());
    std::vector<std::pair<int, int>> hist_z(hist_z_unsorted.begin(), hist_z_unsorted.end());
    // lambda
    std::sort(hist_c.begin(), hist_c.end(), [](const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) {
        return lhs.second < rhs.second;
    });
    std::sort(hist_z.begin(), hist_z.end(), [](const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) {
        return lhs.second < rhs.second;
    });

    // Match the frequency.
    std::map<DataType, DataType> alpha;
    for (int i = 0; i < hist_c.size(); i++) {
        if (hist_z[i].first - 1 < 0 || hist_z[i].first - 1 >= m_space.size()) { continue; }
        alpha[c_space[hist_c[i].first - 1]] = m_space[hist_z[i].first - 1];
    }

    return alpha;
}

/*
 * !! Sort the c_space and m_space before doing l_p optimization.
 * 
 */

template<typename DataType>
std::map<DataType, DataType> l_p_optimization(const std::vector<DataType>& c,
                                              const std::vector<DataType>& z,
                                              const std::vector<DataType>& c_space,
                                              const std::vector<DataType>& m_space,
                                              const int& p) {
    std::map<DataType, DataType> alpha;
    std::map<int, int> hist_c = hist<DataType>(c, c_space);
    std::map<int, int> hist_z = hist<DataType>(z, m_space);

    // Build the graph, and note that c_space >= m_space.
    std::vector<std::vector<double>> cost(c_space.size(), std::vector<double>(c_space.size(), 0));
    for (int i = 0; i < cost.size(); i++) {
        for (int j = 0; j < cost[i].size(); j++) {
            cost[i][j] = calculate_cost_DTE(hist_c[i + 1], hist_z[j + 1], p);
        }
    }
    return hungarian_algorithm(cost, c_space, m_space);
}

template<typename DataType>
std::map<DataType, DataType> sorting_attack(const std::vector<DataType>& c, const std::vector<DataType>& m_space) {
    /* Automatically sorted. */
    std::set<DataType> unique_c(c.begin(), c.end());
    std::set<DataType> mk(m_space.begin(), m_space.end());
    std::map<DataType, DataType> alpha;

    for (int i = 0; i < std::min(mk.size(), unique_c.size()); i++) {
            alpha[unique_c[i]] = mk[i];
    }
    return alpha;
}

template<typename DataType>
std::map<DataType, DataType> cumulative_attack(const std::vector<DataType>& c,
                                               const std::vector<DataType>& z, 
                                               const std::vector<DataType>& c_space,
                                               const std::vector<DataType>& m_space,
                                               const double& delta) {
    std::map<int, int> hist_c = hist<DataType>(c, c_space);
    std::map<int, int> hist_z = hist<DataType>(z, m_space);
    std::vector<int> cdf_c = cumulative_density_function<DataType>(hist_c);
    std::vector<int> cdf_z = cumulative_density_function<DataType>(hist_z);

    std::vector<std::vector<double>> cost(c_space.size(), std::vector<double>(c_space.size(), 0));
    for (int i = 0; i < cost.size(); i++) {
        for (int j = 0; j < cost[i].size(); j++) {
            cost[i][j] = calculate_cost_OPE(hist_c[i], cdf_c[i], hist_z[j], cdf_z[j]);
        }
    }
    return hungarian_algorithm(cost, c_space, m_space);
}

template<typename DataType>
std::map<DataType, DataType> OPE_attack(const std::vector<DataType>& c,
                                        const std::vector<DataType>& z, 
                                        const std::vector<DataType>& c_space,
                                        const std::vector<DataType>& m_space,
                                        const double& delta) {
    /* If the OPE column is dense, then we perform sorting attack, which is trivial. */
    if (fabs(delta - 1.0) < DBL_EPSILON) {
        return sorting_attack(c, m_space);
    } else /* If it is delta-dense where delta < 1.0, then we perform cumulative attack. */ {
        return cumulative_attack(c, z, c_space, m_space, delta);
    }
}


/*
 * This would produce partial mapping between C and M because non-crossing matching algorithm requires it.
 */
template<typename DataType>
std::map<DataType, DataType> non_crossing_attack(const std::vector<DataType>& c,
                                                 const std::vector<DataType>& z, 
                                                 const std::vector<DataType>& c_space,
                                                 const std::vector<DataType>& m_space,
                                                 const double& alpha) {
    std::map<DataType, DataType> ans;
    std::map<int, int> hist_c = hist<DataType>(c, c_space);
    std::map<int, int> hist_z = hist<DataType>(z, m_space);
    
    // Using Maximum Weight Non Crossing Matching alogorithm
    // Step 1: label each edge.
    // Step 2: Connect the origin nodes and destination nodes with edges.

    std::vector<std::vector<double>> cost(c_space.size(), std::vector<double>(c_space.size(), 0));
    for (int i = 0; i < cost.size(); i++) {
        for (int j = 0; j < cost[i].size(); j++) {
            cost[i][j] = calculate_cost_non_crossing(hist_c[i], hist_z[j], alpha);
        }
    }
    return maximum_non_crossing_match(cost, c_space, m_space);
}

/*
 * Binomial attack for frequency-hiding encryption schemes.
 * In this attack, the main target is to find a range in ciphertext sequence c of a certain element z_i in auxiliary information sequence z.
 * 
 * ----------------l--------u------------------
 * c1        <zi      zi            >zi        cn
 * l and u should be the start and end point of zi respectively.
 * 
 * Input: 
 *       k: the cardinality of the element set selected from the top K.
 *       d: the confidence parameter.
 * Output: mapping from z_i to range of its corresponding ciphertexts.
 * 
 */
template<typename DataType>
std::map<DataType, std::pair<double, double>> FH_binominal_attack(const std::vector<DataType>& c,
                                                                  const std::vector<DataType>& z,
                                                                  const std::vector<DataType>& c_space,
                                                                  const std::vector<DataType>& m_space,
                                                                  const int& k,
                                                                  const double& d
                                                                ) {
    std::map<DataType, std::pair<double, double>> alpha;
    std::map<int, int> hist_z_unsorted = hist<DataType>(z, m_space);
    // Compute the first k hightest frequency plaintext elements
    std::vector<int> top_k_elements = top_k(z, m_space, hist_z_unsorted, k);

    /* Do mapping. */
    int n = c.size();
    int nz = z.size();
    double epsilon = sqrt(log((1 - d) / 2) / (-2 * n));
    std::vector<double> f_less_than_zi(top_k_elements.size(), 0);
    std::vector<double> f_zi(top_k_elements.size(), 0);

    // Get the probability of drawing elements less than z_i nad f_zi using dynamic programming method.
    for (int i = 1; i < top_k_elements.size(); i++) {
        f_less_than_zi[i] = f_less_than_zi[i - 1] + (hist_z_unsorted[top_k_elements[i - 1]] * 1.0 / nz) ;
    }
    for (int i = 0; i < top_k_elements.size(); i++) {
        f_zi[i] = hist_z_unsorted[top_k_elements[i]] * 1.0 / nz;
    }

    std::vector<std::pair<double, double>> intervals;
    for (int i = 0; i < top_k_elements.size(); i++) {
        /* Compute the probability of drawing z_i */
        intervals.push_back(std::make_pair((f_less_than_zi[i] - epsilon) * n, (f_less_than_zi[i] + f_zi[i] + 2 * epsilon) * n));
    }

    // remove overlaps.
    remove_overlap(intervals, f_zi);

    // create alpha
    for (int i = 0; i < intervals.size(); i++) {
        alpha[m_space[top_k_elements[i] - 1]] = std::make_pair(intervals[i].first, intervals[i].second);
    }

    return alpha;
}