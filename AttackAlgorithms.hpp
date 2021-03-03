#ifndef _ATTACK_ALGORITHMS_HPP_
#define _ATTACK_ALGORITHMS_HPP_

#include <iostream>
#include <stdlib.h>
#include <cfloat> // for DBL_MAX
#include <cmath>  // for fabs()
#include <vector>
#include <map>

class HungarianAlgorithm {
public:
	HungarianAlgorithm();
	~HungarianAlgorithm();
	double Solve(std::vector<std::vector<double>>& DistMatrix, std::vector<int>& Assignment);
private:
	void assignment_optimal(int *assignment, double *cost, double *distMatrix, int nOfRows, int nOfColumns);
	void buildassignment_vector(int *assignment, bool *starMatrix, int nOfRows, int nOfColumns);
	void compute_assignment_cost(int *assignment, double *cost, double *distMatrix, int nOfRows);
	void step2a(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step2b(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step3(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step4(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col);
	void step5(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
};


HungarianAlgorithm::HungarianAlgorithm(){}
HungarianAlgorithm::~HungarianAlgorithm(){}


//********************************************************//
// A single function wrapper for solving assignment problem.
//********************************************************//
double HungarianAlgorithm::Solve(std::vector <std::vector<double> >& DistMatrix, std::vector<int>& Assignment) {
	unsigned int nRows = DistMatrix.size();
	unsigned int nCols = DistMatrix[0].size();

	double *distMatrixIn = new double[nRows * nCols];
	int *assignment = new int[nRows];
	double cost = 0.0;

	// Fill in the distMatrixIn. The index is "i + nRows * j".
	// Here the cost matrix of size MxN is defined as a double precision array of N*M elements. 
	for (unsigned int i = 0; i < nRows; i++) {
		for (unsigned int j = 0; j < nCols; j++) {
			distMatrixIn[i + nRows * j] = DistMatrix[i][j];
		}
	}
	// call solving function
	assignment_optimal(assignment, &cost, distMatrixIn, nRows, nCols);

	Assignment.clear();
	for (unsigned int r = 0; r < nRows; r++)
		Assignment.push_back(assignment[r]);

	delete[] distMatrixIn;
	delete[] assignment;
	return cost;
}

// Solve optimal solution for assignment problem using Munkres algorithm, also known as Hungarian Algorithm.
void HungarianAlgorithm::assignment_optimal(int *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns) {
	double *distMatrix, *distMatrixTemp, *distMatrixEnd, *columnEnd, value, minValue;
	bool *coveredColumns, *coveredRows, *starMatrix, *newStarMatrix, *primeMatrix;
	int nOfElements, minDim, row, col;

	/* initialization */
	*cost = 0;
	for (row = 0; row<nOfRows; row++)
		assignment[row] = -1;

	/* generate working copy of distance Matrix */
	/* check if all matrix elements are positive */
	nOfElements = nOfRows * nOfColumns;
	distMatrix = (double *)malloc(nOfElements * sizeof(double));
	distMatrixEnd = distMatrix + nOfElements;

	for (row = 0; row<nOfElements; row++) {
		value = distMatrixIn[row];
		if (value < 0)
			std::cerr << "All matrix elements have to be non-negative." << std::endl;
		distMatrix[row] = value;
	}


	/* memory allocation */
	coveredColumns = (bool *)calloc(nOfColumns, sizeof(bool));
	coveredRows = (bool *)calloc(nOfRows, sizeof(bool));
	starMatrix = (bool *)calloc(nOfElements, sizeof(bool));
	primeMatrix = (bool *)calloc(nOfElements, sizeof(bool));
	newStarMatrix = (bool *)calloc(nOfElements, sizeof(bool)); /* used in step4 */

	/* preliminary steps */
	if (nOfRows <= nOfColumns) {
		minDim = nOfRows;

		for (row = 0; row < nOfRows; row++) {
			/* find the smallest element in the row */
			distMatrixTemp = distMatrix + row;
			minValue = *distMatrixTemp;
			distMatrixTemp += nOfRows;
			while (distMatrixTemp < distMatrixEnd) {
				value = *distMatrixTemp;
				if (value < minValue)
					minValue = value;
				distMatrixTemp += nOfRows;
			}

			/* subtract the smallest element from each element of the row */
			distMatrixTemp = distMatrix + row;
			while (distMatrixTemp < distMatrixEnd) {
				*distMatrixTemp -= minValue;
				distMatrixTemp += nOfRows;
			}
		}

		/* Steps 1 and 2a */
		for (row = 0; row<nOfRows; row++)
			for (col = 0; col<nOfColumns; col++)
				if (fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON)
					if (!coveredColumns[col]) {
						starMatrix[row + nOfRows * col] = true;
						coveredColumns[col] = true;
						break;
					}
	}
	else /* if(nOfRows > nOfColumns) */ {
		minDim = nOfColumns;

		for (col = 0; col<nOfColumns; col++) {
			/* find the smallest element in the column */
			distMatrixTemp = distMatrix + nOfRows * col;
			columnEnd = distMatrixTemp + nOfRows;

			minValue = *distMatrixTemp++;
			while (distMatrixTemp < columnEnd) {
				value = *distMatrixTemp++;
				if (value < minValue)
					minValue = value;
			}

			/* subtract the smallest element from each element of the column */
			distMatrixTemp = distMatrix + nOfRows * col;
			while (distMatrixTemp < columnEnd)
				*distMatrixTemp++ -= minValue;
		}

		/* Steps 1 and 2a */
		for (col = 0; col<nOfColumns; col++){
			for (row = 0; row<nOfRows; row++){
				if (fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON) {
					if (!coveredRows[row]) {
						starMatrix[row + nOfRows * col] = true;
						coveredColumns[col] = true;
						coveredRows[row] = true;
						break;
					}
				}
			}
		}
		for (row = 0; row < nOfRows; row++) {
			coveredRows[row] = false;
		}
	}

	/* move to step 2b */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

	/* compute cost and remove invalid assignments */
	compute_assignment_cost(assignment, cost, distMatrixIn, nOfRows);

	/* free allocated memory */
	free(distMatrix);
	free(coveredColumns);
	free(coveredRows);
	free(starMatrix);
	free(primeMatrix);
	free(newStarMatrix);

	return;
}

/********************************************************/
void HungarianAlgorithm::buildassignment_vector(int *assignment, bool *starMatrix, int nOfRows, int nOfColumns) {
	int row, col;

	for (row = 0; row < nOfRows; row++) {
		for (col = 0; col < nOfColumns; col++){
			if (starMatrix[row + nOfRows * col]) {
				assignment[row] = col;
				break;
			}
		}
	}
}

/********************************************************/
void HungarianAlgorithm::compute_assignment_cost(int *assignment, double *cost, double *distMatrix, int nOfRows) {
	int row, col;

	for (row = 0; row < nOfRows; row++) {
		col = assignment[row];
		if (col >= 0) {
			*cost += distMatrix[row + nOfRows * col];
		}
	}
}

/********************************************************/
void HungarianAlgorithm::step2a(int *assignment,
								double *distMatrix,
								bool *starMatrix,
								bool *newStarMatrix,
								bool *primeMatrix,
								bool *coveredColumns,
								bool *coveredRows,
								int nOfRows,
								int nOfColumns,
								int minDim) {
	bool *starMatrixTemp, *columnEnd;
	int col;

	/* cover every column containing a starred zero */
	for (col = 0; col < nOfColumns; col++) {
		starMatrixTemp = starMatrix + nOfRows * col;
		columnEnd = starMatrixTemp + nOfRows;
		while (starMatrixTemp < columnEnd) {
			if (*starMatrixTemp++) {
				coveredColumns[col] = true;
				break;
			}
		}
	}

	/* move to step 3 */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void HungarianAlgorithm::step2b(int *assignment,
							    double *distMatrix,
								bool *starMatrix,
								bool *newStarMatrix,
								bool *primeMatrix,
								bool *coveredColumns,
								bool *coveredRows,
								int nOfRows,
								int nOfColumns,
								int minDim) {
	int col, nOfCoveredColumns;

	/* count covered columns */
	nOfCoveredColumns = 0;
	for (col = 0; col < nOfColumns; col++) {
		if (coveredColumns[col]) {
			nOfCoveredColumns++;
		}
	}

	if (nOfCoveredColumns == minDim) {
		/* algorithm finished */
		buildassignment_vector(assignment, starMatrix, nOfRows, nOfColumns);
	} else {
		/* move to step 3 */
		step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}

}

/********************************************************/
void HungarianAlgorithm::step3(int *assignment,
							   double *distMatrix,
							   bool *starMatrix,
						       bool *newStarMatrix,
							   bool *primeMatrix,
							   bool *coveredColumns,
							   bool *coveredRows,
							   int nOfRows,
							   int nOfColumns,
							   int minDim) {
	bool zerosFound;
	int row, col, starCol;

	zerosFound = true;
	while (zerosFound) {
		zerosFound = false;
		for (col = 0; col < nOfColumns; col++) {
			if (!coveredColumns[col]) {
				for (row = 0; row < nOfRows; row++) {
					if ((!coveredRows[row]) && (fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON)) {
						/* prime zero */
						primeMatrix[row + nOfRows * col] = true;

						/* find starred zero in current row */
						for (starCol = 0; starCol < nOfColumns; starCol++) {
							if (starMatrix[row + nOfRows * starCol]) {
								break;
							}
						}
						if (starCol == nOfColumns) /* no starred zero found */ {
							/* move to step 4 */
							step4(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim, row, col);
							return;
						} else {
							coveredRows[row] = true;
							coveredColumns[starCol] = false;
							zerosFound = true;
							break;
						}
					}
				}
			}
		}
	}

	/* move to step 5 */
	step5(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void HungarianAlgorithm::step4(int *assignment,
							   double *distMatrix,
							   bool *starMatrix,
							   bool *newStarMatrix,
							   bool *primeMatrix,
							   bool *coveredColumns,
							   bool *coveredRows,
							   int nOfRows,
							   int nOfColumns,
							   int minDim,
							   int row,
							   int col) {
	int n, starRow, starCol, primeRow, primeCol;
	int nOfElements = nOfRows * nOfColumns;

	/* generate temporary copy of starMatrix */
	for (n = 0; n < nOfElements; n++) {
		newStarMatrix[n] = starMatrix[n];
	}
		
	/* star current zero */
	newStarMatrix[row + nOfRows * col] = true;

	/* find starred zero in current column */
	starCol = col;
	for (starRow = 0; starRow < nOfRows; starRow++) {
		if (starMatrix[starRow + nOfRows * starCol]) {
			break;
		}
	}

	while (starRow < nOfRows) {
		/* unstar the starred zero */
		newStarMatrix[starRow + nOfRows * starCol] = false;

		/* find primed zero in current row */
		primeRow = starRow;
		for (primeCol = 0; primeCol < nOfColumns; primeCol++) {
			if (primeMatrix[primeRow + nOfRows * primeCol]) {
				break;
			}
		}

		/* star the primed zero */
		newStarMatrix[primeRow + nOfRows * primeCol] = true;

		/* find starred zero in current column */
		starCol = primeCol;
		for (starRow = 0; starRow < nOfRows; starRow++) {
			if (starMatrix[starRow + nOfRows * starCol]) {
				break;
			}
		}
	}

	/* use temporary copy as new starMatrix */
	/* delete all primes, uncover all rows */
	for (n = 0; n < nOfElements; n++) {
		primeMatrix[n] = false;
		starMatrix[n] = newStarMatrix[n];
	}
	for (n = 0; n < nOfRows; n++) {
		coveredRows[n] = false;
	}
		
	/* move to step 2a */
	step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void HungarianAlgorithm::step5(int *assignment,
						       double *distMatrix,
							   bool *starMatrix,
							   bool *newStarMatrix,
							   bool *primeMatrix,
							   bool *coveredColumns,
							   bool *coveredRows,
							   int nOfRows,
							   int nOfColumns,
							   int minDim) {
	double h, value;
	int row, col;

	/* find smallest uncovered element h */
	h = DBL_MAX;
	for (row = 0; row < nOfRows; row++){
		if (!coveredRows[row]){
			for (col = 0; col < nOfColumns; col++){
				if (!coveredColumns[col]) {
					value = distMatrix[row + nOfRows * col];
					if (value < h) {
						h = value;
					}
				}
			}
		}
	}
	/* add h to each covered row */
	for (row = 0; row < nOfRows; row++) {
		if (coveredRows[row]) {
			for (col = 0; col < nOfColumns; col++) {
				distMatrix[row + nOfRows * col] += h;
			}
		}
	}
	/* subtract h from each uncovered column */
	for (col = 0; col < nOfColumns; col++) {
		if (!coveredColumns[col]) {
			for (row = 0; row < nOfRows; row++) {
				distMatrix[row + nOfRows * col] -= h;
				}
		}
	}
	/* move to step 3 */
	step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}


template<typename DataType>
std::map<DataType, DataType> maximum_non_crossing_match(std::vector<std::vector<double>> cost,
                                                        const std::vector<DataType>& c_space,
                                                        const std::vector<DataType>& m_space) {
    std::map<DataType, DataType> alpha;
    
    int remain = cost.size() * cost[0].size();
    while (remain > 0) {
        std::pair<int, int> where_is_k;
        for (int i = 0; i < cost.size(); i++) {
            double k = cost[i][0];
            for (int j = 0; j < cost[i].size(); j++) {
                if (cost[i][j] > k) {
                    k = cost[i][j];
                    where_is_k.first = i;
                    where_is_k.second = j;
                }
            }
        }
        alpha[c_space[where_is_k.first]] = m_space[where_is_k.second];

        for (int j = 0; j < cost[0].size(); j++) {
            if (cost[where_is_k.first][j] != -1) {
                remain --;
                cost[where_is_k.first][j] = -1;
            }
            
        }

        for (int i = 0; i < cost.size(); i++) {
            if (cost[i][where_is_k.second] != -1) {
                remain --;
                cost[i][where_is_k.second] = -1;
            }
            
        }

        // delete crossing edges.
        for (int i = 0; i < where_is_k.first; i++) {
            for (int j = where_is_k.second + 1; j < cost[i].size(); j++) {
                if (cost[i][j] != -1) {
                    remain--;
                    cost[i][j] = -1;
                }
            }
        }
        for (int i = where_is_k.first + 1; i < cost.size(); i++) {
            for (int j = 0; j < where_is_k.second; j++) {
                if (cost[i][j] != -1) {
                    remain--;
                    cost[i][j] = -1;
                }
            }
        }
    }

    return alpha;
}

void remove_overlap(std::vector<std::pair<double, double>>& intervals,
                    const std::vector<double>& f_zi) {
    sort(intervals.begin(), intervals.end(), [](const std::pair<double, double>& lhs, const std::pair<double, double>& rhs) {
        return lhs.first < rhs.first || (lhs.first == rhs.first && lhs.second < rhs.second);
    });
    for (int i = 1; i < intervals.size(); i++) {
        /*
         * ---------------|--fraction_first--/-----------------| i-1
         *                |------------------/-fraction_second-|----------------- i
         */
        if ((intervals[i].first - intervals[i - 1].second) <= DBL_EPSILON) {
            double length = intervals[i - 1].second - intervals[i].first;
            double fraction_first = (f_zi[i - 1] * 1.0 / (f_zi[i] + f_zi[i - 1])) * length;
            double fraction_second = (f_zi[i] * 1.0 / (f_zi[i] + f_zi[i - 1])) * length;
            intervals[i - 1].second -= fraction_second;
            intervals[i].first += fraction_first;
        }
    }
}

/*
 * Input: 
 *      data: the content of the data;
 *      dataset: the space of the data, i.e., the domain.
 * Output: 
 *      histogram:
*/

template<typename DataType>
std::map<int, int> hist(const std::vector<DataType>& data, const std::vector<DataType>& dataset) {
    std::map<int, int> hist;

    for (int i = 0; i < data.size(); i++) {
        // Time complexity O(log(dataset.size()))
        auto position_in_dataset = find(dataset.begin(), dataset.end(), data[i]);
        int index = distance(dataset.begin(), position_in_dataset) + 1;
        hist[index] += 1;
    }

    // Add zero items.
    for (int i = 0; i < dataset.size(); i++) {
        if (hist.end() == hist.find(i + 1)) {
            hist[i + 1] = 0;
        }
    }
    return hist;
}

template<typename DataType>
std::vector<int> cumulative_density_function(std::map<int, int> hist) {
    if (hist.empty()) {
        return {};
    }

    std::vector<int> cdf(hist.size(), 0);
    cdf[0] = hist[1];
    for (int i = 1; i < hist.size(); i++) {
        for (int j = 0; j <= i; j++) {
            cdf[i] = cdf[i - 1] + hist[i + 1];
        }
    }
    return cdf;
}

template<typename DataType>
std::vector<int> top_k(const std::vector<DataType>& z,
                       const std::vector<DataType>& m_space,
                       const std::map<int, int>& hist_z_unsorted,
                       const int& k) {
    std::vector<int> top_k_elements;
    std::vector<std::pair<int, int>> hist_z(hist_z_unsorted.begin(), hist_z_unsorted.end());
    sort(hist_z.begin(), hist_z.end(), [](const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) {
        return lhs.second > rhs.second;
    });
    for (int i = 0; i < k; i++) {
        top_k_elements.push_back(hist_z[i].first); // Stores the index.
    }
    reverse(top_k_elements.begin(), top_k_elements.end());
    return top_k_elements;
}

template<typename DataType>
std::map<DataType, DataType> hungarian_algorithm(std::vector<std::vector<double>>& cost, 
                                                 const std::vector<DataType>& c_space,
                                                 const std::vector<DataType>& m_space) {
    std::vector<int> assignment;
    HungarianAlgorithm ha;
    ha.Solve(cost, assignment);

    std::map<DataType, DataType> alpha;
    for (int i = 0; i < assignment.size(); i++) {
        alpha[c_space[i]] = m_space[assignment[i]];
    }
    return alpha;
}

double calculate_cost_DTE(const double& from, const double& to, const int& p) {
    return pow(fabs(from - to), p);
}

double calculate_cost_OPE(const double& psi, const double& phi, const double& pi, const double& mu) {
    return pow(fabs(psi - pi), 2) + pow(fabs(phi - mu), 2);
}

double calculate_cost_non_crossing(const double& from, const double& to, const double& alpha) {
    return alpha - fabs(from - to);
}

#endif