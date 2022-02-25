"""
 Copyright (c) 2022 Haobin Chen

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
 """
import numpy as np
import Hungarian
import Utils

# If |C_k| != |M_k|, we simply pad M_k with dummy values.
def padAuxiliary(auxiliaryDataset, ciphertext):
  m, n = len(np.unique(auxiliaryDataset)), len(np.unique(ciphertext))
  A, Z = np.array(["A","Z"]).view("int32") 
  frequency = np.bincount(auxiliaryDataset)
  mostFrequent, leastFrequent = np.max(frequency), np.min(frequency)

  NO_CODES = 100
  LEN = 20
  # We assume that C is al least equal to M.
  assert(n >= m)
  if m != n:
    diff = np.abs(m - n)
    
    # Pad
    for i in range(diff):
      randomVal = np.random.randint(low=A,high=Z,size=NO_CODES*LEN,dtype="int32").view(f"U{LEN}")

      randomFrequency = np.random.randint(low=leastFrequent, high=mostFrequent)
      for j in range(randomFrequency):
        auxiliaryDataset.append(randomVal)

def prepareData(x, y):
  # First convert them to strings.
  x = [str(_) for _ in x]
  padAuxiliary(x, y)

  # Then calculate the histograms for auxiliary dataset and the ciphertext column.
  histX = Utils.calculateHistogram(x)
  histY = Utils.calculateHistogram(y)

  return histX, histY

# This function analyzes the frequency information and outputs a mapping from message to the ciphertext
def frequencyAnalysis(auxiliaryDataset, ciphertext):
  histAuxiliary, histCiphertext = prepareData(auxiliaryDataset, ciphertext)
  
  # Map the most frequent ciphertext to the most frequent plaintext.
  # We assume that the support for these two sets are the same.
  guess = []
  for i in range(len(histAuxiliary)):
    item = (histCiphertext[i][0], histAuxiliary[i][0])
    guess.append(item)

  return guess

# This algorithm implements the linear optimization problem defined in NKW's paper.
# This algorithm also utilizes the Hungarian algorithm for bipartite graph matching
# problem that minimized the cost for two nodes.
# The idea behind l_p optimization is to find a minimized cost for matching the cipher-
# text to its corresponding plaintext by the l_p norm.
def LPOpmization(auxiliaryDataset, ciphertext, p=2):
  histAuxiliary, histCiphertext = prepareData(auxiliaryDataset, ciphertext)

  print(histAuxiliary, histCiphertext)

  # Start to compute the l_p optimization.
  dimension = len(histAuxiliary)
  costMatrix = [[0 for i in range(dimension)] for j in range(dimension)]

  # Calculate the cost matrix.
  for i in range(len(costMatrix)):
    for j in range(len(costMatrix[0])):
      LPDistance = np.abs(histCiphertext[i][1] - histAuxiliary[j][1])
      costMatrix[i][j] = LPDistance ** p

  guess = Hungarian.hungarian(np.matrix(np.array(costMatrix).reshape(dimension, dimension)))
  ans = []
  # Transform the guess array.
  for i in range(len(guess[0])):
     ans.append((histCiphertext[guess[0][i]][0], histAuxiliary[guess[0][i]][0]))
  return ans

# This algorithm implements the MLE attack written in the paper The Tao of Inference in 
# Privacy-Preserving Databases by Grubbs et al, in VLDB 2018.
# Other literature has pointed out that frequency analysis is also MLE.
def maximumLikelihoodEstimationAttack(auxiliaryDataset, ciphertext):
  histAuxiliary, histCiphertext = prepareData(auxiliaryDataset, ciphertext)