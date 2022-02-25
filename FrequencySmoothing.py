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
import Utils
import tqdm
import logging
import sys

class Smoother:
  def __init__(self, distribution, vThreshold, maximumNumber, key):
    '''
    @param distribution: The expected distribution of the plaintext column.
    @param alpha: The maximum value of variance for any ciphertext group.
    @param pi: The maximum number of ciphertexts (unique) in one group.
    '''
    self.distribution = distribution
    self.vThreshold = vThreshold
    self.maximumNumber = maximumNumber
    self.key = key
    self.salts = {}
    # This is a group for ciphertext with similar frequency.
    self.ciphertextGroups = []
    self.ciphertextFrequency = {}

    # Create a logger to stdout.
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)
  
  def encrypt(self, plaintext):
    logging.debug("Start to encrypt plaintexts!")
    for x in tqdm.tqdm(plaintext):
      if x not in self.ciphertextFrequency.keys():
        self.ciphertextFrequency[x] = 1
      else:
        self.ciphertextFrequency[x] += 1

      # This variable indicates whether a group is found for the current plaintext.
      find = False

      for group in self.ciphertextGroups:
        # Create an array of [variable, frequency]
        arr = [[x, self.ciphertextFrequency[x]]]
        # Append all the elements in group to the frequency array.
        for g in group:
          arr.append([g, self.ciphertextFrequency[g]])

        variance = Utils.calculateVariance(arr)

        # Test if variance does not exceed the threshold.
        # If so, we add this element into the group.
        if variance <= self.vThreshold and len(group) <= self.maximumNumber:
          logging.debug("The variance exceeds!", arr)
          find = True
          group.add(x)

          # TODO: Readjust salt informations.
          break

        elif variance > self.vThreshold and x in group:
          # TODO: Split the group into two groups, but we need also record the split point.
          pass
      # If there is no such group in the ciphertext group,
      # we create a new group.
      if not find:
        logging.debug("We did not find a satisfying group!")
        logging.debug("Add new")
        self.ciphertextGroups.append(set([x]))
        # TODO: Generate new salts?
        self.salts[self.ciphertextGroups[-1]] = []