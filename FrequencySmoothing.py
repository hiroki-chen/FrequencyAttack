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
import string

'''
salts = <groupID> -> array: tuple <salt, frequency>
plaintextToSetId = <plaintext> -> setId
'''

class Smoother:
  def __init__(self, distribution, vThreshold, maximumNumber, key, p, saltLen):
    '''
    @param distribution: The expected distribution of the plaintext column.
    @param alpha: The maximum value of variance for any ciphertext group.
    @param pi: The maximum number of ciphertexts (unique) in one group.
    '''
    self.distribution = distribution
    self.vThreshold = vThreshold
    self.maximumNumber = maximumNumber
    self.key = key
    self.groupID = 0
    self.p = p
    self.saltLen = saltLen
    # Maps from set id to salt array.
    self.salts = {}
    # This is a group for ciphertext with similar frequency.
    self.plaintextGroups = {}
    self.plaintextFrequency = {}
    self.plaintextToSetID = {}

    # Create a logger to stdout.
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)

    self.ALPHABET = np.array(list(string.ascii_letters))
  
  def encrypt(self, plaintext):
    logging.debug("Start to encrypt plaintexts!")
    ciphertexts = []

    for x in tqdm.tqdm(plaintext):
      if x not in self.plaintextFrequency.keys() and \
         x not in self.plaintextToSetID.keys():
        # If not found we also need to create a new group for x.
        self.plaintextToSetID[x] = self.groupID
        self.plaintextFrequency[x] = 1
        self.plaintextGroups[self.groupID] = []
        self.salts[self.groupID] = []
        self.groupID = self.groupID + 1

        logging.debug("Find a new plaintext: {}".format(x))
      else:
        self.plaintextFrequency[x] += 1

      # Find all groups for this plaintext x.
      print(self.plaintextToSetID)
      groupID = self.plaintextToSetID[x]
      allGroups = self.plaintextGroups[groupID]

      # This variable indicates whether there is a group that can contain x.
      find = False
      for group in allGroups:
        arr = [[x, self.plaintextFrequency[x]]]
        for g in group:
          arr.append([g, self.plaintextFrequency[g]])
        variance = Utils.calculateVariance(arr)
        # Calculate variance and check it.
        if variance <= self.vThreshold:
          if x not in g and len(g) + 1 <= self.maximumNumber:
            g.append([x])
          find = True
          logging.debug("We have found the group!")
          break
      
      if not find:
        allGroups.append(set([x]))
        logging.debug("Not found for {}".format(x))

      salts = self.salts[groupID]
      # Start to encrypt the plaintext.
      sampleNewSalt = np.random.binomial(1, self.p)
      ans = ''
      if sampleNewSalt:
        newSalt = ''.join(np.random.choice(self.ALPHABET, size=self.saltLen))
        salts.append([newSalt, 1])
        logging.debug("Generated a new salt: {}".format(newSalt))

      else:
        np.random.permutation(salts)
        for salt in salts:
          if (Utils.checkSalt(salt, salts)):
            ans = salt
            # Remember to increment the frequency of the selected salt.
            salt[1] = salt[1] + 1
            break
            
        if ans == '':
          newSalt = ''.join(np.random.choice(self.ALPHABET, size=self.saltLen))
          salts.append([newSalt, 1])
          logging.debug("Generated a new salt: {}".format(newSalt))
          ans = newSalt
        
      # Use this salt to encrypt the plaintext.
      ciphertext = Utils.encryptAESECB([x + newSalt], self.key)
      ciphertexts.append(ciphertext)

    return ciphertexts

  def checkSalt(self, salt, salts, plaintextNum):
    saltFrequency = salt[1]
    # We should check if the salt can efficiently 'smooth' the plaintext's frequency.
    
    