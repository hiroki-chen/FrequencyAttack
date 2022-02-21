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
import hashlib
import base64

from Crypto.Cipher import AES
from Crypto import Random
from Crypto.Util.Padding import pad

# Draws a vector of random values from a distribution specified by the type.
def randomMessage(length, type, extraMessage=...):
  if type == 0:
    randomValues = [np.random.uniform(extraMessage['lower'], extraMessage['upper']) for _ in range(length)]
    return randomValues

  elif type == 1:
    randomValues = np.random.normal(extraMessage['mean'], extraMessage['variance'], (length, 1))
    return randomValues

def encryptAESECB(plaintext, masterKey):
  length = len(plaintext)
  # Convert to string array.
  plaintext = [str(_) for _ in plaintext]
  # Get a key.
  hashKey = hashlib.sha256(masterKey.encode()).digest()
  print(hashKey)

  # AES ECB mode encryption
  map = []
  for i in range(length):
    AESContext = AES.new(hashKey, AES.MODE_ECB)
    transformedRaw = pad(plaintext[i].encode(), AES.block_size)
    rawCiphertext = AESContext.encrypt(transformedRaw)
    ciphertext = base64.b64encode(rawCiphertext)
    map.append((ciphertext, plaintext[i]))
  
  return map

def encryptAESCBC(plaintext, masterKey):
  length = len(plaintext)
  # Convert to string array.
  plaintext = [str(_) for _ in plaintext]
  # Get a key.
  hashKey = hashlib.sha256(masterKey.encode()).digest()
  print(hashKey)
  IV = [Random.new().read(AES.block_size) for _ in range(length)]
  
  # AES CBC mode encryption
  map = []
  for i in range(length):
    AESContext = AES.new(hashKey, AES.MODE_CBC, IV[i])
    transformedRaw = pad(plaintext[i].encode(), AES.block_size)
    rawCiphertext = AESContext.encrypt(transformedRaw)
    ciphertext = base64.b64encode(IV[i] + rawCiphertext)
    map.append((ciphertext, plaintext[i]))
  
  return map

def calculateAccuracy(guess, real):
  intersect = set.intersection(set(guess), set(real))
  return float(len(intersect) / len(guess))

if __name__ == '__main__':
  print(encryptAESECB(['12jfisa', 'fas', 'fas'], 'fuck'))