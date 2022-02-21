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
import InferenceAttack as ia
import Utils

if __name__ == '__main__':
  real = Utils.encryptAESECB([1,1,2,5,5,5], 'fuck')
  print(real)
  ciphertext = [x[0] for x in real]
  guess = ia.frequencyAnalysis([1,1,1,2,5,5,5,5,5], ciphertext)
  print(guess)
  accuracy = Utils.calculateAccuracy(guess, real)
  print(accuracy)
  ia.LPOpmization([1,1,1,2,5,5,5,5,5], ciphertext)