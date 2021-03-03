# FrequencyAttack
The attack algorithms are based on Naveed et al's paper regarding inference attack.
In order to evaluate the efficiency as well as the security of frequency-smoothing technique against traditional *inference attack*, we run experiments to see if distortion works for the traditional frequency analysis and also newly proposed algorithms including linear opitmization methods.
Maximun likelihood estimate attack is for cross-column attack which may exploit correlation between columns but is equivalent for the case of single column attack, so we do not perform MLE attack on single DE column.
Any determinsitc encryption is valid for the test, and for convenience, therefore, we use self-made encryption algorithm to see how frequency leakage will help an adversary to recover plaintexts of a certain column.

All the parameters are described in the comments.
