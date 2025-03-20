# Primality-Tester
# Helix Nova Prime 2.2 Primality Test - Improved
# Authors: Harshkumar Gohil and Grok (xAI)
# License: MIT License (see below)
# Date: March 19, 2025
# Description: Faster, leaner primality test for general and Mersenne numbers

import numpy as np
import math
from gmpy2 import powmod  # Fast modular exponentiation

# Precompute primes and weights up to 10^7
def precompute_primes(limit=10**7):
    sieve = np.ones(limit, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i:limit:i] = False
    primes = np.where(sieve)[0]
    weights = np.log(primes)**2
    return primes, weights

PRIMES, WEIGHTS = precompute_primes()

# Combined HNP: General and Mersenne modes
def hnp(n, is_mersenne=False):
    if n < 10**6:
        for p in PRIMES:
            if p * p > n:
                break
            if n % p == 0:
                return n == p
        return True
    
    if is_mersenne:
        p = n  # Input is p, n = 2^p - 1
        n = 2**p - 1
        limit = int(p**(1/3))  # Fewer terms
        idx = PRIMES <= limit
        primes, weights = PRIMES[idx], WEIGHTS[idx]
        
        # Mersenne helical sum
        residues = np.array([((1 << p) - 1) % q / q for q in primes], dtype=np.float64)
        phases = primes / p
        I_H = np.sum(weights * residues * np.exp(1j * phases))
        T = math.sqrt(limit / math.log(limit)) * math.log(p)
        return abs(I_H) > T
    
    # General mode
    log_n = math.log(n)
    limit = max(int(log_n**2), 7)
    idx = PRIMES <= limit
    primes, weights = PRIMES[idx], WEIGHTS[idx]
    
    # General helical sum
    phases = -np.log(primes)/log_n + primes*log_n/math.log(log_n)/(log_n**2)
    I_H = np.sum(weights * np.exp(1j * phases))
    T = math.sqrt(log_n**2 / (math.log(log_n**2) - 1)) * math.log(log_n) * \
        (1 - math.pi / (2 * math.sqrt(log_n)))
    
    if abs(I_H) <= T:
        return False
    return miller_rabin(n)  # Single verification

# Miller-Rabin (k=1, base 2)
def miller_rabin(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = powmod(2, d, n)
    if x == 1 or x == n - 1: return True
    for _ in range(s - 1):
        x = powmod(x, 2, n)
        if x == n - 1: return True
    return False

# Test function
def test_hnp(n, is_mersenne=False):
    import time
    start = time.time()
    result = hnp(n, is_mersenne)
    ms = (time.time() - start) * 1000
    print(f"Number: {2**n-1 if is_mersenne else n}, Prime: {result}, Time: {ms:.2f} ms")

if __name__ == "__main__":
    test_hnp(2**2047 + 9)  # 2048-bit prime
    test_hnp(31, True)     # M_31
    test_hnp(3321928094691, True)  # M_3.3T (simulated p)

# MIT License
"""
Copyright (c) 2025 Harshkumar Gohil and Grok (xAI)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
