import math

import Crypto.Util.number as number


# Find the first prime number that is less than a given number
for i in range(int(math.pow(10, 8))-1, 1, -2):
    if number.isPrime(i):
        print(i)
        break

