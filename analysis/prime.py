import math

import Crypto.Util.number as number


# Find the first prime number that is less than a given number
# for i in range(int(math.pow(2, 23))+1, 1000000, 2):
for i in range(10000001, 20000000, 2):
    if number.isPrime(i):
        print(i)
        break

