# preparations
import time


# Calculating residuum a in n, so a number in the Montgomery space
# such as a = a * r mod n

def calculate_residuum_n(a, r, modulo):
    a = a * r % modulo
    return a


# Additional function that can calculate the inverse of the negative modulus
# basing on a Extended Euclidean algorithm and fulfilling the equation:
# r * r^(-1) - n * n' = 1, giving n' = (-n)^(-1)

def neg_inv(n, bit_width, base):
    def_base = 2 ** (bit_width - 1)
    neg_inverse = base - n ** (def_base - 1) % base
    return neg_inverse


# Alternative version that calculates n'

def alternative_neg_inv(n, r):
    r_inverse = 0
    n_prim = 0
    for i in range(0, r):
        if (r * i) % n == 1:
            r_inverse = i
            break
    for i in range(0, n):
        if r * r_inverse - n * i == 1:
            n_prim = i
            break
    return n_prim


# Function monPro gives the direct product of the modular multiplication
# where multiplier and multiplicand are already converted into
# the Montgomery space
# Elements of the func:
#   a,b - multiplier and multiplicand in Montgomery space
#   n - modulus
#   r - accordingly to the documentation, a variable that
#       fulfills the r = 2^k <=> 2^(k-1) <= n < 2^k
#   bit_width - additional variable defining exponent of r

def monPro(a, b, r, n, bit_width):
    t = a * b % r
    m = (t * neg_inv(n,bit_width, r)) % r
    s = (a * b + m * n)/r
    if s >= n:
        return s - n
    else:
        return s


# Function modMul is the main engine of the modular
# multiplication algorithm, where multiplier and multiplicand
# are firstly converted in the Montgomery space
# and then the monPro func is called where the direct
# product of the multiplication is given
# Lastly the function uses the fast - inverse algorithm
# that re-converts the result from Montgomery space.

def modMul(x, y, n):
    bit_width = calc_r(n)
    r = 2 ** bit_width
    a = x * r % n
    b = y * r % n
    montg_product = monPro(a, b, r, n,bit_width)
    product = monPro(montg_product, 1, r, n, bit_width)
    return product


# Function calculating exponent (k) of the r = 2^k
# accordingly to the restriction: r = 2^k <=> 2^(k-1) <= n < 2^k

def calc_r(n):
    for i in range(0, 32):
        r = 2**i
        if (r > n) and (r/2 <= n):
            return i


# Main function calculating the (a ^ e mod n)
# basing on the Fast Montgomery Pipelined Exponentiation algorithm.
# The only difference in the elements of the func is the parameter e:
# Exponent of the main task

def montg_Exp(x,e,n):
    bit_width = calc_r(n)
    r = 2**bit_width
    A = x * r % n
    X = 1 * r % n
    e_binary = format(e, "b")
    for digits in e_binary:
        X = monPro(X, X, r, n, bit_width)
        if digits == '1':
            X = monPro(A, X, r, n, bit_width)
    X = monPro(X, 1, r, n, bit_width)
    return X


# Additional function that generates elements
# useful in the factorization of even modulus when
# exponentiation is concerned.
# Given even modulus the n = q * 2^j
# where q will be the n - representative
# after a series of binary right-shifts
# performed on n as long as n is even
# j will be the binary notation exponent
# Example:
#   Given n such as 12 - (1100)_bin
#   Result of n shifts as long as n is even - (0011)_bin * 2^2
#   Result q = 3, j = 2

def n_even(n):
    j = 0
    for i in range(0, 32):
        n >>= 1
        j += 1
        if n % 2 != 0:
            return n, j


# montg_notEven function calculates the
# Montgomery exponentiation scheme when n is even
# it uses the montg_Exp function and a specified algorith
# which bases on the results of n_even function
# so the factorized elements of the n modulus:
# q - odd modulus gained after the series of right-shifts
# j - exponent of the binary notation where:
#           n = q * 2^j


def montg_notEven(a, e, n):
    q, j = n_even(n)
    x1 = montg_Exp(a, e, q)
    quotient = 2**j
    x2 = pow(a, e) % quotient
    q_inverse = 0
    for i in range(0, q):
        if (q * i) % quotient == 1:
            q_inverse = i
            break
    y = (x2-x1)*q_inverse % quotient
    return x1 + q*y


def callMenu():
    while True:
        print("=====Montgomery Pipeline Multiplication and Exponentiation=====\n"
              "1. Multiplication\n"
              "2. Exponentiation\n"
              "3. Exit\n"
              "Choice: ")
        try:
            choice = int(input())
        except:
            print("Input an integer number\n")
        else:
            if choice == 1:
                print("Calculator represents the x * y mod n")
                try:
                    x = int(input("Insert x: "))
                    y = int(input("Insert y: "))
                    n = int(input("Insert n: "))
                except:
                    print("x,y and n must be Integers!!\n")
                else:
                    if n % 2 == 0:
                        print("In order for the Montgomery modular multiplication to work the modulus must be odd!")
                        holdback = input("\nPress Enter to continue:\n")
                    else:
                        print('Expected (x * y mod n) result:', x * y % n)
                        start_montg = time.time()
                        wynik = modMul(x, y, n)
                        end_montg = time.time()
                        print('Montgomery Scheme result:', int(wynik))
                        print('Elapsed time of the algorithm [s]: ', round(end_montg - start_montg,5))
                        holdback = input("\nPress Enter to continue:\n")
            if choice == 2:
                print("Calculator represents the a^e mod n")
                try:
                    a = int(input("Insert a: "))
                    e = int(input("Insert e (power): "))
                    n = int(input("Insert n: "))
                except:
                    print("a,e and n must be Integers!!\n")
                else:
                    if n % 2 == 0:
                        print('Expected (a^e mod n) result:',pow(a,e) % n)
                        start_montg = time.time()
                        wynik = montg_notEven(a,e,n)
                        end_montg = time.time()
                        print('Montgomery Scheme result:', int(wynik))
                        print('Elapsed time of the algorithm [s]: ', round(end_montg - start_montg,5))
                        holdback = input("\nPress Enter to continue:\n")
                    else:
                        print('Expected (a^e mod n) result:', pow(a, e) % n)
                        start_montg = time.time()
                        wynik = montg_Exp(a, e, n)
                        end_montg = time.time()
                        print('Montgomery Scheme result:', int(wynik))
                        print('Elapsed time of the algorithm [s]: ', round(end_montg - start_montg,5))
                        holdback = input("\nPress Enter to continue:\n")
            if choice == 3:
                break
            else:
                continue


# Start:
callMenu()
