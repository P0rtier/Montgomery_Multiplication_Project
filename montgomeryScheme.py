# preparations
import time
import math


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


# Additional function returning coefficient m
def monProM(a, b, r, n, bit_width):
    t = a * b
    m = ((t % r) * neg_inv(n, bit_width, r)) % r
    return m


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


# Alternative function that calculates n'
# via equation represented in the article:
# p' = -p^-1 = p + 2 + x mod 2^(k+2)
# where x = 0 if k+2 <= 2e_2
# or x = 2^(2_e2)*3^(2_e3) mod 2^(k+2) if 2e_2 < k+2 <= 3e_2

def neg_invAlternate(n,k,e2,e3):
    x = 0
    if k+2 <= 2*e2:
        x = 0
    elif 2*e2 < k+2 <= 3*e2:
        x = (pow(2,2*e2) * pow(3,2*e3)) % pow(2,k+2)
    return n + 2 + x % pow(2,k+2)


# This function encodes the given y into the notation
# presented in the main task article, knowing that 2e_2
# bits of P' (n') and P + 2 are equal, we can
# define particular bits of elements as:
# Y_a ^ (B) = (Y mod B) >> a
# where >> is interpreted as a-sequence right shift (binary)

def give_newY(y,a,b):
    y1 = (y % pow(2,b))
    if(y1!=0):
        y1 >>= a
    return y1


# Main function of pipelined multiplication with coefficients,
# that are calculated via equations given in the main task article.
# Function returns result in the montgomery domain, hence for the appropriate
# interpretation, result must be further re-calculated to the standard domain.
# The coefficients of the multiplication are printed out on the terminal.

def alternativeModularMultiplication3(e2,e3,x,y):
    n = (pow(2,e2) * pow(3,e3)) - 1
    bit_width = calc_r(n)
    r = 2**calc_r(n)
    x = x * r % n
    y = y * r % n
    k = e2 + math.ceil(math.log2(pow(3,e3)))
    t = x*y
    n_prim = neg_invAlternate(n,k,e2,e3)
    #n_prim = neg_inv(n,bit_width,r)

    # coefficient C:
    t_temp = give_newY(t,e2,k+2-e2)
    n_prim_temp = give_newY(n_prim,e2,k+2-e2)
    C = t_temp * n_prim_temp

    # coefficient D:
    t_temp = give_newY(t,0,k+2-(2*e2))
    n_prim_temp = give_newY(n_prim,2*e2,k+2)
    D = t_temp * n_prim_temp

    # coefficient E:
    t_temp = give_newY(t,0,e2)
    n_prim_temp = give_newY(n_prim,e2,2*e2)
    E = t_temp * n_prim_temp
    m = (pow(2,2*e2)*((C+D) % pow(2,k+2-(2*e2))) + pow(2,e2) * give_newY(t,e2,k+2)
         + pow(2,e2)*E) % (pow(2,k+2) + give_newY(t,0,e2))
    m = monProM(x, y, r, n, bit_width)
    m = pow(2,e2) * give_newY(m,e2,k+2) + give_newY(t,0,e2)

    # coefficient F:
    m_temp = give_newY(m,e2,k+2)
    n_temp = give_newY(n,2*e2,k+2)
    n_prim_temp = give_newY(n_prim,e2,2*e2)
    F = m_temp * (pow(2,e2) * n_temp + n_prim_temp)

    # coefficient G:
    t_temp = give_newY(t,0,e2)
    n_temp = give_newY(n,2*e2,k+2)
    G = t_temp * n_temp

    MP = pow(2,2*e2) * F + pow(2,2*e2) * G + pow(2,e2) * E - m
    s = (t + MP) / r
    print("Coefficients of the multiplication:\n"
          "C: " + str(C) +
          "\nD: " + str(D) +
          "\nE: " + str(E) +
          "\nG: " + str(G) +
          "\nF: " + str(F) +
          "\nM: " + str(m))
    if s >= n:
        return s - n
    else:
        return s


def callMenu():
    while True:
        print("=====Montgomery Pipeline Multiplication and Exponentiation=====\n"
              "1. Multiplication\n"
              "2. Exponentiation\n"
              "3. Pipelined coefficient multiplication\n"
              "4. Exit\n"
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
                print("Calculator represents the x * y mod n\n"
                      "n = (2^e2)*(3^e3) - 1\n"
                      "Multiplication is executed via pipelined coefficients")
                try:
                    x = int(input("Insert x: "))
                    y = int(input("Insert y: "))
                    e2 = int(input("Insert e2: "))
                    e3 = int(input("Insert e3: "))
                except:
                    print("x,y,e2 and e3 must be Integers!!\n")
                else:
                    n = (pow(2,e2) * pow(3,e3)) - 1
                    bit_width = calc_r(n)
                    r = 2**bit_width
                    mid = alternativeModularMultiplication3(e2,e3,x,y)
                    print('Expected (x * y mod n) result:', x * y % n)
                    print('Montgomery Scheme result:', int(monPro(mid,1,r,n,bit_width)))
                    holdback = input("\nPress Enter to continue:\n")
            if choice == 4:
                break
            else:
                continue


# Start:
callMenu()
