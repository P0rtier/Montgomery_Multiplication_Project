# Montgomery Multiplication Project
  This project revolves around the implementation of the Fast Modular Pipelined multiplication,  
  firstly introduced by Peter Montgomery. The project consist of simple menu for user to import  
  their own factors of either direct multiplication or exponantiation.  

# Multiplication scheme
```python
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
    m = (t * neg_inv(n, bit_width, r)) % r
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
```
# Exponantiation scheme:
```python
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
```
