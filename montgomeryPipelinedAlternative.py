import math

def neg_inv(n, bit_width, base):
    def_base = 2 ** (bit_width - 1)
    neg_inverse = base - n ** (def_base - 1) % base
    return neg_inverse


def calc_r(n):
    for i in range(0, 32):
        r = 2**i
        if (r > n) and (r/2 <= n):
            return i


def monPro(a, b, r, n, bit_width):
    t = a * b
    m = ((t % r) * neg_inv(n, bit_width, r)) % r
    s = (a * b + m * n)/r
    if s >= n:
        return s - n
    else:
        return s

def monProM(a, b, r, n, bit_width):
    t = a * b
    m = ((t % r) * neg_inv(n, bit_width, r)) % r
    return m


def modMul(x, y, n):
    bit_width = calc_r(n)
    r = 2 ** bit_width
    a = x * r % n
    b = y * r % n
    montg_product = monPro(a, b, r, n,bit_width)
    product = monPro(montg_product, 1, r, n, bit_width)
    return product





def neg_invAlternate(n,k,e2,e3):
    x = 0
    if k+2 <= 2*e2:
        x = 0
    elif 2*e2 < k+2 <= 3*e2:
        x = (pow(2,2*e2) * pow(3,2*e3)) % pow(2,k+2)
    return n + 2 + x % pow(2,k+2)



def give_newY(y,a,b):
    y1 = (y % pow(2,b)) >> a
    return y1



# n = 2^e2 * 3^e3 - 1
# k = e2 + ceil(log_2(3^e3))

def alternativeModularMultiplication1(e2,e3,x,y):
    n = (pow(2,e2) * pow(3,e3)) - 1
    bit_width = calc_r(n)
    r = 2**bit_width
    k = e2 + math.ceil(math.log2(pow(3,e3)))
    t = x*y
    #n_prim = neg_invAlternate(n,k,e2,e3)
    n_prim = neg_inv(n,bit_width,r)
    #C:
    t = give_newY(t,e2,k+2-e2)
    n_prim = give_newY(n_prim,e2,k+2-e2)
    C = t * n_prim
    #D:
    t = give_newY(t,0,k+2-2*e2)
    n_prim = give_newY(n_prim,2*e2,k+2)
    D = t * n_prim


    first_tToM = give_newY(t,e2,k+2)
    #E:
    t = give_newY(first_tToM,0,e2)
    n_prim = give_newY(n_prim,e2,2*e2)
    E = t * n_prim

    #m1:
    t = give_newY(t,0,e2)
    m = (pow(2,e2)*((C+D) % pow(2,k+2-2*e2)) + pow(2,e2) * first_tToM+ pow(2,e2)*E) % pow(2,k+2) + t

    m = pow(2,e2) * give_newY(m,e2,k+2) + give_newY(t,0,e2)
    #F:
    m = give_newY(m,e2,k+2)
    n = give_newY(n,2*e2,k+2)
    n_prim = give_newY(n_prim,e2,2*e2)
    F = m * (pow(2,e2) * n + n_prim)
    #G:
    t = give_newY(t,0,e2)
    n = give_newY(n,2*e2,k+2)
    G = t * n
    MP = pow(2,2*e2) * F + pow(2,2*e2) * G + pow(2,e2) * E - m
    S = (t + MP)/r
    return S

def alternativeModularMultiplication2(e2,e3,x,y):
    n = (pow(2,e2) * pow(3,e3)) - 1
    bit_width = calc_r(n)
    r = 2**calc_r(n)
    x = x * r % n
    y = y * r % n
    k = e2 + math.ceil(math.log2(pow(3,e3)))
    t = x*y
    #n_prim = neg_invAlternate(n,k,e2,e3)
    n_prim = neg_inv(n,bit_width,r)
    #C:
    t_temp = give_newY(t,e2,k+2-e2)
    n_prim_temp = give_newY(n_prim,e2,k+2-e2)
    C = t_temp * n_prim_temp
    #D:
    t_temp = give_newY(t,0,k+2-2*e2)
    n_prim_temp = give_newY(n_prim,2*e2,k+2)
    D = t_temp * n_prim_temp


    first_tToM = give_newY(t,e2,k+2)
    #E:
    t_temp = give_newY(first_tToM,0,e2)
    n_prim_temp = give_newY(n_prim,e2,2*e2)
    E = t_temp * n_prim_temp

    #m1:
    t_temp = give_newY(t,0,e2)
    m = (pow(2,e2)*((C+D) % pow(2,k+2-2*e2)) + pow(2,e2) * first_tToM
         + pow(2,e2)*E) % pow(2,k+2) + t_temp

    m = pow(2,e2) * give_newY(m,e2,k+2) + give_newY(t,0,e2)
    #F:
    m_temp = give_newY(m,e2,k+2)
    n_temp = give_newY(n,2*e2,k+2)
    n_prim_temp = give_newY(n_prim,e2,2*e2)
    F = m_temp * (pow(2,e2) * n_temp + n_prim_temp)
    #G:
    t_temp = give_newY(t,0,e2)
    n_temp = give_newY(n,2*e2,k+2)
    G = t_temp * n_temp
    MP = pow(2,2*e2) * F + pow(2,2*e2) * G + pow(2,e2) * E - m
    m = monProM(x,y,r,n,bit_width)
    S = (t + m*n)/r
    return S

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
    #C:
    t_temp = give_newY(t,e2,k+2-e2)
    n_prim_temp = give_newY(n_prim,e2,k+2-e2)
    C = t_temp * n_prim_temp
    #D:
    t_temp = give_newY(t,0,k+2-(2*e2))
    n_prim_temp = give_newY(n_prim,2*e2,k+2)
    D = t_temp * n_prim_temp


    #E:
    t_temp = give_newY(t,0,e2)
    n_prim_temp = give_newY(n_prim,e2,2*e2)
    E = t_temp * n_prim_temp

    #m1:
    m = (pow(2,2*e2)*((C+D) % pow(2,k+2-(2*e2))) + pow(2,e2) * give_newY(t,e2,k+2)
         + pow(2,e2)*E) % pow(2,k+2) + give_newY(t,0,e2)

    m = monProM(x, y, r, n, bit_width)
    m = pow(2,e2) * give_newY(m,e2,k+2) + give_newY(t,0,e2)
    #F:
    m_temp = give_newY(m,e2,k+2)
    n_temp = give_newY(n,2*e2,k+2)
    n_prim_temp = give_newY(n_prim,e2,2*e2)
    F = m_temp * (pow(2,e2) * n_temp + n_prim_temp)
    #G:
    t_temp = give_newY(t,0,e2)
    n_temp = give_newY(n,2*e2,k+2)
    G = t_temp * n_temp
    MP = pow(2,2*e2) * F + pow(2,2*e2) * G + pow(2,e2) * E - m
    #m = monProM(x, y, r, n, bit_width)
    #m = pow(2,e2)*give_newY(m,e2,k+2) + give_newY(t,0,e2)
    s = (t + MP) / r
    if s >= n:
        return s - n
    else:
        return s

e2 = 3
e3 = 2
n = (pow(2,e2) * pow(3,e3)) - 1
x = 11
y = 12
xy = x*y
bit_width = calc_r(n)
r = 2 ** calc_r(n)
print('Oczekiwana:',xy % n)

print('Montg:',modMul(x,y,n))

mid_result = (alternativeModularMultiplication3(e2,e3,x,y))
print('Montg',int(monPro(mid_result, 1, r, n, bit_width)))


# Start:
#callMenu()
