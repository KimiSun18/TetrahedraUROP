import numpy as np
import math


def minor(M, i, j):
    return [row[:j] + row[j+1:] for row in (M[:i]+M[i+1:])]


def val(n, p):
    if n == 0:
        return(float("inf"), 0)
    c = 0
    while n % p == 0:
        c += 1
        n = n / p
    return (int(c), int(n))


def primes(n):
    n = abs(n)
    S = set()
    d = 2
    while d*d <= n:
        if n % d == 0:
            S.add(d)
            while (n % d) == 0:
                n //= d
        d += 1
    if n > 1:
        S.add(n)
    return S


def st(i, j, k=None):
    if k is not None:
        return tuple(sorted((i, j, k)))
    else:
        return tuple(sorted((i, j)))


def tetra(e12, e34, e13, e24, e14, e23):
    e = {
        (1, 2): e12,
        (3, 4): e34,
        (1, 3): e13,
        (2, 4): e24,
        (1, 4): e14,
        (2, 3): e23,
    }
    E12 = pow(e12, 2)
    E34 = pow(e34, 2)
    E13 = pow(e13, 2)
    E24 = pow(e24, 2)
    E14 = pow(e14, 2)
    E23 = pow(e23, 2)
    f = {
        (1, 2): E12,
        (3, 4): E34,
        (1, 3): E13,
        (2, 4): E24,
        (1, 4): E14,
        (2, 3): E23,
    }
    M = [[0, E12, E13, E14, 1], [E12, 0, E23, E24, 1], [
        E13, E23, 0, E34, 1], [E14, E24, E34, 0, 1], [1, 1, 1, 1, 0]]
    # Cayley-Menger
    D = {}  # Determinants and minors
    D[()] = int(round(np.linalg.det(np.array(M))))
    for (i, j, k) in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)):
        D[(i, j, k)] = -(e[st(i, j)]+e[st(i, k)]+e[st(j, k)])*(-e[st(i, j)]+e[st(i, k)] +
                                                               e[st(j, k)])*(e[st(i, j)]-e[st(i, k)]+e[st(j, k)])*(e[st(i, j)]+e[st(i, k)]-e[st(j, k)])
    for (i, j) in ((1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3)):
        L = [1, 2, 3, 4]
        L.remove(i)
        L.remove(j)
        k = L[0]
        l = L[1]
        D[(i, j)] = -1*pow(f[st(i, j)], 2)+(f[st(i, k)]+f[st(i, l)]+f[st(j, k)]+f[st(j, l)
                                                                                  ]-2*f[st(k, l)])*f[st(i, j)]+(f[st(i, k)]-f[st(j, k)])*(f[st(j, l)]-f[st(i, l)])
        # print(e[(i, j)], '* arccos(',
        #       D[(i, j)], '/sqrt(', D[st(i, j, k)], '*', D[st(i, j, l)], '))')
    return e, D


def padic(e12, e34, e13, e24, e14, e23):
    e, D = tetra(e12, e34, e13, e24, e14, e23)
    Valid = True
    for (i, j, k) in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)):
        if 2*max(e[(i, j)], e[(j, k)], e[(i, k)]) >= e[(i, j)]+e[(j, k)]+e[(i, k)]:
            Valid = False
    if not Valid or D[()] <= 0:
        return ('N', (e12, e34, e13, e24, e14, e23))
    P = set()  # list of prime divisors of Dijk
    for l in range(4):
        L = [1, 2, 3, 4]
        L.remove(l+1)
        P = P | primes(D[tuple(L)])
    # print('primes', P)
    # P = {17}
    # print(e)
    # print(D)

    def helper(p, e, D):
        # print()
        # print(p)
        V = {}
        T = 0
        for (i, j) in ((1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3)):
            L = [1, 2, 3, 4]
            L.remove(i)
            L.remove(j)
            k = L[0]
            l = L[1]
            t = val(D[st(i, j, k)], p)[0] + val(D[st(i, j, l)], p)[0]
            if t > T:
                T = t
            V[(i, j)] = val(D[(i, j)], p)[0] - t / 2
            if p == 2:
                V[(i, j)] += 1
        # COMPUTE SQRT(-2D) up to T digits
        r, u = val(-2*D[()], p)
        if r % 2 != 0:
            # print('power invalid')
            return (False, V)
        if p == 2:
            if u % 8 != 1:
                # print('2 invalid')
                return (False, V)
            else:
                x = 1
        else:
            x = 0
            for y in range(p):
                if (pow(y, 2)-u) % p == 0:
                    x = y
                    break
            if x == 0:
                # print(p, 'invalid')
                return (False, V)

        for r in range(1, T):
            for y in range(p):
                if (pow(x+y*pow(p, r), 2) + 2*D[()]) % pow(p, r+1) == 0:
                    x += y*pow(p, r)
                    break

        # print(p, x)

        for (i, j) in ((1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3)):
            L = [1, 2, 3, 4]
            L.remove(i)
            L.remove(j)
            k = L[0]
            l = L[1]
            if V[(i, j)] >= 0:
                V[(i, j)] = 0
                # print('zero', V[(i, j)])
            else:
                t = val(D[st(i, j, k)], p)[0] + val(D[st(i, j, l)], p)[0]
                if (2*pow(D[(i, j)], 2)+2*e[(i, j)]*D[(i, j)]*x) % pow(p, t) == 0:
                    sgn = 1
                else:
                    sgn = -1
                V[(i, j)] = sgn * abs(V[(i, j)])
            # print(p)
        print(p, V)

        res = 0  # compute linear relation
        for (i, j) in ((1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3)):
            res += e[(i, j)]*V[(i, j)]
        return (res != 0, V)

    for p in P:
        res = helper(p, e, D)
        if res[0]:
            # print(p, (e12, e34, e13, e24, e14, e23))
            # print(res[1])
            return (p, (e12, e34, e13, e24, e14, e23))
    print((e12, e34, e13, e24, e14, e23))
    # print(D)
    return ('F', (e12, e34, e13, e24, e14, e23))


def FourParam(l):
    f = 0
    for a in range(1, l+1):
        for b in range(1, a+1):
            for c in range(1, l+1):
                for d in range(1, l+1):
                    # ADD: TO REDUCE CASES BY SYMMETRIES
                    # print(e12, e34, e13, e24, e14, e23)
                    g = 1  # math.gcd(math.gcd(a, b), math.gcd(c, d))
                    if g == 1:
                        res = padic(a, b, a, b, c, d)
                        if res[0] == 'F':
                            f += 1
    return f


def ThreeParam(l):  # Turns out to be Hill's first family :(
    for a in range(1, l+1):
        for c in range(a, l+1):
            for b in range(1, l+1):
                # ADD: TO REDUCE CASES BY SYMMETRIES
                # print(e12, e34, e13, e24, e14, e23)
                g = math.gcd(math.gcd(a, b), c)
                if g == 1:
                    res = padic(a, a, a, b, c, c)
                    if res[0] == 'F':
                        print((a, a, a, b, c, c), b*b, 3*(c-a)*(c+a))
    return


def CheckCube(l):
    F = []
    for e12 in range(1, l+1):
        for e34 in range(1, e12+1):  # WLOG
            for e13 in range(1, l+1):
                for e24 in range(1, e13+1):  # WLOG
                    for e14 in range(1, l+1):
                        for e23 in range(1, l+1):
                            g = math.gcd(math.gcd(math.gcd(e12, e34), math.gcd(
                                e13, e24)), math.gcd(e14, e23))
                            isF = (e12 == e13 and e24 == e34) or (e12 == e14 and e23 == e34) or (e13 == e14 and e23 == e24) or (
                                e12 == e24 and e13 == e34) or (e12 == e23 and e14 == e34) or (e13 == e23 and e14 == e24)
                            s = sorted([e12, e34, e13, e24, e14, e23])
                            isH = (s[1] == s[2] and s[2] == s[3] and s[4] == s[5]) or (
                                s[0] == s[1] and s[1] == s[2] and s[4] == s[5]) or (s[0] == s[1] and s[1] == s[2] and s[3] == s[4])
                            if g == 1 and (not isF) and (not isH):
                                res = padic(e12, e34, e13, e24, e14, e23)
                                # if res[0] == 'F':
                                #     print(res)
    return


def tetra2(l):
    print(l)
    e12, e34, e13, e24, e14, e23 = l
    e = {
        (1, 2): e12,
        (3, 4): e34,
        (1, 3): e13,
        (2, 4): e24,
        (1, 4): e14,
        (2, 3): e23,
    }
    E12 = pow(e12, 2)
    E34 = pow(e34, 2)
    E13 = pow(e13, 2)
    E24 = pow(e24, 2)
    E14 = pow(e14, 2)
    E23 = pow(e23, 2)
    f = {
        (1, 2): E12,
        (3, 4): E34,
        (1, 3): E13,
        (2, 4): E24,
        (1, 4): E14,
        (2, 3): E23,
    }
    M = [[0, E12, E13, E14, 1], [E12, 0, E23, E24, 1], [
        E13, E23, 0, E34, 1], [E14, E24, E34, 0, 1], [1, 1, 1, 1, 0]]
    # Cayley-Menger
    D = {}  # Determinants and minors
    D[()] = int(round(np.linalg.det(np.array(M))))
    for (i, j, k) in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)):
        D[(i, j, k)] = -(e[st(i, j)]+e[st(i, k)]+e[st(j, k)])*(-e[st(i, j)]+e[st(i, k)] +
                                                               e[st(j, k)])*(e[st(i, j)]-e[st(i, k)]+e[st(j, k)])*(e[st(i, j)]+e[st(i, k)]-e[st(j, k)])
    for (i, j) in ((1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3)):
        L = [1, 2, 3, 4]
        L.remove(i)
        L.remove(j)
        k = L[0]
        l = L[1]
        D[(i, j)] = -1*pow(f[st(i, j)], 2)+(f[st(i, k)]+f[st(i, l)]+f[st(j, k)]+f[st(j, l)
                                                                                  ]-2*f[st(k, l)])*f[st(i, j)]+(f[st(i, k)]-f[st(j, k)])*(f[st(j, l)]-f[st(i, l)])
        # print(e[(i, j)], '* arccos(',
        #       D[(i, j)], '/sqrt(', D[st(i, j, k)], '*', D[st(i, j, l)], '))')
    # print(D)
    Valid = True
    for (i, j, k) in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)):
        if 2*max(e[(i, j)], e[(j, k)], e[(i, k)]) >= e[(i, j)]+e[(j, k)]+e[(i, k)]:
            Valid = False
    if not Valid or D[()] <= 0:
        return 0.5
    A = {}
    B = {}
    C = {}
    E = {}
    for (i, j) in ((1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3)):
        L = [1, 2, 3, 4]
        L.remove(i)
        L.remove(j)
        k = L[0]
        l = L[1]
        C[(i, j)] = D[(i, j)] / math.pow(D[st(i, j, k)]*D[st(i, j, l)], 0.5)
        A[(i, j)] = math.acos(C[(i, j)])
        B[(i, j)] = A[(i, j)]/math.pi
        E[(i, j)] = (D[(i, j)], D[st(i, j, k)]*D[st(i, j, l)])
    # print(C)
    print(B)
    # print(E)
    s = 0
    for (i, j) in ((1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3)):
        s += e[(i, j)]*A[(i, j)]
    print('res:', s, s/math.pi)
    print()
    return s/math.pi

# for i in tetra(13, 9, 13, 9, 11, 12):
#     print(i)
# padic(17, 15, 17, 15, 16, 6)


def isF(e):
    (e12, e34, e13, e24, e14, e23) = e
    return (e12 == e13 and e24 == e34) or (e12 == e14 and e23 == e34) or (e13 == e14 and e23 == e24) or (
        e12 == e24 and e13 == e34) or (e12 == e23 and e14 == e34) or (e13 == e23 and e14 == e24)


def isH(e):
    (e12, e34, e13, e24, e14, e23) = e
    s = sorted([e12, e34, e13, e24, e14, e23])
    return (s[1] == s[2] and s[2] == s[3] and s[4] == s[5]) or (
        s[0] == s[1] and s[1] == s[2] and s[4] == s[5]) or (s[0] == s[1] and s[1] == s[2] and s[3] == s[4])


def t1(e):
    (e12, e34, e13, e24, e14, e23) = e
    s = (e13+e24+e14+e23)/2
    return (e12, e34, s-e24, s-e13, s-e23, s-e14)


def t2(e):
    (e12, e34, e13, e24, e14, e23) = e
    s = (e12+e34+e14+e23)/2
    return (s-e34, s-e12, e13, e24, s-e23, s-e14)


def t3(e):
    (e12, e34, e13, e24, e14, e23) = e
    s = (e12+e34+e13+e24)/2
    return (s-e34, s-e12, s-e24, s-e13, e14, e23)


def allsym(t):
    e = {}
    f = {}
    e[(1, 2)], e[(3, 4)], e[(1, 3)], e[(2, 4)], e[(1, 4)], e[(2, 3)] = t
    L = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 5):
                    if i != j and i != k and i != l and j != k and j != l and k != l:
                        f[st(i, j)], f[st(k, l)], f[st(i, k)], f[st(j, l)], f[st(i, l)], f[st(
                            j, k)] = e[(1, 2)], e[(3, 4)], e[(1, 3)], e[(2, 4)], e[(1, 4)], e[(2, 3)]
                        L.append((f[(1, 2)], f[(3, 4)], f[(1, 3)],
                                  f[(2, 4)], f[(1, 4)], f[(2, 3)]))
    return L


def CheckCube2(l):
    F = []
    for e12 in range(1, l+1):
        for e34 in range(1, e12+1):  # WLOG
            for e13 in range(1, l+1):
                for e24 in range(1, e13+1):  # WLOG
                    for e14 in range(1, l+1):
                        for e23 in range(1, l+1):
                            g = math.gcd(math.gcd(math.gcd(e12, e34), math.gcd(
                                e13, e24)), math.gcd(e14, e23))
                            if g == 1 and (not isF((e12, e34, e13, e24, e14, e23))) and (not isH((e12, e34, e13, e24, e14, e23))):
                                res = tetra2((e12, e34, e13, e24, e14, e23))
                                if abs(res - round(res)) < 0.000001:
                                    print((e12, e34, e13, e24, e14, e23))
                                # if res[0] == 'F':
                                #     print(res)
    return


# for i in tetra(7, 5, 7, 5, 6, 6):
#     print(i)
# print(padic(7, 5, 7, 5, 6, 6))
# tetra2((7, 5, 7, 5, 6, 6))

# print(CheckCube2(30))

# Not in FourParam or Hill First

# L = [(10, 8, 13, 11, 15, 9), (11, 9, 12, 8, 9, 13), (11, 11, 13, 13,
#                                                      12, 16), (13, 12, 16, 8, 14, 16), (17, 15, 19, 13, 14, 20), (19, 11, 13, 11, 12, 13)]

# X = [(9, 7, 11, 9, 12, 8),
#      (9, 7, 12, 8, 11, 9),
#      (10, 8, 13, 11, 15, 9),
#      (10, 8, 15, 9, 13, 11),
#      (11, 9, 9, 7, 12, 8),
#      (11, 9, 12, 8, 9, 7),
#      (11, 9, 12, 8, 9, 13),
#      (11, 9, 13, 9, 8, 12),
#      (11, 9, 15, 13, 18, 10),
#      (11, 9, 18, 10, 15, 13),
#      (11, 9, 26, 22, 30, 16),
#      (11, 9, 30, 16, 26, 22),
#      (11, 10, 24, 22, 30, 14),
#      (11, 10, 30, 14, 24, 22),
#      (11, 11, 13, 13, 12, 16),
#      (11, 11, 13, 13, 16, 12),
#      (11, 11, 16, 12, 13, 13),
#      (12, 8, 9, 7, 11, 9),
#      (12, 8, 11, 9, 9, 7),
#      (12, 8, 11, 9, 9, 13),
#      (12, 8, 13, 9, 9, 11),
#      (12, 8, 25, 23, 30, 16),
#      (12, 8, 30, 16, 25, 23),
#      (12, 10, 17, 15, 21, 11),
#      (12, 10, 21, 11, 17, 15),
#      (13, 7, 22, 21, 23, 24),
#      (13, 7, 24, 23, 21, 22),
#      (13, 9, 11, 9, 8, 12),
#      (13, 9, 12, 8, 9, 11),
#      (13, 11, 10, 8, 15, 9),
#      (13, 11, 13, 12, 11, 19),
#      (13, 11, 15, 9, 10, 8),
#      (13, 11, 19, 11, 12, 13),
#      (13, 11, 19, 17, 24, 12),
#      (13, 11, 21, 17, 27, 12),
#      (13, 11, 24, 12, 19, 17),
#      (13, 11, 27, 12, 21, 17),
#      (13, 12, 13, 11, 11, 19),
#      (13, 12, 16, 8, 14, 16),
#      (13, 12, 16, 14, 8, 16),
#      (13, 12, 19, 11, 11, 13),
#      (13, 13, 11, 11, 12, 16),
#      (13, 13, 11, 11, 16, 12),
#      (13, 13, 16, 12, 11, 11),
#      (14, 10, 20, 18, 27, 12),
#      (14, 10, 27, 12, 20, 18),
#      (14, 12, 21, 19, 27, 13),
#      (14, 12, 27, 13, 21, 19),
#      (15, 9, 10, 8, 13, 11),
#      (15, 9, 13, 11, 10, 8),
#      (15, 11, 15, 13, 16, 12),
#      (15, 11, 16, 12, 15, 13),
#      (15, 13, 11, 9, 18, 10),
#      (15, 13, 15, 11, 16, 12),
#      (15, 13, 16, 12, 15, 11),
#      (15, 13, 18, 10, 11, 9),
#      (15, 13, 18, 10, 19, 9),
#      (15, 13, 19, 9, 18, 10),
#      (15, 13, 23, 21, 30, 14),
#      (15, 13, 30, 14, 23, 21),
#      (16, 4, 27, 25, 30, 24),
#      (16, 4, 28, 26, 29, 23),
#      (16, 4, 29, 23, 28, 26),
#      (16, 4, 30, 24, 27, 25),
#      (16, 8, 13, 12, 14, 16),
#      (16, 8, 16, 14, 12, 13),
#      (16, 12, 11, 11, 13, 13),
#      (16, 12, 13, 13, 11, 11),
#      (16, 12, 15, 11, 15, 13),
#      (16, 12, 15, 13, 15, 11),
#      (16, 12, 17, 11, 13, 19),
#      (16, 12, 19, 13, 11, 17),
#      (16, 14, 13, 12, 8, 16),
#      (16, 14, 16, 8, 12, 13),
#      (16, 16, 29, 25, 16, 30),
#      (16, 16, 29, 25, 30, 16),
#      (16, 16, 30, 16, 25, 29),
#      (16, 16, 30, 16, 29, 25),
#      (17, 3, 24, 22, 26, 22),
#      (17, 3, 25, 21, 25, 23),
#      (17, 3, 25, 23, 25, 21),
#      (17, 3, 26, 22, 24, 22),
#      (17, 11, 16, 12, 13, 19),
#      (17, 11, 19, 13, 12, 16),
#      (17, 13, 20, 16, 21, 15),
#      (17, 13, 21, 15, 20, 16),
#      (17, 15, 12, 10, 21, 11),
#      (17, 15, 19, 13, 14, 20),
#      (17, 15, 20, 14, 13, 19),
#      (17, 15, 21, 11, 12, 10),
#      (18, 10, 11, 9, 15, 13),
#      (18, 10, 15, 13, 11, 9),
#      (18, 10, 15, 13, 19, 9),
#      (18, 10, 19, 9, 15, 13),
#      (18, 14, 27, 27, 16, 30),
#      (18, 14, 27, 27, 30, 16),
#      (18, 14, 30, 16, 27, 27),
#      (19, 9, 15, 13, 18, 10),
#      (19, 9, 18, 10, 15, 13),
#      (19, 11, 13, 11, 12, 13),
#      (19, 11, 13, 12, 11, 13),
#      (19, 13, 16, 12, 11, 17),
#      (19, 13, 17, 11, 12, 16),
#      (19, 13, 17, 15, 14, 20),
#      (19, 13, 20, 14, 15, 17),
#      (19, 15, 24, 20, 27, 17),
#      (19, 15, 27, 17, 24, 20),
#      (19, 17, 13, 11, 24, 12),
#      (19, 17, 21, 15, 7, 29),
#      (19, 17, 24, 12, 13, 11),
#      (19, 17, 29, 7, 15, 21),
#      (19, 19, 27, 22, 14, 30),
#      (19, 19, 27, 22, 30, 14),
#      (19, 19, 30, 14, 22, 27),
#      (19, 19, 30, 14, 27, 22),
#      (20, 7, 27, 25, 28, 29),
#      (20, 7, 29, 28, 25, 27),
#      (20, 11, 21, 11, 30, 9),
#      (20, 11, 30, 9, 21, 11),
#      (20, 14, 17, 15, 13, 19),
#      (20, 14, 19, 13, 15, 17),
#      (20, 15, 23, 19, 17, 29),
#      (20, 15, 25, 21, 15, 27),
#      (20, 15, 27, 15, 21, 25),
#      (20, 15, 29, 17, 19, 23),
#      (20, 16, 17, 13, 21, 15),
#      (20, 16, 21, 15, 17, 13),
#      (20, 18, 14, 10, 27, 12),
#      (20, 18, 27, 12, 14, 10),
#      (21, 11, 12, 10, 17, 15),
#      (21, 11, 17, 15, 12, 10),
#      (21, 11, 20, 11, 30, 9),
#      (21, 11, 24, 15, 26, 5),
#      (21, 11, 26, 5, 24, 15),
#      (21, 11, 30, 9, 20, 11),
#      (21, 15, 17, 13, 20, 16),
#      (21, 15, 19, 17, 7, 29),
#      (21, 15, 20, 16, 17, 13),
#      (21, 15, 22, 14, 17, 25),
#      (21, 15, 25, 17, 14, 22),
#      (21, 15, 26, 18, 30, 26),
#      (21, 15, 29, 7, 17, 19),
#      (21, 15, 30, 26, 26, 18),
#      (21, 17, 13, 11, 27, 12),
#      (21, 17, 27, 12, 13, 11),
#      (21, 19, 14, 12, 27, 13),
#      (21, 19, 22, 16, 23, 17),
#      (21, 19, 23, 17, 22, 16),
#      (21, 19, 26, 14, 26, 29),
#      (21, 19, 27, 13, 14, 12),
#      (21, 19, 29, 26, 14, 26),
#      (22, 14, 21, 15, 17, 25),
#      (22, 14, 25, 17, 15, 21),
#      (22, 14, 25, 19, 30, 26),
#      (22, 14, 27, 14, 17, 27),
#      (22, 14, 27, 17, 14, 27),
#      (22, 14, 30, 26, 25, 19),
#      (22, 16, 21, 19, 23, 17),
#      (22, 16, 23, 17, 21, 19),
#      (22, 21, 13, 7, 23, 24),
#      (22, 21, 24, 23, 7, 13),
#      (22, 22, 23, 23, 25, 27),
#      (22, 22, 23, 23, 27, 25),
#      (22, 22, 27, 25, 23, 23),
#      (23, 9, 23, 23, 25, 29),
#      (23, 9, 23, 23, 29, 25),
#      (23, 9, 25, 21, 27, 27),
#      (23, 9, 27, 27, 21, 25),
#      (23, 9, 27, 27, 25, 21),
#      (23, 9, 29, 25, 23, 23),
#      (23, 13, 26, 18, 14, 27),
#      (23, 13, 27, 14, 18, 26),
#      (23, 17, 21, 19, 22, 16),
#      (23, 17, 22, 16, 21, 19),
#      (23, 17, 24, 20, 25, 19),
#      (23, 17, 25, 19, 24, 20),
#      (23, 19, 20, 15, 17, 29),
#      (23, 19, 29, 17, 15, 20),
#      (23, 21, 15, 13, 30, 14),
#      (23, 21, 26, 18, 19, 27),
#      (23, 21, 27, 19, 18, 26),
#      (23, 21, 30, 13, 29, 30),
#      (23, 21, 30, 14, 15, 13),
#      (23, 21, 30, 29, 13, 30),
#      (23, 23, 22, 22, 25, 27),
#      (23, 23, 22, 22, 27, 25),
#      (23, 23, 23, 9, 25, 29),
#      (23, 23, 23, 9, 29, 25),
#      (23, 23, 27, 25, 22, 22),
#      (23, 23, 29, 25, 9, 23),
#      (23, 23, 29, 25, 23, 9),
#      (24, 12, 13, 11, 19, 17),
#      (24, 12, 19, 17, 13, 11),
#      (24, 15, 21, 11, 26, 5),
#      (24, 15, 26, 5, 21, 11),
#      (24, 20, 19, 15, 27, 17),
#      (24, 20, 23, 17, 25, 19),
#      (24, 20, 25, 19, 23, 17),
#      (24, 20, 27, 17, 19, 15),
#      (24, 22, 11, 10, 30, 14),
#      (24, 22, 17, 3, 26, 22),
#      (24, 22, 26, 22, 17, 3),
#      (24, 22, 30, 14, 11, 10),
#      (24, 23, 13, 7, 21, 22),
#      (24, 23, 22, 21, 7, 13),
#      (25, 17, 21, 15, 14, 22),
#      (25, 17, 22, 14, 15, 21),
#      (25, 19, 22, 14, 30, 26),
#      (25, 19, 23, 17, 24, 20),
#      (25, 19, 24, 20, 23, 17),
#      (25, 19, 29, 23, 30, 22),
#      (25, 19, 30, 22, 29, 23),
#      (25, 19, 30, 26, 22, 14),
#      (25, 21, 17, 3, 25, 23),
#      (25, 21, 20, 15, 15, 27),
#      (25, 21, 23, 9, 27, 27),
#      (25, 21, 25, 23, 17, 3),
#      (25, 21, 27, 15, 15, 20),
#      (25, 21, 27, 27, 9, 23),
#      (25, 21, 27, 27, 23, 9),
#      (25, 22, 29, 9, 30, 27),
#      (25, 22, 30, 27, 29, 9),
#      (25, 23, 12, 8, 30, 16),
#      (25, 23, 17, 3, 25, 21),
#      (25, 23, 25, 21, 17, 3),
#      (25, 23, 30, 16, 12, 8),
#      (26, 5, 21, 11, 24, 15),
#      (26, 5, 24, 15, 21, 11),
#      (26, 14, 21, 19, 26, 29),
#      (26, 14, 29, 26, 19, 21),
#      (26, 18, 21, 15, 30, 26),
#      (26, 18, 23, 13, 14, 27),
#      (26, 18, 23, 21, 19, 27),
#      (26, 18, 27, 14, 13, 23),
#      (26, 18, 27, 19, 21, 23),
#      (26, 18, 30, 26, 21, 15),
#      (26, 22, 11, 9, 30, 16),
#      (26, 22, 17, 3, 24, 22),
#      (26, 22, 24, 22, 17, 3),
#      (26, 22, 30, 16, 11, 9),
#      (27, 12, 13, 11, 21, 17),
#      (27, 12, 14, 10, 20, 18),
#      (27, 12, 20, 18, 14, 10),
#      (27, 12, 21, 17, 13, 11),
#      (27, 13, 14, 12, 21, 19),
#      (27, 13, 21, 19, 14, 12),
#      (27, 14, 22, 14, 17, 27),
#      (27, 14, 23, 13, 18, 26),
#      (27, 14, 26, 18, 13, 23),
#      (27, 14, 27, 17, 14, 22),
#      (27, 15, 20, 15, 21, 25),
#      (27, 15, 25, 21, 15, 20),
#      (27, 17, 19, 15, 24, 20),
#      (27, 17, 22, 14, 14, 27),
#      (27, 17, 24, 20, 19, 15),
#      (27, 17, 27, 14, 14, 22),
#      (27, 19, 23, 21, 18, 26),
#      (27, 19, 26, 18, 21, 23),
#      (27, 22, 19, 19, 14, 30),
#      (27, 22, 19, 19, 30, 14),
#      (27, 22, 30, 14, 19, 19),
#      (27, 25, 16, 4, 30, 24),
#      (27, 25, 20, 7, 28, 29),
#      (27, 25, 22, 22, 23, 23),
#      (27, 25, 23, 23, 22, 22),
#      (27, 25, 29, 21, 30, 22),
#      (27, 25, 29, 28, 7, 20),
#      (27, 25, 30, 22, 29, 21),
#      (27, 25, 30, 24, 16, 4),
#      (27, 27, 18, 14, 16, 30),
#      (27, 27, 18, 14, 30, 16),
#      (27, 27, 23, 9, 21, 25),
#      (27, 27, 23, 9, 25, 21),
#      (27, 27, 25, 21, 9, 23),
#      (27, 27, 25, 21, 23, 9),
#      (27, 27, 30, 16, 14, 18),
#      (27, 27, 30, 16, 18, 14),
#      (28, 26, 16, 4, 29, 23),
#      (28, 26, 29, 23, 16, 4),
#      (29, 7, 19, 17, 15, 21),
#      (29, 7, 21, 15, 17, 19),
#      (29, 9, 25, 22, 30, 27),
#      (29, 9, 30, 27, 25, 22),
#      (29, 17, 20, 15, 19, 23),
#      (29, 17, 23, 19, 15, 20),
#      (29, 21, 27, 25, 30, 22),
#      (29, 21, 30, 22, 27, 25),
#      (29, 23, 16, 4, 28, 26),
#      (29, 23, 25, 19, 30, 22),
#      (29, 23, 28, 26, 16, 4),
#      (29, 23, 30, 22, 25, 19),
#      (29, 25, 16, 16, 16, 30),
#      (29, 25, 16, 16, 30, 16),
#      (29, 25, 23, 9, 23, 23),
#      (29, 25, 23, 23, 9, 23),
#      (29, 25, 23, 23, 23, 9),
#      (29, 25, 30, 16, 16, 16),
#      (29, 26, 21, 19, 14, 26),
#      (29, 26, 26, 14, 19, 21),
#      (29, 28, 20, 7, 25, 27),
#      (29, 28, 27, 25, 7, 20),
#      (30, 9, 20, 11, 21, 11),
#      (30, 9, 21, 11, 20, 11),
#      (30, 13, 23, 21, 29, 30),
#      (30, 13, 30, 29, 21, 23),
#      (30, 14, 11, 10, 24, 22),
#      (30, 14, 15, 13, 23, 21),
#      (30, 14, 19, 19, 22, 27),
#      (30, 14, 19, 19, 27, 22),
#      (30, 14, 23, 21, 15, 13),
#      (30, 14, 24, 22, 11, 10),
#      (30, 14, 27, 22, 19, 19),
#      (30, 16, 11, 9, 26, 22),
#      (30, 16, 12, 8, 25, 23),
#      (30, 16, 16, 16, 25, 29),
#      (30, 16, 16, 16, 29, 25),
#      (30, 16, 18, 14, 27, 27),
#      (30, 16, 25, 23, 12, 8),
#      (30, 16, 26, 22, 11, 9),
#      (30, 16, 27, 27, 14, 18),
#      (30, 16, 27, 27, 18, 14),
#      (30, 16, 29, 25, 16, 16),
#      (30, 22, 25, 19, 29, 23),
#      (30, 22, 27, 25, 29, 21),
#      (30, 22, 29, 21, 27, 25),
#      (30, 22, 29, 23, 25, 19),
#      (30, 24, 16, 4, 27, 25),
#      (30, 24, 27, 25, 16, 4),
#      (30, 26, 21, 15, 26, 18),
#      (30, 26, 22, 14, 25, 19),
#      (30, 26, 25, 19, 22, 14),
#      (30, 26, 26, 18, 21, 15),
#      (30, 27, 25, 22, 29, 9),
#      (30, 27, 29, 9, 25, 22),
#      (30, 29, 23, 21, 13, 30),
#      (30, 29, 30, 13, 21, 23)]
# Y = []
# for x in X:
#     # print(x)
#     b = isF(x) or isF(t1(x)) or isF(t2(x)) or isF(t3(x)) or isF(t1(t2(x))) or isF(t2(t3(x))) or isF(t3(t1(x))) or isH(
#         x) or isH(t1(x)) or isH(t2(x)) or isH(t3(x)) or isH(t1(t2(x))) or isH(t2(t3(x))) or isH(t3(t1(x)))
#     S = allsym(x)
#     for s in S:
#         if s in Y:
#             b = True
#     if not b:
#         Y.append(x)
# for y in Y:
#     print(y)  # , st(y[0]+y[1], y[2]+y[3], y[4]+y[5]))

Y = [(11, 9, 26, 22, 30, 16),
     (11, 10, 24, 22, 30, 14),
     (12, 8, 25, 23, 30, 16),
     (13, 7, 22, 21, 23, 24),
     (13, 11, 21, 17, 27, 12),
     (13, 12, 16, 8, 14, 16),
     (14, 10, 20, 18, 27, 12),
     (15, 13, 18, 10, 19, 9),
     (16, 4, 27, 25, 30, 24),
     (16, 4, 28, 26, 29, 23),
     (16, 16, 29, 25, 16, 30),
     (17, 3, 24, 22, 26, 22),
     (17, 3, 25, 21, 25, 23),
     (18, 14, 27, 27, 16, 30),
     (19, 17, 21, 15, 7, 29),
     (20, 7, 27, 25, 28, 29),
     (20, 11, 21, 11, 30, 9),
     (20, 15, 23, 19, 17, 29),
     (20, 15, 25, 21, 15, 27),
     (21, 11, 24, 15, 26, 5),
     (21, 15, 26, 18, 30, 26),
     (21, 19, 26, 14, 26, 29),
     (22, 14, 25, 19, 30, 26),
     (22, 14, 27, 14, 17, 27),
     (23, 9, 23, 23, 25, 29),
     (23, 9, 25, 21, 27, 27),
     (23, 13, 26, 18, 14, 27),
     (23, 21, 30, 13, 29, 30),
     (25, 22, 29, 9, 30, 27),
     (13, 11, 13, 12, 11, 19)]  # H2

for y in Y:
    tetra2(y)

# L = []
# S = []

# for x in X:
#     s = sorted(x)
#     if s not in S:
#         S.append(s)
#         L.append(x)

# # for l in L:
# #     print(l)

# for l in L:
#     tetra2(l)

# x = (11, 9, 26, 22, 30, 16)
# print(allsym(x))
# print(t1(x))
# print(t2(x))
# print(t3(x))
# print(t1(t2(x)))
# print(t2(t1(x)))
# print(t2(t3(x)))
# print(t3(t2(x)))
# print(t1(t3(x)))
# print(t3(t1(x)))
# print(t1(t2(t3(x))))
