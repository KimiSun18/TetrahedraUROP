import warnings
import numpy as np
from scipy.optimize import linprog
import math
from SymbolicTetra import *
f = open("output.txt", "w")

warnings.filterwarnings("ignore")  # CRTL has unecessary warnings due to LP


def SUB(l):  # outputs subsets of a set
    base = []
    lists = [base]
    for i in range(len(l)):
        orig = lists[:]
        new = l[i]
        for j in range(len(lists)):
            lists[j] = lists[j] + [new]
        lists = orig + lists
    return lists


P = SUB([0, 1, 2, 3, 4, 5])  # power sets of index set of edges

# simple fractions


def GCD(x, y):
    while(y):
        x, y = y, x % y
    return x


def VAL(R):
    return R[0]/R[1]


def RED(R):
    g = GCD(R[0], R[1])
    R[0] = int(R[0]/g)
    R[1] = int(R[1]/g)
    return R


def ADD(R, S):
    T = [R[0]*S[1]+R[1]*S[0], R[1]*S[1]]
    return RED(T)


def SUB(R, S):
    T = [R[0]*S[1]-R[1]*S[0], R[1]*S[1]]
    return RED(T)


def MUL(R, S):
    T = [R[0]*S[0], R[1]*S[1]]
    return RED(T)


def DIV(R, S):
    T = [R[0]*S[1], R[1]*S[0]]
    return RED(T)


def DTP(R):
    T = DIV([2, 1], R)
    return T[1] == 1


def GEN():  # returns the pi or 2pi combinations for a general member of family 2
    M = []
    for a in range(13):
        for b in range(13):
            for c in range(13):
                for d in range(13):
                    if b+d == a+c and (5*a+b+4*c == 6 or 5*a+b+4*c == 12):
                        M.append([a, b, c, d, 5*a+b+4*c])
    return M


def SPE():  # returns all specific members of family 2
    A = [5, 6]
    B = [1, 6]
    C = [2, 3]
    L = []
    for a in range(13):
        for b in range(13):
            for c in range(13):
                for d in range(13):
                    if b+d != a+c:
                        x = DIV(SUB([1, 1], ADD(ADD(MUL([a, 1], A), MUL(
                            [b, 1], B)), MUL([c, 1], C))), [b+d-a-c, 1])
                        y = DIV(SUB([2, 1], ADD(ADD(MUL([a, 1], A), MUL(
                            [b, 1], B)), MUL([c, 1], C))), [b+d-a-c, 1])
                        if VAL(x) > 1/6 and VAL(x) <= 1/3:
                            l = [SUB(A, x), ADD(B, x), SUB(
                                C, x), SUB(C, x), x, x]
                            if l not in L:
                                L.append(l)
                        if VAL(y) > 1/6 and VAL(y) <= 1/3:
                            l = [SUB(A, y), ADD(B, y), SUB(
                                C, y), SUB(C, y), y, y]
                            if l not in L:
                                L.append(l)
    return L


def CON(S):  # converts list of tetrahedra in KKPR form to a list of those in fraction form
    T = []
    for i in range(len(S)):
        t = []
        for j in range(6):
            t.append(RED([S[i][j+1], S[i][0]]))
        T.append(t)
    for i in range(len(T)):  # assuming KKPR switches 13 and 24
        T[i][2], T[i][3] = T[i][3], T[i][2]
    return T


def DIH(T):  # puts tetra in form as in Anas' program
    return [VAL(T[0])*math.pi, VAL(T[2])*math.pi, VAL(T[4])*math.pi, VAL(T[5])*math.pi, VAL(T[3])*math.pi, VAL(T[1])*math.pi]


def EDG(T):  # returns edgelengths removing duplicates
    D = calculate_lengths(DIH(T))
    for i in range(4):
        for j in range(4):
            if i > j:
                D.pop((j+1, i+1))
    # print(D)
    return D


def MR(l, x):  # max range to check for combinations of l for an angle of size x
    return math.floor(l/x)+1


def PID(T):  # returns pi combinations of dihedral angles
    L = []
    r0 = 1.0001
    for a in range(MR(r0, VAL(T[0]))):
        r1 = r0-a*VAL(T[0])
        for b in range(MR(r1, VAL(T[1]))):
            r2 = r1-b*VAL(T[1])
            for c in range(MR(r2, VAL(T[2]))):
                r3 = r2-c*VAL(T[2])
                for d in range(MR(r3, VAL(T[3]))):
                    r4 = r3-d*VAL(T[3])
                    for e in range(MR(r4, VAL(T[4]))):
                        r5 = r4-e*VAL(T[4])
                        for f in range(MR(r5, VAL(T[5]))):
                            r = r5-f*VAL(T[5])
                            if abs(r-0.0001) < 0.0001:
                                L.append([a, b, c, d, e, f])
    return L


def ADC(T):  # returns and prints pi and 2pi combinations of dihedral angles
    L = []
    r0 = 2.0001
    for a in range(MR(r0, VAL(T[0]))):
        r1 = r0-a*VAL(T[0])
        for b in range(MR(r1, VAL(T[1]))):
            r2 = r1-b*VAL(T[1])
            for c in range(MR(r2, VAL(T[2]))):
                r3 = r2-c*VAL(T[2])
                for d in range(MR(r3, VAL(T[3]))):
                    r4 = r3-d*VAL(T[3])
                    for e in range(MR(r4, VAL(T[4]))):
                        r5 = r4-e*VAL(T[4])
                        for f in range(MR(r5, VAL(T[5]))):
                            r = r5-f*VAL(T[5])
                            if abs(r-1.0001) < 0.0001 or abs(r-0.0001) < 0.0001:
                                L.append([a, b, c, d, e, f])
    return L


def EUD(L):  # Given dihedral angle combinations, returns the number of times each edge is used in KKPR order
    D = [0, 0, 0, 0, 0, 0]
    for i in range(6):
        for l in L:
            if l[i] > 0:
                D[i] += 1
    return D


def CRTN12(T):  # true if no pi combination or only one and it involes only a pair of opposite edges
    L = PID(T)
    C = [0, 0, 0, 0, 0, 0]
    for l in L:
        for i in range(6):
            if l[i] > 0 and C[i] == 0:
                C[i] = 1
    return (C == [0, 0, 0, 0, 0, 0] or C == [1, 1, 0, 0, 0, 0] or C == [0, 0, 1, 1, 0, 0] or C == [0, 0, 0, 0, 1, 1])


def PIP(T):  # Returns pi combinations of planar angles
    L = []
    P = overall(DIH(T))
    r0 = 2*math.pi+0.0001
    for a in range(MR(r0, P[(2, 1, 3)])):
        r1 = r0 - a*P[(2, 1, 3)]
        for b in range(MR(r1, P[(2, 1, 4)])):
            r2 = r1-b*P[(2, 1, 4)]
            for c in range(MR(r2, P[(3, 1, 4)])):
                r3 = r2-c*P[(3, 1, 4)]
                for d in range(MR(r3, P[(1, 2, 3)])):
                    r4 = r3-d*P[(1, 2, 3)]
                    for e in range(MR(r4, P[(1, 2, 4)])):
                        r5 = r4-e*P[(1, 2, 4)]
                        for f in range(MR(r5, P[(3, 2, 4)])):
                            r6 = r5-f*P[(3, 2, 4)]
                            for g in range(MR(r6, P[(1, 3, 2)])):
                                r7 = r6-g*P[(1, 3, 2)]
                                for h in range(MR(r7, P[(1, 3, 4)])):
                                    r8 = r7-h*P[(1, 3, 4)]
                                    for i in range(MR(r8, P[(2, 3, 4)])):
                                        r9 = r8-i*P[(2, 3, 4)]
                                        for j in range(MR(r9, P[(1, 4, 2)])):
                                            r10 = r9-j*P[(1, 4, 2)]
                                            for k in range(MR(r10, P[(1, 4, 3)])):
                                                r11 = r10-k*P[(1, 4, 3)]
                                                for l in range(MR(r11, P[(2, 4, 3)])):
                                                    r = r11-l*P[(2, 4, 3)]
                                                    if abs(r-0.0001) < 0.0001 or abs(r-math.pi-0.0001) < 0.0001:
                                                        L.append(
                                                            [a, b, c, d, e, f, g, h, i, j, k, l])
    return L


def EUP(L):  # Given planar angle combinations, returns the number of times each edge is used in KKPR order
    D = [L[0]+L[1]+L[3]+L[4], L[7]+L[8]+L[9]+L[10], L[0]+L[2]+L[6]+L[7],
         L[4]+L[5]+L[9]+L[11], L[1]+L[9]+L[10]+L[11], L[3]+L[5]+L[6]+L[8]]
    return D


def CRTN3(T):  # true if every way planar angles make pi or 2pi involves some angle in no pi combination
    D = EUD(PID(T))
    P = PIP(T)
    Q = []
    for i in range(6):
        if D[i] == 0:
            for j in range(len(P)):
                if EUP(P[j])[i] > 0 and j not in Q:
                    Q.append(j)
    return(len(Q) == len(P))


def CRTF0(T):  # true if some dihedral angle does not divide 2pi
    for t in T:
        if DIV([2, 1], t)[1] != 1:
            return True
    return False


def CRTF1(T):  # true if all edges distinct
    D = EDG(T)
    for key1, value1 in D.items():
        for key2, value2 in D.items():
            if key1 != key2 and abs(value1 - value2) < 0.0001:
                return False
    return True


def CRTF2(T):  # true if in "non symmetric parallelogram form", e.g. as those in family 2
    D = calculate_lengths(DIH(T))
    if abs(D[(1, 2)]-D[(3, 4)]) < 0.0001 and abs(D[(2, 3)]-D[(1, 4)]) < 0.0001 and abs(D[(1, 3)]-D[(2, 4)]) > 0.0001 and abs(D[(1, 2)]-D[(2, 3)]) > 0.0001:
        return True
    elif abs(D[(1, 2)]-D[(3, 4)]) < 0.0001 and abs(D[(2, 3)]-D[(1, 4)]) > 0.0001 and abs(D[(1, 3)]-D[(2, 4)]) < 0.0001 and abs(D[(1, 2)]-D[(1, 3)]) > 0.0001:
        return True
    elif abs(D[(1, 2)]-D[(3, 4)]) > 0.0001 and abs(D[(2, 3)]-D[(1, 4)]) < 0.0001 and abs(D[(1, 3)]-D[(2, 4)]) < 0.0001 and abs(D[(1, 3)]-D[(2, 3)]) > 0.0001:
        return True
    else:
        return False


def PRI(L):  # prints out list L of tetrahedra in latex form to file output.txt
    for j in range(len(L)):
        T = L[j]
        print("\\left(", file=f)
        for i in range(len(T)):
            print("\\frac{", T[i][0], "}{", T[i][1], "}", file=f)
            if i < len(T)-1:
                print(", ", file=f)
        print("\\right), ", file=f)
        if j % 4 == 3:
            print("\\\\ ", file=f)


S = [  # the 1 from family 2 + all 59 sporadic tetrahedra
    [12, 3, 4, 3, 4, 6, 8],
    [24, 5, 9, 6, 8, 13, 15],
    [12, 3, 6, 4, 6, 4, 6],
    [24, 7, 11, 7, 13, 8, 12],
    [15, 3, 3, 3, 5, 10, 10],
    [15, 2, 4, 4, 4, 10, 10],
    [15, 3, 3, 4, 4, 9, 11],
    [15, 3, 3, 5, 5, 9, 9],
    [15, 5, 5, 5, 9, 6, 6],
    [15, 3, 7, 6, 6, 7, 7],
    [15, 4, 8, 5, 5, 7, 7],
    [21, 3, 9, 7, 7, 12, 12],
    [21, 4, 10, 6, 6, 12, 12],
    [21, 6, 6, 7, 7, 9, 15],
    [30, 6, 12, 10, 15, 10, 20],
    [30, 4, 14, 10, 15, 12, 18],
    [60, 8, 28, 19, 31, 25, 35],
    [60, 12, 24, 15, 35, 25, 35],
    [60, 13, 23, 15, 35, 24, 36],
    [60, 13, 23, 19, 31, 20, 40],
    [30, 6, 18, 10, 10, 15, 15],
    [30, 4, 16, 12, 12, 15, 15],
    [30, 9, 21, 10, 10, 12, 12],
    [30, 6, 6, 10, 12, 15, 20],
    [30, 5, 7, 11, 11, 15, 20],
    [60, 7, 17, 20, 24, 35, 35],
    [60, 7, 17, 22, 22, 33, 37],
    [60, 10, 14, 17, 27, 35, 35],
    [60, 12, 12, 17, 27, 33, 37],
    [30, 6, 10, 10, 15, 12, 18],
    [30, 5, 11, 10, 15, 13, 17],
    [60, 10, 22, 21, 29, 25, 35],
    [60, 11, 21, 19, 31, 26, 34],
    [60, 11, 21, 21, 29, 24, 36],
    [60, 12, 20, 19, 31, 25, 35],
    [30, 6, 10, 6, 10, 15, 24],
    [60, 7, 25, 12, 20, 35, 43],
    [30, 6, 12, 6, 12, 15, 20],
    [60, 12, 24, 13, 23, 29, 41],
    [30, 6, 12, 10, 10, 15, 18],
    [30, 7, 13, 9, 9, 15, 18],
    [60, 12, 24, 17, 23, 33, 33],
    [60, 14, 26, 15, 21, 33, 33],
    [60, 15, 21, 20, 20, 27, 39],
    [60, 17, 23, 18, 18, 27, 39],
    [30, 6, 15, 6, 18, 10, 20],
    [30, 6, 15, 7, 17, 9, 21],
    [60, 9, 33, 14, 34, 21, 39],
    [60, 9, 33, 15, 33, 20, 40],
    [60, 11, 31, 12, 36, 21, 39],
    [60, 11, 31, 15, 33, 18, 42],
    [30, 6, 15, 10, 15, 12, 15],
    [30, 6, 15, 11, 14, 11, 16],
    [30, 8, 13, 8, 17, 12, 15],
    [30, 8, 13, 9, 18, 11, 14],
    [30, 8, 17, 9, 12, 11, 16],
    [30, 9, 12, 9, 18, 10, 15],
    [30, 10, 12, 10, 12, 15, 12],
    [60, 19, 25, 20, 24, 29, 25]
]


def MXD(T):  # compute the 6 by s D-matrix D_i (e) of T
    E = calculate_lengths(DIH(T))
    C = ADC(T)
    D = []
    for c in C:
        d = [0, 0, 0, 0, 0, 0]
        d[0] = c[0]/E[(1, 2)]
        d[1] = c[1]/E[(3, 4)]
        d[2] = c[2]/E[(1, 3)]
        d[3] = c[3]/E[(2, 4)]
        d[4] = c[4]/E[(1, 4)]
        d[5] = c[5]/E[(2, 3)]
        D.append(d)
    return D

# weak linear programming criterion

# def CRTL1(T):  # true if criteria D holds and prints out the index sets if exists
#     D = MXD(T)
#     for I in P:
#         for J in P:
#             if len(I) > 0 and len(J) > 0:
#                 B = True
#                 for d in D:
#                     Si = 0
#                     for i in I:
#                         Si += d[i]
#                     Ai = Si / len(I)
#                     Sj = 0
#                     for j in J:
#                         Sj += d[j]
#                     Aj = Sj / len(J)
#                     if Ai <= Aj:
#                         B = False
#                 if B:
#                     return (True, I, J)
#     return (False,)

# strong linear programming criterion


def CRTL2(T):  # true if the dual problem has a solution
    D = MXD(T)
    A = np.array(D)*-1
    c = []
    for i in range(6):
        c.append([1])  # c = (1, ... ,1)^t in R^6
    b = []
    for i in range(len(D)):
        b.append([0])  # b = 0 vecor in R^s
    # minimize c^tx over x such that Ax<=0 with no bounds
    l = linprog(c, A, b, None, None, (None, None))
    return (l.fun < -0.001, l)

# # The one from the second family and the 42 sporadics left undecided
# L = [
#     # [[7, 12], [5, 12], [5, 12], [5, 12], [1, 4], [1, 4]], # FROM SECOND FAMILY, DONE by LP
#     [[1, 4], [1, 3], [1, 3], [1, 4], [1, 2], [2, 3]],
#     [[5, 24], [3, 8], [1, 3], [1, 4], [13, 24], [5, 8]],
#     [[1, 4], [1, 2], [1, 2], [1, 3], [1, 3], [1, 2]],
#     [[7, 24], [11, 24], [13, 24], [7, 24], [1, 3], [1, 2]],
#     [[1, 5], [1, 5], [1, 3], [1, 5], [2, 3], [2, 3]],
#     [[1, 5], [1, 5], [4, 15], [4, 15], [3, 5], [11, 15]],
#     [[1, 5], [1, 5], [1, 3], [1, 3], [3, 5], [3, 5]],
#     [[1, 3], [1, 3], [3, 5], [1, 3], [2, 5], [2, 5]],
#     [[4, 15], [8, 15], [1, 3], [1, 3], [7, 15], [7, 15]],
#     [[1, 7], [3, 7], [1, 3], [1, 3], [4, 7], [4, 7]],
#     [[2, 7], [2, 7], [1, 3], [1, 3], [3, 7], [5, 7]],
#     [[1, 5], [2, 5], [1, 2], [1, 3], [1, 3], [2, 3]],
#     [[2, 15], [7, 15], [1, 2], [1, 3], [2, 5], [3, 5]],
#     [[2, 15], [7, 15], [31, 60], [19, 60], [5, 12], [7, 12]],
#     [[1, 5], [2, 5], [7, 12], [1, 4], [5, 12], [7, 12]],
#     [[13, 60], [23, 60], [7, 12], [1, 4], [2, 5], [3, 5]],
#     [[1, 5], [3, 5], [1, 3], [1, 3], [1, 2], [1, 2]],
#     [[3, 10], [7, 10], [1, 3], [1, 3], [2, 5], [2, 5]],
#     [[1, 5], [1, 5], [2, 5], [1, 3], [1, 2], [2, 3]],
#     [[1, 6], [7, 30], [11, 30], [11, 30], [1, 2], [2, 3]],
#     [[1, 5], [1, 5], [9, 20], [17, 60], [11, 20], [37, 60]],
#     [[1, 5], [1, 3], [1, 2], [1, 3], [2, 5], [3, 5]],
#     [[1, 6], [11, 30], [1, 2], [1, 3], [13, 30], [17, 30]],
#     [[1, 6], [11, 30], [29, 60], [7, 20], [5, 12], [7, 12]],
#     [[1, 5], [1, 3], [1, 3], [1, 5], [1, 2], [4, 5]],
#     [[7, 60], [5, 12], [1, 3], [1, 5], [7, 12], [43, 60]],
#     [[1, 5], [2, 5], [2, 5], [1, 5], [1, 2], [2, 3]],
#     [[1, 5], [2, 5], [1, 3], [1, 3], [1, 2], [3, 5]],
#     [[7, 30], [13, 30], [3, 10], [3, 10], [1, 2], [3, 5]],
#     [[1, 4], [7, 20], [1, 3], [1, 3], [9, 20], [13, 20]],
#     [[1, 5], [1, 2], [3, 5], [1, 5], [1, 3], [2, 3]],
#     [[1, 5], [1, 2], [17, 30], [7, 30], [3, 10], [7, 10]],
#     [[3, 20], [11, 20], [17, 30], [7, 30], [7, 20], [13, 20]],
#     # [[3, 20], [11, 20], [11, 20], [1, 4], [1, 3], [2, 3]] , # DONE by LP
#     # [[11, 60], [31, 60], [11, 20], [1, 4], [3, 10], [7, 10]] , # DONE by LP
#     [[1, 5], [1, 2], [1, 2], [1, 3], [2, 5], [1, 2]],
#     [[1, 5], [1, 2], [7, 15], [11, 30], [11, 30], [8, 15]],
#     [[4, 15], [13, 30], [17, 30], [4, 15], [2, 5], [1, 2]],
#     [[4, 15], [13, 30], [3, 5], [3, 10], [11, 30], [7, 15]],
#     [[4, 15], [17, 30], [2, 5], [3, 10], [11, 30], [8, 15]],
#     [[3, 10], [2, 5], [3, 5], [3, 10], [1, 3], [1, 2]],
#     [[1, 3], [2, 5], [2, 5], [1, 3], [1, 2], [2, 5]]]


def EXC(T):
    P = overall(DIH(T))
    E1 = [P[(2, 1, 3)], P[(2, 1, 4)], P[(3, 1, 4)]]
    E2 = [P[(1, 2, 3)], P[(1, 2, 4)], P[(3, 2, 4)]]
    E3 = [P[(1, 3, 4)], P[(1, 3, 4)], P[(2, 3, 4)]]
    E4 = [P[(1, 4, 2)], P[(1, 4, 3)], P[(2, 4, 3)]]
    A = [E1, E2, E3, E4]
    EX = []
    for E in A:
        s = 0
        for e in E:
            s += e/2
        p = 1
        for i in range(3):
            p *= np.tan(float((s-E[i])/2))
        EX.append(4*np.arctan(np.sqrt(p)))
    return EX


# for T in CON(S):
#     print('EXC', EXC(T))
#     d = [math.pi*VAL(t) for t in T]
#     e = [d[0]+d[2]+d[4]-math.pi, d[0]+d[3]+d[5]-math.pi,
#          d[1]+d[2]+d[5]-math.pi, d[1]+d[3]+d[4]-math.pi]
#     print('New', e)
#     print()


def COR(T):
    # L = EXC(T)+[2*math.pi*VAL(t) for t in T]
    d = [math.pi*VAL(t) for t in T]
    e = [d[0]+d[2]+d[4]-math.pi, d[0]+d[3]+d[5]-math.pi,
         d[1]+d[2]+d[5]-math.pi, d[1]+d[3]+d[4]-math.pi]
    L = e+[2*math.pi*VAL(t) for t in T]
    print(L)
    R1 = []
    R2 = []
    r0 = 4*math.pi+0.0001
    for a0 in range(MR(r0, L[0])):
        r1 = r0-a0*L[0]
        for a1 in range(MR(r1, L[1])):
            r2 = r1-a1*L[1]
            for a2 in range(MR(r2, L[2])):
                r3 = r2-a2*L[2]
                for a3 in range(MR(r3, L[3])):
                    r4 = r3-a3*L[3]
                    for a4 in range(MR(r4, L[4])):
                        r5 = r4-a4*L[4]
                        for a5 in range(MR(r5, L[5])):
                            r6 = r5-a5*L[5]
                            for a6 in range(MR(r6, L[6])):
                                r7 = r6-a6*L[6]
                                for a7 in range(MR(r7, L[7])):
                                    r8 = r7-a7*L[7]
                                    for a8 in range(MR(r8, L[8])):
                                        r9 = r8-a8*L[8]
                                        for a9 in range(MR(r9, L[9])):
                                            r = r9 - a9 * L[9]
                                            if abs(r-0.0001) < 0.0001:
                                                R2.append(
                                                    [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9])
                                            if abs(r-0.0001-2*math.pi) < 0.0001:
                                                R1.append(
                                                    [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9])
    return (R1, R2)


def CRTE(T):
    a = False
    b = False
    R1, R2 = COR(T)
    I = [0, 0, 0, 0]
    for r in R1:
        for i in range(4):
            if r[i] > 0:
                I[i] += 1
    f0 = CRTF0(T)
    f1 = CRTF1(T)
    f2 = CRTF2(T)
    no_f2f = f0 and (f1 or f2)
    if I == [0, 0, 0, 0] and no_f2f:
        a = True
    J = [0, 0, 0, 0]
    for r in R1+R2:
        for i in range(4):
            if r[i] > 0:
                J[i] += 1
    if J[0] == 0 or J[1] == 0 or J[2] == 0 or J[3] == 0:
        b = True
    return (a, b)


def PNT(T):  # can be proven to not tile
    L = CRTL2(T)[0]
    n12 = CRTN12(T)
    n3 = CRTN3(T)
    N = n12 or n3
    f0 = CRTF0(T)
    f1 = CRTF1(T)
    f2 = CRTF2(T)
    F = f0 and (f1 or f2)
    e = CRTE(T)
    E = (e[0] and F) or e[1]
    #print('L', L)
    # print('N', n12 or n3)
    # print('F', F)
    # print('     e[0]', e[0])
    # print('     e[1]', e[1])
    # print('E', E)
    return ((N and F) or L or E, (N and F), L, E)


a = [0, 0, 0, 0]
for T in CON(S):
    res = PNT(T)
    print(res[0], T)
    print(res[1:])
    for i in range(4):
        if res[i]:
            a[i] += 1
    print()
print(a[0], "out of 59 sporadics untilable")
print('Due to [OG, L, E]', a[1:])
print()


a = [0, 0, 0, 0]
for T in SPE():
    res = PNT(T)
    print(res[0], T)
    print(res[1:])
    for i in range(4):
        if res[i]:
            a[i] += 1
    print()
print(a[0], "out of 23 specifics untilable")
print('Due to [OG, L, E]', a[1:])
print()
