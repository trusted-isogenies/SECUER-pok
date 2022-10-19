import random
from scipy import stats
from scipy import sparse, linalg


# This function comes from original IsogenyGraph repo
# (https://github.com/gfinol/IsogenyGraph/)
def elliptic_curve(F2):
    """ Returns a supersingular elliptic curve defined
    over F2, which can be used a starting node to
    generate the entire graph."""

    F = F2.subfield(1)
    R = F['s']
    (s,) = R._first_ngens(1)
    p = F2.characteristic()

    # Bröker's algorithm
    if p % 3 == 2:
        first_node = EllipticCurve(F2, [0, 1])
    elif p % 4 == 3:
        first_node = EllipticCurve(F2, [1, 0])
    else:
        q = 3
        while not (is_prime(q) and kronecker(-q, p) == -1):
            q += 4

        if q == 3:
            first_node = EllipticCurve(F2, [0, 1])
        else:
            K = QuadraticField(-q, names=('a',))
            (a,) = K._first_ngens(1)
            P_K = K.hilbert_class_polynomial()
            j = R(P_K).roots(multiplicities=False)[0]
            first_node = EllipticCurve(F2, j=j)

    return first_node

def equiv(j0, P0, j1, P1):
    """ Defines the equivalence relationship used
    to identify nodes in the graph with level structure.
    The functions assumes that if j0 == j1, then
    P0 and P1 must lie on the same curve."""

    if j0 != j1:
        return False

    assert P0.curve() == P1.curve()

    wp = P0.weil_pairing(P1, N)
    return wp == 1


# p is the characteristic of the field,
# while N represents the order of the subgroup
# used to define the level structure
# The code fixes the isogeny degree to be 2, but
# it can be easily generalized by modifying lines 88-89

p = 109
N = 3

print("p =", p)

adj = {}
visited = set()
to_visit = set()

# We're working over F_p^4 to get more torsion (as in Arpin's note)
# and because isomorphisms may only be defined over F_p^4
FF.<x> = GF(p)[]
F.<i> = GF(p^4)

j0 = elliptic_curve(F).j_invariant()
print("j0 =", j0)

E0 = EllipticCurve(F, j=j0)
P, Q = E0.gens()
# print(factor(P.order()))

assert P.order() % N == 0
P *= P.order() // N
Q *= P.order() // N



to_visit.add((j0, P))



while len(to_visit) != 0:
    (j, P) = to_visit.pop()
    visited.add((j, P))

    E = EllipticCurve(F, j=j)
    assert P.curve() == E

    nhbs = []

    two_torsion = E(0).division_points(2)
    two_torsion.remove(E(0))

    for K in two_torsion:
        phi = E.isogeny(K)
        Enew = phi.codomain()
        j_nhb = Enew.j_invariant()

        Eiso = EllipticCurve(F, j=j_nhb)
        iso = Enew.isomorphism_to(Eiso)

        pair = (j_nhb, iso(phi(P)))
        nhbs.append(pair)

        for v in visited.union(to_visit):
            if equiv(*pair, *v):
                break
        else:
            to_visit.add(pair)

    adj[(j, P)] = nhbs


visited = list(visited)

print("Curves:", len(set([x[0] for x in visited])))
print("Nodes:", len(visited))



# This may not be the most efficient way to do it,
# but it works and it "only" roughly doubles runtime

A = sparse.dok_matrix((len(visited), len(visited)))

for k, v in adj.items():
    x = -1
    for i, vv in enumerate(visited):
        if equiv(*vv, *k):
            x = i
            break

    for vv in v:
        y = -1
        for i, vvv in enumerate(visited):
            if equiv(*vvv, *vv):
                y = i
                break

        count = 0

        for vvv in v:
            count += equiv(*vv, *vvv)
        A[x, y] = count / len(v)


sparse.save_npz("%d-%d.npz" % (p, N), A.tocsc())

