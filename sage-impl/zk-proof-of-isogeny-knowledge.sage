#!/usr/bin/env sage

# Written by Lorenz Panny

proof.all(False)

W, H = 2, 4             # glueing W × H SIDH squares together

iterations = 9          #XXX >=438, maybe round up to 448 or 512?
                        # ...actually, do we really need 2^256 security here?

e2, e3 = 8, 5
#e2, e3 = 43, 26
#e2, e3 = 216, 137   # SIKEp434

################################################################

p = 2^e2 * 3^e3 - 1
F.<ii> = GF(p^2, modulus=polygen(ZZ)^2+1)

# requires Sage ≥ 9.5
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

# random point of order n on E
def rpt(E,n):
    assert n.divides(E.order())
    while (P := E.random_point()).order() % n: pass
    return P.order()//n * P

# isogeny from kernel point
@cached_function
def isog(K):
    return EllipticCurveHom_composite(K.curve(), K)

# given (kernel of top, kernel of left), return (kernel of bottom, kernel of right)
def complete_sidh_square(Kt, Kl):
    Kb = isog(Kl)(Kt)
    Kr = isog(Kt)(Kl)
    assert isog(Kb).codomain() == isog(Kr).codomain()
    return (Kb, Kr)

# pick a chain of L random degree-n isogenies
def generate(E, n, L):
    Ks = []
    for i in range(L):
        Ks.append(K := rpt(E, n))
        E = isog(K).codomain()
    return Ks, E

# produce a single commitment for the interactive proof system
def commit(Ks_t):
    Ks_l, E2 = generate(Ks_t[0].curve(), 3^e3, H)

    Ks_b = Ks_t[:]
    Ks_r = []
    for row in range(H):
        Kl = Ks_l[row]
        for col in range(W):
            Ks_b[col], Kl = complete_sidh_square(Ks_b[col], Kl)
        Ks_r.append(Kl)

    E3 = isog(Ks_b[-1]).codomain()

    # challenge 0: reveal left
    # challenge 1: reveal bottom
    # challenge 2: reveal right
    resp = ((E3,Ks_l), (E2,Ks_b), (E2,Ks_r))

    return resp, (E2,E3)

# hash inputs & commitments to produce challenges (Fiat-Shamir)
def challenges(E0, E1, comms):
    transcript = (E0.a_invariants(), E1.a_invariants())
    for E2,E3 in comms:
        transcript += (E2.a_invariants(), E3.a_invariants())

    from hashlib import shake_256
    num_byts = 40 + iterations//5
    xof = shake_256(repr(transcript).encode())
    hash_ = int.from_bytes(xof.digest(num_byts), 'big')

    r = []
    for _ in range(iterations):
        r.append(hash_ % 3)
        hash_ //= 3
    return tuple(r)

# the main function: random walk + proof of knowledge
def walk(E0):
    phi_Ks, E1 = generate(E0, 2^e2, W)

    answs, comms = [], []
    for it in range(iterations):
        print('  ' + 'o'*it + '·'*(iterations-it), end='\r', file=sys.stderr)
        answ, comm = commit(phi_Ks)
        answs.append(answ)
        comms.append(comm)
    print('\r\x1b[K', file=sys.stderr)

    chall = challenges(E0, E1, comms)

    resps = []
    for idx,answ in zip(chall, answs):
        resps.append(answ[idx])
    proof = (chall, resps)

    return E1, proof

# the other main function: verifiying the proof
def verify(E0, E1, proof):
    chall, resps = proof

    comms = []

    for idx,answ in zip(chall, resps):

        Eo, Ks = answ

        E = Ks[0].curve()
        for K in Ks:
            assert K in E
            E = isog(K).codomain()

        if idx == 0:
            E2, E3 = E, Eo
        else:
            E2, E3 = Eo, E

        comms.append((E2, E3))

    return chall == challenges(E0, E1, comms)

################################################################

if __name__ == '__main__' and '__file__' in globals():

    print()

    E0 = EllipticCurve(F, [1,0])    # from previous contributor

    print('input curve:')
    print(f'    {E0}')
    print()

    E1, proof = walk(E0)

    print('output curve:')
    print(f'    {E1}')
    print()

    print('proof:')
    print('    ' + '\n    '.join(map(str, zip(*proof))))
    print()

    assert verify(E0, E1, proof)

