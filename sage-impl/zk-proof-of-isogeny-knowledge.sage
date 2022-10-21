#!/usr/bin/env sage
proof.all(False)
warnings.simplefilter('ignore', FutureWarning)

####################################################################################
# ⚠ This implementation is out of date from the paper and the C implementation. ⚠ #
####################################################################################

# toy parameters
a, b = 8, 5
reps = 11
W, H = 2, 3

################################################################

p = 2^a * 3^b - 1
F.<i> = GF(p^2, modulus=[1,0,1])

# random point of order n on E
def random_point(E, n):
    assert n.divides(E.order())
    while (P := E.random_point()).order() % n:
        pass
    return P.order()//n * P

# isomorphism to deterministic representative
def canonical_model(E):
    EE = E.short_weierstrass_model()
    Es = {EE.change_weierstrass_model(u,0,0,0)
          for u in E.base_field().one().nth_root(6,all=True)}
    EE = min(Es, key=lambda C: C.a_invariants())
    return E.isomorphism_to(EE)

# isogeny from kernel point
@cached_function
def isogeny(K):
    E = K.curve()
    assert E == canonical_model(E).domain()
    phi = E.isogeny(K, algorithm='factored')
    return canonical_model(phi.codomain()) * phi

# given (kernel of top, kernel of left), return (kernel of bottom, kernel of right)
def complete_sidh_square(Kt, Kl):
    Kb = isogeny(Kl)(Kt)
    Kr = isogeny(Kt)(Kl)
    assert isogeny(Kb).codomain() == isogeny(Kr).codomain()
    return (Kb, Kr)

# pick a non-backtracking chain of L random degree-n isogenies
def generate(E, n, L):
    Ks = []
    back = None
    for _ in range(L):
        while True:
            K = random_point(E, n)
            if back is None or back(K).order() == n:
                break
        Ks.append(K)
        E = (psi := isogeny(K)).codomain()
        back = psi.dual()
    return Ks, E

# produce a single commitment for the interactive proof system
def commit(Ks_t):
    Ks_l, E2 = generate(Ks_t[0].curve(), 3^b, H)

    Ks_b = Ks_t[:]
    Ks_r = []
    for row in range(H):
        Kl = Ks_l[row]
        for col in range(W):
            Ks_b[col], Kl = complete_sidh_square(Ks_b[col], Kl)
        Ks_r.append(Kl)

    E3 = isogeny(Ks_b[-1]).codomain()

    # only reveal images of subgroups, not points
    Ks_b = [randrange(1,2^a,2)*K for K in Ks_b]

    res = {
            -1: (E3, Ks_l),     # note: each Ks[0] remembers its curve
             0: (E2, Ks_b),
            +1: (E2, Ks_r),
         }

    return res, (E2,E3)

# hash inputs & commitments to produce challenges (Fiat-Shamir)
def challenges(E0, E1, comms):
    transcript = (E0.a_invariants(), E1.a_invariants())

    # note: in the real implementation, we use a hiding commitment
    # here instead of simply revealing the curves!
    for E2,E3 in comms:
        transcript += (E2.a_invariants(), E3.a_invariants())

    from hashlib import shake_256
    num_byts = 40 + reps//5
    xof = shake_256(repr(transcript).encode())
    hash_ = int.from_bytes(xof.digest(num_byts), 'big')

    r = []
    for _ in range(reps):
        r.append(hash_ % 3 - 1)
        hash_ //= 3
    return tuple(r)

# the main function: random walk + proof of knowledge
def walk(E0):
    phi_Ks, E1 = generate(E0, 2^a, W)

    answs, comms = [], []
    for it in range(reps):
        print('  ' + 'o'*it + '·'*(reps-it), end='\r', file=sys.stderr)
        answ, comm = commit(phi_Ks)
        answs.append(answ)
        comms.append(comm)
    print('\r\x1b[K', file=sys.stderr)

    chall = challenges(E0, E1, comms)

    ress = []
    for idx,answ in zip(chall, answs):
        ress.append(answ[idx])
    proof = (chall, ress)

    return E1, proof

# the other main function: verifiying the proof
def verify(E0, E1, proof):
    chall, ress = proof

    if not len(chall) == len(ress) == reps:
        return False

    comms = []

    for it,(idx,answ) in enumerate(zip(chall, ress)):
        print('  ' + 'o'*it + '·'*(reps-it), end='\r', file=sys.stderr)

        Eo, Ks = answ

        E = Ks[0].curve()
        for K in Ks:
            assert K in E
            E = isogeny(K).codomain()

        if idx == -1:
            E2, E3 = E, Eo
        else:
            E2, E3 = Eo, E

        comms.append((E2, E3))
    print('\r\x1b[K', file=sys.stderr)

    return chall == challenges(E0, E1, comms)

################################################################

if __name__ == '__main__' and '__file__' in globals():

    print()

    E0 = EllipticCurve(F, [1,0])    # from previous contributor
    E0 = canonical_model(E0).codomain()

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

