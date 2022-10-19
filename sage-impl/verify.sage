#!/usr/bin/env sage
proof.all(False)
import sage.schemes.elliptic_curves.hom_composite

#TODO input validation

################################################################

#TODO should store these things in the proof file?

p = 2^216 * 3^137 - 1           # SIKEp434  #TODO other sizes

reps = 219
width = 4
height = 7

################################################################

F.<i> = GF((p,2), modulus=polygen(ZZ)^2+1)

def read_bytes():
    return bytes.fromhex(input())

def read_element():
    re,im = input().split(',')
    return int(re,0) + int(im,0)*i

def invalid():
    print('\x1b[31minvalid proof\x1b[0m', file=sys.stderr)
    exit(1)

E0 = read_element()
E1 = read_element()
tags = []
response = []
while len(response) < reps:
    tag = int(input())
    tags.append(tag)
    if tag == -1:               # left
        r2 = read_bytes()
        h3 = read_bytes()
        chain = [read_element() for _ in range(height)]
        response.append((tag, r2, h3, chain))
    elif tag == 0:              # bottom
        r2 = read_bytes()
        r3 = read_bytes()
        E2 = read_element()
        chain = [read_element() for _ in range(width)]
        response.append((tag, r2, r3, E2, chain))
    elif tag == +1:             # right
        h2 = read_bytes()
        r3 = read_bytes()
        chain = [read_element() for _ in range(height)]
        response.append((tag, h2, r3, chain))
    else:
        invalid()

################################################################

def pack(el):
    assert el in F
    plen = (p.bit_length() + 7) // 8
    return b''.join(int(v).to_bytes(plen,'little') for v in el)

from Crypto.Hash import SHAKE256

def hash_hiding(r, v):
    h = SHAKE256.new()
    h.update(b'\x02')       # RO_COMMIT
    h.update(r)
    h.update(pack(v))
    return h.read(int(32))  #XXX 32 is HASH_BYTES from poik.h

def curve(A):
    E = EllipticCurve([0,A,0,1,0])
    E.set_order((p+1)^2)    # is a SIKE curve
    return E

def challenge(hs):
    h = SHAKE256.new()
    h.update(b'\x01')       # RO_CHALLENGE
    h.update(pack(curve(E0).j_invariant()))
    h.update(pack(curve(E1).j_invariant()))
    for h2,h3 in hs:
        h.update(h2 + h3)

    trits = []
    while len(trits) < reps:
        v, = h.read(int(1))
        if v >= 243:
            continue
        for j in range(5):
            v,r = divmod(v, 3)
            trits.append({2:-1}.get(r,r))
    return trits[:reps]

def canonical(E):
    R.<X> = FunctionField(F)
    eq = E.j_invariant() - EllipticCurve([0,X,0,1,0]).j_invariant()
    As = eq.numerator().roots(multiplicities=False)
    A = min(As, key = lambda el: (int(el[0]), int(el[1])))
    return E.isomorphism_to(curve(A))

def unpack_chain(A, chain):
    E = curve(A)
    for i,x in enumerate(chain):
        K = E.lift_x(x)
        E = E.isogeny(K, algorithm='factored').codomain()
        E = canonical(E).codomain()
    return E

def unpack_response(resp):
    if resp[0] == -1:
        r2, h3, chain = resp[1:]
        j2 = unpack_chain(E0, chain).j_invariant()
        h2 = hash_hiding(r2, j2)
    elif resp[0] == 0:
        r2, r3, E2, chain = resp[1:]
        j2 = curve(E2).j_invariant()
        j3 = unpack_chain(E2, chain).j_invariant()
        h2 = hash_hiding(r2, j2)
        h3 = hash_hiding(r3, j3)
    elif resp[0] == +1:
        h2, r3, chain = resp[1:]
        j3 = unpack_chain(E1, chain).j_invariant()
        h3 = hash_hiding(r3, j3)
    else:
        invalid()
    return h2, h3

################################################################

@parallel(os.cpu_count())
def work(idx):
    return unpack_response(response[idx])

hs = [None]*len(response)
for ((idx,),_),(h2,h3) in work(list(range(len(response)))):
    print(f'{idx:3} {h2.hex()} {h3.hex()}', file=sys.stderr)
    hs[idx] = (h2,h3)

if challenge(hs) == tags:
    print(f'\x1b[32msuccess!\x1b[0m', file=sys.stderr)
    print(E0)
    print(E1)
    exit()

invalid()

