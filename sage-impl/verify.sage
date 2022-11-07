#!/usr/bin/env sage
proof.all(False)
warnings.simplefilter('ignore', FutureWarning)

import sys, re
from typing import NamedTuple
from Crypto.Hash import SHAKE256

################################################################

def invalid():
    print('\x1b[31minvalid proof\x1b[0m', file=sys.stderr)
    exit(1)

def read_int():
    s = input()
    if not re.match('^-?[0-9]+$', s):
        invalid()
    return int(s)

def read_bytes():
    s = input()
    if not re.match('^([0-9a-f]{2})+$', s, re.IGNORECASE):
        invalid()
    return bytes.fromhex(s)

################################################################

def pack(el):
    p = parent(el).characteristic()
    plen = (p.bit_length() + 7) // 8
    return b''.join(int(v).to_bytes(plen,'little') for v in el)

def hash_hiding(r, v):
    h = SHAKE256.new()
    h.update(b'\x02')       # RO_COMMIT
    h.update(r)
    h.update(pack(v))
    return h.read(int(32))  #XXX 32 is HASH_BYTES from poik.h

def curve(A):
    E = EllipticCurve([0,A,0,1,0])
    p = parent(A).characteristic()
    E.set_order((p+1)^2)    # is a SIKE curve
    return E

def challenge(E0, E1, hs):
    h = SHAKE256.new()
    h.update(b'\x01')       # RO_CHALLENGE
    h.update(pack(curve(E0).j_invariant()))
    h.update(pack(curve(E1).j_invariant()))
    for h2,h3 in hs:
        h.update(h2 + h3)

    trits = []
    while len(trits) < len(hs):
        v, = h.read(int(1))
        if v >= 243:
            continue
        for _ in range(5):
            v,r = divmod(v, 3)
            trits.append({2:-1}.get(r,r))
    return trits[:len(hs)]

def canonical(E):
    R.<X> = FunctionField(E.base_field())
    eq = E.j_invariant() - EllipticCurve([0,X,0,1,0]).j_invariant()
    As = eq.numerator().roots(multiplicities=False)
    A = min(As, key = lambda el: (int(el[0]), int(el[1])))
    return E.isomorphism_to(curve(A))

def unpack_chain(A, chain):
    E = curve(A)
    for x in chain:
        K = E.lift_x(x)
        E = E.isogeny(K, algorithm='factored').codomain()
        E = canonical(E).codomain()
    return E

################################################################

class Parameters(NamedTuple):
    a: int
    b: int
    reps: int
    w: int
    h: int

    @property
    @cached_method
    def p(self):
        return 2^self.a * 3^self.b - 1

    @property
    @cached_method
    def i(self):
        return GF((self.p,2), 'i', modulus=[1,0,1]).gen()

    def read_element(self):
        s = input()
        if not re.match('^0x[0-9a-f]+,0x[0-9a-f]+$', s, re.IGNORECASE):
            invalid()
        c1, ci = (int(t,0) for t in s.split(','))
        return c1 + ci*self.i

class LeftResponse(NamedTuple):
    r2: bytes
    h3: bytes
    chain: list[RingElement]

    @staticmethod
    def read(param):
        return LeftResponse(read_bytes(), read_bytes(), [param.read_element() for _ in range(param.h)])

    def unpack(self, E0, _):
        j2 = unpack_chain(E0, self.chain).j_invariant()
        h2 = hash_hiding(self.r2, j2)
        return h2, self.h3

class BottomResponse(NamedTuple):
    r2: bytes
    r3: bytes
    E2: RingElement
    chain: list[RingElement]

    @staticmethod
    def read(param):
        return BottomResponse(read_bytes(), read_bytes(), param.read_element(), [param.read_element() for _ in range(param.w)])

    def unpack(self, *_):
        j2 = curve(self.E2).j_invariant()
        j3 = unpack_chain(self.E2, self.chain).j_invariant()
        h2 = hash_hiding(self.r2, j2)
        h3 = hash_hiding(self.r3, j3)
        return h2, h3

class RightResponse(NamedTuple):
    h2: bytes
    r3: bytes
    chain: list[RingElement]

    @staticmethod
    def read(param):
        return RightResponse(read_bytes(), read_bytes(), [param.read_element() for _ in range(param.h)])

    def unpack(self, _, E1):
        j3 = unpack_chain(E1, self.chain).j_invariant()
        h3 = hash_hiding(self.r3, j3)
        return self.h2, h3

################################################################

class Proof(NamedTuple):
    param: Parameters
    E0: RingElement
    E1: RingElement
    tags: list[int]
    responses: list[LeftResponse | BottomResponse | RightResponse]

    @staticmethod
    def read(param):
        E0 = param.read_element()
        E1 = param.read_element()
        tags = []
        responses = []
        while len(responses) < param.reps:
            tag = read_int()
            tags.append(tag)
            if tag == -1:
                responses.append(LeftResponse.read(param))
            elif tag == 0:
                responses.append(BottomResponse.read(param))
            elif tag == +1:
                responses.append(RightResponse.read(param))
            else:
                return None
        return Proof(param, E0, E1, tags, responses)

    def verify(proof):

        @parallel(os.cpu_count())
        def work(idx):
            return proof.responses[idx].unpack(proof.E0, proof.E1)

        hs = [None]*proof.param.reps

        for ((idx,),_),(h2,h3) in work(list(range(proof.param.reps))):
            print(f'{idx:3} {h2.hex()} {h3.hex()}', file=sys.stderr)
            hs[idx] = (h2,h3)

        if challenge(proof.E0, proof.E1, hs) != proof.tags:
            return None

        return proof.E0, proof.E1

##############################################################################################

#TODO should store these things in the proof file?
params = {
        'p434': Parameters(216, 137, 219, 4, 7),
        'p503': Parameters(250, 159, 219, 4, 7),
        'p610': Parameters(305, 192, 329, 4, 7),
        'p751': Parameters(372, 239, 438, 4, 7),
    }

if __name__ == '__main__' and '__file__' in globals():
    param = [v for k,v in params.items() if '--'+k in sys.argv[1:]]
    if len(param) != 1:
        print("Select the parameter set by passing exactly one of the following arguments:",
              ', '.join('--'+k for k in params), file=sys.stderr)
        exit(-1)
    param, = param

    proof = Proof.read(param)
    if proof is None:
        invalid()

    result = proof.verify()
    if result is None:
        invalid()

    print('\x1b[32mvalid proof\x1b[0m', file=sys.stderr)
    E0, E1 = result
    print(','.join(map(hex, E0)))
    print(','.join(map(hex, E1)))

