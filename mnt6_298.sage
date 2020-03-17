import hashlib

from utils import *

# Fr
r = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
fr = constants(r, 17, True, 'Fr')
Fr = FiniteField(r)
# Fq
q = 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
fq = constants(q, 17, True, 'Fq')
Fq = FiniteField(q)
# Non-residue
non_residue = Fq(BigInteger([
    0x58eefd67fea995ca,
    0x12f14affbb33a004,
    0x4780323da44ac69b,
    0x88acf9bea707eed9,
    0x14bbbb859e8,
])) / Fq(fq['r'])
non_residue = Integer(non_residue)
print('non_residue', non_residue)

# Fq3
Fq3.<A3> = GF(q**3, name='A3', modulus=x^3 - non_residue)
A3s = A3*A3

a = Fq(BigInteger([
    0xb9b2411bfd0eafef,
    0xc61a10fadd9fecbd,
    0x89f128e59811f3fb,
    0x980c0f780adadabb,
    0x9ba1f11320,
])) / Fq(fq['r'])
b = Fq(BigInteger([
	0xa94cb16ed8e733b,
	0xe1ed15e8119bae6,
	0xae927592157c8121,
	0x990dbcbc6661cf95,
	0xecff0892ef,
])) / Fq(fq['r'])


g1_coeff_a = Fq(a)
print('G1 coeff a:', g1_coeff_a)
print(to_montgomery_vec(g1_coeff_a, fq))

g1_coeff_b = Fq(b)
print('G1 coeff b:', g1_coeff_b)
print(to_montgomery_vec(g1_coeff_b, fq))

g2_coeff_a = Fq3(g1_coeff_a * A3s)
print('G2 coeff a', vector(g2_coeff_a))
print(tuple(to_montgomery_vec(i, fq) for i in vector(g2_coeff_a)))

g2_coeff_b = Fq3(non_residue * g1_coeff_b)
print('G2 coeff b', vector(g2_coeff_b))
print(tuple(to_montgomery_vec(i, fq) for i in vector(g2_coeff_b)))

G1 = EllipticCurve(Fq, [g1_coeff_a, g1_coeff_b])
G2 = EllipticCurve(Fq3, [g2_coeff_a, g2_coeff_b])

G1_cofactor = 1
print('G1 cofactor', G1_cofactor)
G1_cofactor_inv = Fr(G1_cofactor).inverse_of_unit()
print('G1 cofactor inv', G1_cofactor_inv)

# Computed via G2.order() and put here for speed
G2_order = 107797360357109903430794490309592072278927783178002957258386272738249831275482847380369870112378590485473667408581494432820181343609082473450588258860186349169969034947844658697099498402986681525178913967480073918525963503551219889381053812647119750622329389425224581120
G2_cofactor = int(G2_order / r)
assert G2_cofactor == 226502022472576270196498690498308461791828762732602586162207535351960270082712694977333372361549082214519252261735048131889018501404377856786623430385820659037970876666767495659520
print('G2 cofactor', G2_cofactor)
print(to_montgomery_vec(Fr(G2_cofactor), fr))
G2_cofactor_inv = Fr(G2_cofactor).inverse_of_unit()
print('G2 cofactor inv', G2_cofactor_inv)
print(to_montgomery_vec(G2_cofactor_inv, fr))