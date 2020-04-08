import hashlib

from utils import *

# Fr
r = 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
fr = constants(r, 17, True, 'Fr')
Fr = FiniteField(r)
# Fq
q = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
fq = constants(q, 17, True, 'Fq')
Fq = FiniteField(q)
# Non-residue
non_residue = 17

# Fq2
Fq2.<A2> = GF(q**2, name='A2', modulus=x^2 - non_residue)

a = 2
b = 423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685


g1_coeff_a = Fq(a)
print('G1 coeff a:', g1_coeff_a)
print(to_montgomery_vec(g1_coeff_a, fq))

g1_coeff_b = Fq(b)
print('G1 coeff b:', g1_coeff_b)
print(to_montgomery_vec(g1_coeff_b, fq))

g2_coeff_a = Fq2(non_residue * g1_coeff_a)
print('G2 coeff a', vector(g2_coeff_a))
print(tuple(to_montgomery_vec(i, fq) for i in vector(g2_coeff_a)))

g2_coeff_b = Fq2(non_residue * g1_coeff_b * A2)
print('G2 coeff b', vector(g2_coeff_b))
print(tuple(to_montgomery_vec(i, fq) for i in vector(g2_coeff_b)))

G1 = EllipticCurve(Fq, [g1_coeff_a, g1_coeff_b])
G2 = EllipticCurve(Fq2, [g2_coeff_a, g2_coeff_b])

G1_cofactor = 1
print('G1 cofactor', G1_cofactor)
G1_cofactor_inv = Fr(G1_cofactor).inverse_of_unit()
print('G1 cofactor inv', G1_cofactor_inv)

# Computed via G2.order() and put here for speed
G2_order = 226502022472576270196498690498308461791828763717577836213694428486838284845848266994921634460577423496873194324622171087135390176549150010950410620094633001115253228252990372839425
G2_cofactor = int(G2_order / r)
print('G2 cofactor', G2_cofactor)
print(biguint_to_u64_vec(G2_cofactor))
G2_cofactor_inv = Fr(G2_cofactor).inverse_of_unit()
print('G2 cofactor inv', G2_cofactor_inv)
print(to_montgomery_vec(G2_cofactor_inv, fr))