import hashlib

from utils import *

# Fr
r = 0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB26C5C28C859A99B3EEBCA9429212636B9DFF97634993AA4D6C381BC3F0057974EA099170FA13A4FD90776E240000001
fr = constants(r, 17, True, 'Fr')
Fr = FiniteField(r)
# Fq
q = 0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB117E776F218059DB80F0DA5CB537E38685ACCE9767254A4638810719AC425F0E39D54522CDD119F5E9063DE245E8001
fq = constants(q, 17, True, 'Fq')
Fq = FiniteField(q)
# Non-residue
non_residue = 13

# Fq2
Fq2.<A2> = GF(q**2, name='A2', modulus=x^2 - non_residue)

a = 2
b = 0x01373684A8C9DCAE7A016AC5D7748D3313CD8E39051C596560835DF0C9E50A5B59B882A92C78DC537E51A16703EC9855C77FC3D8BB21C8D68BB8CFB9DB4B8C8FBA773111C36C8B1B4E8F1ECE940EF9EAAD265458E06372009C9A0491678EF4


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
G2_order = 1755483545388786116744270475466687259186947712032004459714210070280389500116987496124098574823389466285978151140129779880473941566986064523874032812554001388754406832203678555787480693761694100603526668957377103500617895632245595170616424076193013069812655758706590804941865916654955161997668012879640089484267319662922323100918044574497880923438974476743556836643394916941627035296377155451512923542277406843394136761402218172188221567330050123715379201
G2_cofactor = int(G2_order / r)
print('G2 cofactor', G2_cofactor)
print(to_montgomery_vec(Fr(G2_cofactor), fr))
G2_cofactor_inv = Fr(G2_cofactor).inverse_of_unit()
print('G2 cofactor inv', G2_cofactor_inv)
print(to_montgomery_vec(G2_cofactor_inv, fr))