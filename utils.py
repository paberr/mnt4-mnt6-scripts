from sage.all import *

def biguint_num_bits(v):
    bits = 0

    while v != 0:
        v = v >> 1
        bits += 1

    return bits


### Convert BigUint into a vector of 64-bit limbs.
def biguint_to_u64_vec(v, limbs=0):
    v = Integer(v)
    m = 1 << 64
    ret = []

    while v > 0:
        ret.append(v % m)
        v = v >> 64

    while len(ret) < limbs:
        ret.append(0)

    return ret

def BigInteger(limbs):
    v = 0

    for limb in reversed(limbs):
        v = v << 64
        v |= limb
    
    return v

def vec_to_hex(v):
    return '[' + ', '.join([hex(val) for val in v]) + ']'

def constants(modulus, generator, should_print=False, name='Fq'):
    # The arithmetic in this library only works if the modulus*2 is smaller than the backing
    # representation. Compute the number of limbs we need.
    limbs = 1
    mod2 = modulus << 1 # modulus * 2
    cur = 1 << 64 # always 64-bit limbs for now
    while cur < mod2:
        limbs += 1
        cur = cur << 64

    modulus_num_bits = biguint_num_bits(modulus)

    # The number of bits we should "shave" from a randomly sampled reputation, i.e.,
    # if our modulus is 381 bits and our representation is 384 bits, we should shave
    # 3 bits from the beginning of a randomly sampled 384 bit representation to
    # reduce the cost of rejection sampling.
    repr_shave_bits = (64 * limbs) - modulus_num_bits

    # Compute R = 2**(64 * limbs) mod m
    r = (1 << (limbs * 64)) % modulus

    # modulus - 1 = 2^s * t
    s = 0
    t = modulus - 1;
    while t % 2 == 0:
        t = t >> 1
        s += 1

    # Compute 2^s root of unity given the generator
    root_of_unity = Integer(pow(generator, t, modulus))
    root_of_unity_v = biguint_to_u64_vec(to_montgomery(root_of_unity, modulus, r), limbs)
    generator_v = biguint_to_u64_vec(to_montgomery(generator, modulus, r), limbs)

    mod_minus_1_over_2 = (modulus - 1) >> 1
    mod_minus_1_over_2_v = biguint_to_u64_vec(mod_minus_1_over_2, limbs)

    # Compute R^2 mod m
    r2 = (r * r) % modulus
    r2_v = biguint_to_u64_vec(r2, limbs)

    r_v = biguint_to_u64_vec(r, limbs)
    modulus_v = biguint_to_u64_vec(modulus, limbs)

    # Compute -m^-1 mod 2**64 by exponentiating by totient(2**64) - 1
    u64_mod = 1 << 64
    inv = 1
    for _ in range(63):
        inv = (inv * inv) % u64_mod
        inv = (inv * modulus_v[0]) % u64_mod
    inv = (-inv) % u64_mod

    t_minus_1_over_2 = (t - 1) >> 1
    t_minus_1_over_2_v = biguint_to_u64_vec(t_minus_1_over_2, limbs)
    t_v = biguint_to_u64_vec(t, limbs)

    if should_print:
        print("""
impl FpParameters for {}Parameters {{
    type BigInt = BigInteger;

    // MODULUS = {}
    const MODULUS: BigInteger = BigInteger({});

    const MODULUS_BITS: u32 = {};

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = {};

    const R: BigInteger = BigInteger({});

    const R2: BigInteger = BigInteger({});

    const INV: u64 = {};

    const GENERATOR: BigInteger = BigInteger({});

    const TWO_ADICITY: u32 = {};

    const ROOT_OF_UNITY: BigInteger = BigInteger({});

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger({});

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    // T = (MODULUS - 1) / 2^S =
    // {}
    const T: BigInteger = BigInteger({});

    // (T - 1) / 2 =
    // {}
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger({});
}}""".format(name, modulus_v, vec_to_hex(modulus_v), modulus_num_bits, repr_shave_bits, vec_to_hex(r_v), vec_to_hex(r2_v), hex(inv), vec_to_hex(generator_v), s, vec_to_hex(root_of_unity_v), vec_to_hex(mod_minus_1_over_2_v), t_v, vec_to_hex(t_v), t_minus_1_over_2_v, vec_to_hex(t_minus_1_over_2_v)))

    return {
        'limbs': limbs,
        'modulus': modulus,
        'modulus_bits': modulus_num_bits,
        'repr_shave_bits': repr_shave_bits,
        'r': r,
        'r2': r2,
        'inv': inv,
        'generator': generator,
        'two_adicity': s,
        'root_of_unity': root_of_unity,
        'modulus_minus_one_div_two': mod_minus_1_over_2,
        't': t,
        't_minus_1_div_two': t_minus_1_over_2,
    }

bls12_381 = {
    'limbs': 4,
    'modulus': BigInteger([
        0xffffffff00000001,
        0x53bda402fffe5bfe,
        0x3339d80809a1d805,
        0x73eda753299d7d48,
    ]),
    'modulus_bits': 255,
    'repr_shave_bits': 1,
    'r': BigInteger([
        0x1fffffffe,
        0x5884b7fa00034802,
        0x998c4fefecbc4ff5,
        0x1824b159acc5056f,
    ]),
    'r2': BigInteger([
        0xc999e990f3f29c6d,
        0x2b6cedcb87925c23,
        0x5d314967254398f,
        0x748d9d99f59ff11,
    ]),
    'inv': 0xfffffffeffffffff,
    'generator': BigInteger([
        0xefffffff1,
        0x17e363d300189c0f,
        0xff9c57876f8457b0,
        0x351332208fc5a8c4,
    ]),
    'two_adicity': 32,
    'root_of_unity': BigInteger([
        0xb9b58d8c5f0e466a,
        0x5b1b4c801819d7ec,
        0xaf53ae352a31e64,
        0x5bf3adda19e9b27b,
    ]),
    'modulus_minus_one_div_two': BigInteger([
        0x7fffffff80000000,
        0xa9ded2017fff2dff,
        0x199cec0404d0ec02,
        0x39f6d3a994cebea4,
    ]),
    't': BigInteger([
        0xfffe5bfeffffffff,
        0x9a1d80553bda402,
        0x299d7d483339d808,
        0x73eda753,
    ]),
    't_minus_1_div_two': BigInteger([
        0x7fff2dff7fffffff,
        0x4d0ec02a9ded201,
        0x94cebea4199cec04,
        0x39f6d3a9,
    ]),
}

def to_montgomery(x, m, R):
    return (x * R) % m

def to_montgomery_vec(x, fp):
    mont = (x * fp['r']) % fp['modulus']
    return biguint_to_u64_vec(mont, fp['limbs'])

def constants_to_montgomery(c):
    c['generator'] = to_montgomery(c['generator'], c['modulus'], c['r'])
    c['root_of_unity'] = to_montgomery(c['root_of_unity'], c['modulus'], c['r'])
    return c


#########
# Tests #
#########
bls12_381_computed = constants_to_montgomery(constants(52435875175126190479447740508185965837690552500527637822603658699938581184513, 7, False))
assert bls12_381_computed == bls12_381
