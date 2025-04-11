from sage.stats.distributions.discrete_gaussian_integer import \
DiscreteGaussianDistributionIntegerSampler

from sage.modules.free_module_integer import IntegerLattice

def MatInfNorm(A):
    return max([c.norm(2) for c in A.columns()])

#A = random_matrix(ZZ, 3, 3)

n = 3 
m = 4

# beta should be at least sqrt(m) to make sure the SIS instance has solution
beta = ceil(sqrt(m)) * m 

# We took aprox factor of 2 for now
q = next_prime(ceil(4*n*sqrt(m)* beta))

A = matrix(ZZ, n,n, [0,1,2,3,4,5,6,7,9])
#B = matrix(ZZ, n,n, [1,2,3,0,1,4,5,6,0])
B = random_matrix(ZZ, n, n)
A_inv = A.inverse()

s = q*MatInfNorm(B)/(4*sqrt(n*m) * beta)


def random_guassian_vector(s):
    sigma = s 
    D = DiscreteGaussianDistributionIntegerSampler(sigma)
    return tuple(D() for _ in range(n))

# Vodoo in the proff, some implementation that works for small dimensions
def sis_oracle(a_i): 
    A = matrix(IntegerModRing(q), a_i)
    ns = A.kernel()

    L = IntegerLattice(ns)

    v = L.shortest_vector()
    print(f"Oracle O solved SIS! shortest vector is {v}")
    return v

#we do it with a trapdoor now, so we start with small x and build the problem that this SIS instance will follow
trapdoor = vector(ZZ, [1, 0, 0])
sis_cols = A.columns()
sis_cols.append(-A * trapdoor + vector(ZZ , [q, 2 *q, 3 * q]))
sis_instance = matrix(ZZ, sis_cols).transpose()


def trapdoor_guasian_vector(s): 
    return [B * c  for c in sis_instance.columns()]

print(f"s,q is {s,q}")

def core_step(B):
    x_i = [random_guassian_vector(s) for i in range(m)]
    #x_i = trapdoor_guasian_vector(s)

    ## Generate a SIS instance now
    c_i = []
    v_i = []
    for i in range(m): 
        c_ii = (B.inverse() * vector(x_i[i])).apply_map(lambda x: round(x))
        c_i.append(c_ii)
        v_i.append(B* c_ii)

    a_i = [x.apply_map(lambda y: y % q) for x in c_i]

    if matrix(ZZ, a_i).rank() < m:
        ## Now look for SIS solutions with this problem
        print("We have found a successfull  SIS instance!")
        z = sis_oracle(a_i)
        if (z.norm(2) > beta):
            return 0

        print(f"Found decent SIS solution to approximation factor {beta}! continuing reduction {z}")

        L = IntegerLattice(B)
        new_lat_vec = vector(ZZ, [0] * n)
        c_i_z_i = vector(ZZ, [0] *n)
        for i in range(m): 
            new_lat_vec += 1/q * (v_i[i] * z[i])
            c_i_z_i += c_i[i] * z[i]

        print(f"generating a new shorter lattice vector: {new_lat_vec}")
        coeff = B.solve_right(new_lat_vec)

        print(f"V is in the Lattice becuase:\n {B}*{coeff}=V")
        return 1

while core_step(B) != 1:
    pass
    
print(IntegerLattice(B).shortest_vector())
