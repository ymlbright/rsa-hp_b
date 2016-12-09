import time

debug = True

# display matrix picture with 0 and X
def matrix_overview(BB, bound):
    for ii in range(BB.dimensions()[0]):
        a = ('%02d ' % ii)
        for jj in range(BB.dimensions()[1]):
            a += '0' if BB[ii,jj] == 0 else 'X'
            a += ' '
        if BB[ii, ii] >= bound:
            a += '~'
        print a

def coppersmith_howgrave_univariate(pol, modulus, beta, mm, tt, XX):
    """
    Coppersmith revisited by Howgrave-Graham
    
    finds a solution if:
    * b|modulus, b >= modulus^beta , 0 < beta <= 1
    * |x| < XX
    """
    #
    # init
    #
    dd = pol.degree()
    nn = dd * mm + tt

    #
    # checks
    #
    if not 0 < beta <= 1:
        raise ValueError("beta should belongs in (0, 1]")

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    #
    # calculate bounds and display them
    #
    """
    * we want to find g(x) such that ||g(xX)|| <= b^m / sqrt(n)
    * we know LLL will give us a short vector v such that:
    ||v|| <= 2^((n - 1)/4) * det(L)^(1/n)
    * we will use that vector as a coefficient vector for our g(x)
    
    * so we want to satisfy:
    2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)
    
    so we can obtain ||v|| < N^(beta*m) / sqrt(n) <= b^m / sqrt(n)
    (it's important to use N because we might not know b)
    """
    if debug:
        # t optimized?
        print "\n# Optimized t?\n"
        print "we want X^(n-1) < N^(beta*m) so that each vector is helpful"
        cond1 = RR(XX^(nn-1))
        print "* X^(n-1) = ", cond1
        cond2 = pow(modulus, beta*mm)
        print "* N^(beta*m) = ", cond2
        print "* X^(n-1) < N^(beta*m) \n-> GOOD" if cond1 < cond2 else "* X^(n-1) >= N^(beta*m) \n-> NOT GOOD"
        
        # bound for X
        print "\n# X bound respected?\n"
        print "we want X <= N^(((2*beta*m)/(n-1)) - ((delta*m*(m+1))/(n*(n-1)))) / 2 = M"
        print "* X =", XX
        cond2 = RR(modulus^(((2*beta*mm)/(nn-1)) - ((dd*mm*(mm+1))/(nn*(nn-1)))) / 2)
        print "* M =", cond2
        print "* X <= M \n-> GOOD" if XX <= cond2 else "* X > M \n-> NOT GOOD"

        # solution possible?
        print "\n# Solutions possible?\n"
        detL = RR(modulus^(dd * mm * (mm + 1) / 2) * XX^(nn * (nn - 1) / 2))
        print "we can find a solution if 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)"
        cond1 = RR(2^((nn - 1)/4) * detL^(1/nn))
        print "* 2^((n - 1)/4) * det(L)^(1/n) = ", cond1
        cond2 = RR(modulus^(beta*mm) / sqrt(nn))
        print "* N^(beta*m) / sqrt(n) = ", cond2
        print "* 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n) \n-> SOLUTION WILL BE FOUND" if cond1 < cond2 else "* 2^((n - 1)/4) * det(L)^(1/n) >= N^(beta*m) / sqroot(n) \n-> NO SOLUTIONS MIGHT BE FOUND (but we never know)"

        # warning about X
        print "\n# Note that no solutions will be found _for sure_ if you don't respect:\n* |root| < X \n* b >= modulus^beta\n"
    
    #
    # Coppersmith revisited algo for univariate
    #

    # change ring of pol and x
    polZ = pol.change_ring(ZZ)
    x = polZ.parent().gen()

    # compute polynomials
    gg = []
    for ii in range(mm):
        for jj in range(dd):
            gg.append((x * XX)**jj * modulus**(mm - ii) * polZ(x * XX)**ii)
    for ii in range(tt):
        gg.append((x * XX)**ii * polZ(x * XX)**mm)
    
    # construct lattice B
    BB = Matrix(ZZ, nn)

    for ii in range(nn):
        for jj in range(ii+1):
            BB[ii, jj] = gg[ii][jj]

    # display basis matrix
    if debug:
        matrix_overview(BB, modulus^mm)

    # LLL
    BB = BB.LLL()

    # transform shortest vector in polynomial    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[0, ii] / XX**ii

    # factor polynomial
    potential_roots = new_pol.roots()
    print "potential roots:", potential_roots

    # test roots
    roots = []
    for root in potential_roots:
        if root[0].is_integer():
            result = polZ(ZZ(root[0]))
            if gcd(modulus, result) >= modulus^beta:
                roots.append(ZZ(root[0]))

    # 
    return roots

############################################
# Test on Stereotyped Messages
##########################################    

# print "//////////////////////////////////"
# print "// TEST 1"
# print "////////////////////////////////"

# # RSA gen options (for the demo)
# length_N = 1024  # size of the modulus
# Kbits = 200      # size of the root
# e = 3

# # RSA gen (for the demo)
# p = next_prime(2^int(round(length_N/2)))
# q = next_prime(p)
# N = p*q
# ZmodN = Zmod(N);

# # Create problem (for the demo)
# K = ZZ.random_element(0, 2^Kbits)
# Kdigits = K.digits(2)
# M = [0]*Kbits + [1]*(length_N-Kbits); 
# for i in range(len(Kdigits)):
#     M[i] = Kdigits[i]
# M = ZZ(M, 2)
# C = ZmodN(M)^e

# # Problem to equation (default)
# P.<x> = PolynomialRing(ZmodN) #, implementation='NTL')
# pol = (2^length_N - 2^Kbits + x)^e - C
# dd = pol.degree()

# # Tweak those
# beta = 1                                # b = N
# epsilon = beta / 7                      # <= beta / 7
# mm = ceil(beta**2 / (dd * epsilon))     # optimized value
# tt = floor(dd * mm * ((1/beta) - 1))    # optimized value
# XX = ceil(N**((beta**2/dd) - epsilon))  # optimized value

# # Coppersmith
# start_time = time.time()
# roots = coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX)

# # output
# print "\n# Solutions"
# print "we want to find:",str(K)
# print "we found:", str(roots)
# print("in: %s seconds " % (time.time() - start_time))
# print "\n"

# ############################################
# # Test on Factoring with High Bits Known
# ##########################################
# print "//////////////////////////////////"
# print "// TEST 2"
# print "////////////////////////////////"

# # RSA gen
# length_N = 1024;
# p = next_prime(2^int(round(length_N/2)));
# q = next_prime( round(pi.n()*p) );
# N = p*q;

# # qbar is q + [hidden_size_random]
# hidden = 200;
# diff = ZZ.random_element(0, 2^hidden-1)
# qbar = q + diff; 

# F.<x> = PolynomialRing(Zmod(N), implementation='NTL'); 
# pol = x - qbar
# dd = pol.degree()

# # PLAY WITH THOSE:
# beta = 0.5                             # we should have q >= N^beta
# epsilon = beta / 7                     # <= beta/7
# mm = ceil(beta**2 / (dd * epsilon))    # optimized
# tt = floor(dd * mm * ((1/beta) - 1))   # optimized
# XX = ceil(N**((beta**2/dd) - epsilon)) # we should have |diff| < X

# # Coppersmith
# start_time = time.time()
# roots = coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX)

# # output
# print "\n# Solutions"
# print "we want to find:", qbar - q
# print "we found:", roots
# print("in: %s seconds " % (time.time() - start_time))

print "//////////////////////////////////"
print "// TEST 3"
print "////////////////////////////////"

# RSA key
p = 0xfc3e98501e287508d9382d22d698842da518d1719f4af707a55977204bf7500c08bab424ed905f8948b8703071efc3112225698523cc9c7dccf1224fe5e8ea18e491ca42131d5776d8f2c20a1ddf561a9373b8ff5af3cd1408c44d070fc606027e39ca804a973454b8a7fcadb8d49783acee50f3ab8f16660d32c9f28799833b
q = 0xc7c90b300ab1cb7bb5f2e19e2fd7749a260703f50525452689f4435daa90f1cba8c8deb81bf93c58e9f85023f53265709593eb05fb061d6164cd4070a73b596b0264d32d4328626141827896d89f00cec4d2ded954411bb2d173b5eeb23651ee2e36f1da69c9ab88c3a7cc6fa8effc4e81112021be0507e4860283ee8fde6f39
pbar = 0xfc3e98501e287508d9382d22d698842da518d1719f4af707a55977204bf7500c08bab424ed905f8948b8703071efc3112225698523cc9c7dccf1224fe5e8ea18ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
N = 0xc4dac091e52ac2c9221f434f022ff5b9797243685ae620d8a91935195402c08afc3e98501e287508d9382d22d698842da518d1719f4af707a55977204bf7500c08bab424ed905f8948b8703071efc3112225698523cc9c7dccf1224fe5e8ea18e69042fc0f671476c834491d59bb28c853fc1c62e6fadf30d94e9956a2da1c438718c58e9cf0217360193397c91a265cb080263b1c619ee0a224bc83299015630fef3ab9f33bf7a831390928502edd2370a5a63b52a1b7ecbbc6c563ff65e54e1693687c4b50c792352f489ff81f59d80a5d523b2cdd7ecef1952effb762ea486919ecd444fe59e21b1be76c669a90fda9bd7aaec655cf991744a913833ecd23

F.<x> = PolynomialRing(Zmod(N), implementation='NTL'); 
pol = x - pbar
dd = pol.degree()

# PLAY WITH THOSE:
beta = 0.5                             # we should have q >= N^beta
epsilon = beta / 35                    # <= beta/35
mm = ceil(beta**2 / (dd * epsilon))    # optimized
tt = floor(dd * mm * ((1/beta) - 1))   # optimized
XX = ceil(N**((beta**2/dd) - epsilon)) # we should have |diff| < X

# Coppersmith
start_time = time.time()
roots = coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX)

# output
print "\n# Solutions"
print "we want to find:", p - pbar
print "we found:", roots
print("in: %s seconds " % (time.time() - start_time))

"""
Setting permissions of DOT_SAGE directory so only you can read and write it.
//////////////////////////////////
// TEST 1
////////////////////////////////

# Optimized t?

we want X^(n-1) < N^(beta*m) so that each vector is helpful
* X^(n-1) =  5.26588450218277e469
* N^(beta*m) =  5809605995369958062859502533304574370686975176362895236661486152287203730997110225737336044533118407251326157754980517443990529594540047121662885672187318379170893091380779314421226637138275349470290853160784521096650764992953261892033155088964701773548941402064393684794946750211086736290671577104533659642057118917143584402358425763453268564512919742021495562856804237181830944118597652633110777860501395758734009382103277519347965512270318311650775740609245712375989187523162575910893401371870662641830320681166427150335438145263362770625995772981051308652613471165711681334529841626638582686839433337262143010780233438621153258840860120619228367830876704643942342141894676446795044841146659969772821383582256097597403390152196838898111302490307579997628617454989618216755310133769183576086459390337713555959065986325737177623859347069743089552382311434116472575356994213750090774944532102495872113514499697576541330780931
* X^(n-1) < N^(beta*m) 
-> GOOD

# X bound respected?

we want X <= N^(((2*beta*m)/(n-1)) - ((delta*m*(m+1))/(n*(n-1)))) / 2 = M
* X = 51901978826686593576182134394350049796327272495710410629191
* M = 5.78960446186581e76
* X <= M 
-> GOOD

# Solutions possible?

we can find a solution if 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)
* 2^((n - 1)/4) * det(L)^(1/n) =  9.38051702189965e851
* N^(beta*m) / sqrt(n) =  1.93653533178999e924
* 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n) 
-> SOLUTION WILL BE FOUND

# Note that no solutions will be found _for sure_ if you don't respect:
* |root| < X 
* b >= modulus^beta

00 X 0 0 0 0 0 0 0 0 ~
01 0 X 0 0 0 0 0 0 0 ~
02 0 0 X 0 0 0 0 0 0 ~
03 X X X X 0 0 0 0 0 
04 0 X X X X 0 0 0 0 
05 0 0 X X X X 0 0 0 
06 X X X X X X X 0 0 
07 0 X X X X X X X 0 
08 0 0 X X X X X X X 
potential roots: [(400069133414591588306107259527152610081926740160271826127377, 1)]

# Solutions
we want to find: 400069133414591588306107259527152610081926740160271826127377
we found: [400069133414591588306107259527152610081926740160271826127377]
in: 0.146455049515 seconds 


//////////////////////////////////
// TEST 2
////////////////////////////////

# Optimized t?

we want X^(n-1) < N^(beta*m) so that each vector is helpful
* X^(n-1) =  8.70626317072222e385
* N^(beta*m) =  3.18956065351443e617
* X^(n-1) < N^(beta*m) 
-> GOOD

# X bound respected?

we want X <= N^(((2*beta*m)/(n-1)) - ((delta*m*(m+1))/(n*(n-1)))) / 2 = M
* X = 13622652689555162275328740864963325782847721740190089216
* M = 7.24576134690965e65
* X <= M 
-> GOOD

# Solutions possible?

we can find a solution if 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)
* 2^((n - 1)/4) * det(L)^(1/n) =  2.73243713251504e579
* N^(beta*m) / sqrt(n) =  1.12767998355292e617
* 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n) 
-> SOLUTION WILL BE FOUND

# Note that no solutions will be found _for sure_ if you don't respect:
* |root| < X 
* b >= modulus^beta

00 X 0 0 0 0 0 0 0 ~
01 X X 0 0 0 0 0 0 
02 X X X 0 0 0 0 0 
03 X X X X 0 0 0 0 
04 X X X X X 0 0 0 
05 0 X X X X X 0 0 
06 0 0 X X X X X 0 
07 0 0 0 X X X X X 
potential roots: [(42121870893450634577463914985889299119866228583627912396576170307551916037987547771260822964333724590046064811464937900519371437792981097422481012078311459, 1), (1183580336003521515246309350207881012079049629237975884396458, 3), (333148247511781621682357144377220120467138093012617253664782517629913979199/281474976710656, 3)]

# Solutions
we want to find: 1183580336003521515246309350207881012079049629237975884396458
we found: [42121870893450634577463914985889299119866228583627912396576170307551916037987547771260822964333724590046064811464937900519371437792981097422481012078311459, 1183580336003521515246309350207881012079049629237975884396458]
in: 0.0174250602722 seconds 
"""