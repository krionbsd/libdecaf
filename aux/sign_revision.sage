F = GF(2^255-19)
dM = F(-121665)
d = F(-121665/121666)
d_2 = -d
ii = sqrt(F(-1))
def lobit(x): return int(x) & 1
def hibit(x): return lobit(2*x)
cofactor = 8
a = -1

magic = sqrt(dM-1)
magic_2 = sqrt(1-dM)
if lobit(magic): magic = -magic

class DecafException(Exception): pass
class NotOnCurveException(DecafException): pass
class InvalidDecafException(DecafException): pass
class UnimplementedException(DecafException): pass


def decode_to_isomorphic_edwards(s):
    if s == 0: return (0,1)
    if hibit(s): raise InvalidDecafException("s has high bit")
    if not is_square(s^4 + (2-4*dM)*s^2 + 1): raise NotOnCurveException("Not on curve")
    if s^2+1 == 0: raise InvalidDecafException("s = i")
    
    t = sqrt(s^4 + (2-4*dM)*s^2 + 1)
    # t=0 -> (s^2+1)^2 = 4dM -> dM is square -> incomplete

    y = (1-s^2)/(1+s^2)
    x = 2*magic_2*s/t
    if hibit(x): x = -x
    if cofactor == 8 and hibit(y): raise InvalidDecafException("y must be positive")
    
    assert y^2 + x^2 == 1 + d_2*x^2*y^2
    return (x,y)

def decode_to_isogenous_edwards(s,a=1):
    if s == 0: return (0,1)
    if hibit(s): raise InvalidDecafException("s has high bit")
    if not is_square(s^4 + (2-4*dM)*s^2 + 1): raise NotOnCurveException("Not on curve")
    if s^2+1 == 0: raise InvalidDecafException("s = i")
    
    t = sqrt(s^4 + (2-4*dM)*s^2 + 1)
    
    x = 2*s / (1+a*s^2)
    y = (1-a*s^2) / t

    xIso = 2*magic_2*s/t
    if hibit(xIso): x = -x
    if cofactor == 8:
        yIso = (1-a*s^2)/(1+a*s^2)
        if hibit(yIso): raise InvalidDecafException("y must be positive")
    dM1 = dM
    if a == -1: dM1 -= 1
    assert y^2 + a*x^2 == 1 + dM1*x^2*y^2
    return (x,y)
    
    
    





def eddsa_decode(y_ser):
    """
    Recover x,y corresponding to the eddsa standard representation.
    """
    hi = int(y_ser) & 2^255
    y = F(y_ser-hi)
    
    x = sqrt((y^2-1)/(d*y^2+1))                                                               
    if int(x) & 1: x = -x
    if hi: x = -x                                                                    
    assert y^2 - x^2 == 1+d*x^2*y^2
    return (x,y)

def eddsa_to_decaf(x,y):
    """
    Converts an EdDSA point to a Decaf representation, in a manner compatible
    with libdecaf.
    
    The input point must be even.
    
    Note well!  Decaf does not represent the cofactor information of a point.
    So e2d(d2e(s)) = s, but d2e(e2d(x,y)) might not be (x,y).
    """
    if x*y == 0: return 0 # This will happen anyway with straightforward square root trick
    if not is_square((1-y)/(1+y)): raise Exception("Unimplemented: odd point in eddsa_to_decaf")
    if hibit(magic/(x*y)): (x,y) = (ii*y,ii*x)
    if hibit(2*magic/x): y = -y
    s = sqrt((1-y)/(1+y))
    if hibit(s): s = -s
    return s

def isqrt_trick(to_isr,to_inv):
    """
    The "inverse square root" trick:
    Return 1/sqrt(to_isr), 1/to_inv.
    Also return their product because that turns out to be useful.
    """
    to_sqrt = to_isr*to_inv^2
    
    if to_sqrt == 0: return 0,0,0 # This happens automatically in C; just to avoid problems in SAGE
    if not is_square(to_sqrt): raise Exception("Not square in isqrt_trick!")
    
    isr_times_inv = 1/sqrt(to_sqrt)
    isr = isr_times_inv * to_inv
    inv = isr_times_inv * isr * to_isr
    
    assert isr^2 == 1/to_isr
    assert inv == 1/to_inv
    return isr, inv, isr_times_inv
    

def eddsa_to_decaf_opt(x,y,z=None):
    """
    Optimized version of eddsa_to_decaf.   Uses only one isqrt.
    There's probably some way to further optimize if you have a T-coord,
    but whatever.
    """
    if z is None:
        # Pretend that we're in projective coordinates
        z = F.random_element()
        x *= z
        y *= z
    
    isr,inv,isr_times_inv = isqrt_trick(z^2-y^2,x*y)
    minv = inv*magic*z
    
    rotate = hibit(minv*z)
    if rotate:
        isr = isr_times_inv*(z^2-y^2)*magic
        y = ii*x
    
    if hibit(2*minv*y) != rotate: y = -y
    s = (z-y) * isr
    
    if hibit(s): s = -s
    return s


def decaf_to_eddsa(s):
    """
    Convert a Decaf representation to an EdDSA point, in a manner compatible
    with libdecaf.
    
    Note well!  Decaf does not represent the cofactor information of a point.
    So e2d(d2e(s)) = s, but d2e(e2d(x,y)) might not be (x,y).
    """
    if s == 0: return (0,1)
    if s < 0 or s >= F.modulus(): raise Exception("out of field!")
    s = F(s)
    if hibit(s): raise Exception("invalid: s has high bit")
    if not is_square(s^4 + (2-4*dM)*s^2 + 1): raise Exception("invalid: not on curve")
    
    t = sqrt(s^4 + (2-4*dM)*s^2 + 1)/s
    if hibit(t): t = -t
    y = (1-s^2)/(1+s^2)
    x = 2*magic/t
    
    if y == 0 or lobit(t/y): raise Exception("invalid: t/y has high bit")
    
    assert y^2 - x^2 == 1+d*x^2*y^2
    return (x,y)

def decaf_to_eddsa_opt(s):
    """
    Convert a Decaf representation to an EdDSA point, in a manner compatible
    with libdecaf.  Optimized to use only one invsqrt.
    
    This function would be slightly simpler if we didn't want to decode to affine.
    """
    if s < 0 or s >= F.modulus(): raise Exception("out of field!")
    s = F(s)
    if hibit(s): raise Exception("Invalid: s has high bit")
    if s == 0: return (0,1)
    
    curve_eqn = s^4 + (2-4*dM)*s^2 + 1
    isr,inv,isr_times_inv = isqrt_trick(curve_eqn,s*(1-s^2)*(1+s^2))
    if isr == 0: raise Exception("Invalid: nonstandard encoding of zero")
    
    x = 2 * magic * s * isr
    y = (1-s^2)^2 * s * inv
    
    tmp = isr_times_inv * curve_eqn * (1+s^2) # = sqrt(curve_eqn) / s / (1-s^2)
    hibit_t = hibit(tmp * (1-s^2))
    if hibit_t: x = -x
    if lobit(tmp * (1+s^2)) != hibit_t: raise Exception("invalid: bits don't match")
    
    assert y^2 - x^2 == 1+d*x^2*y^2
    return (x,y)
