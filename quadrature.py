import numpy as np

"Numerical integration using Simpson's rule with n steps."
    
def quadSimpson(f, a, b, n):

    h = (b - a) / n #step size
    s = f(a) + f(b)
    
    # such as for f(x)= x,x^2,x^3,x^4,x^5
    #iterates over odd-indexed points in the interval.(x^1,x^3,x^5...)
    # step size 2 , iterate for odd values starts from 1,3,5..
    for i in range(1, n, 2):
        s += 4 * f(a + i * h)
        
    #iterates over even-indexed points in the interval.(x^2,x^4...)
    # step size 2 , iterate for even values starts from 2,4,6...

    for i in range(2, n, 2):
        s += 2 * f(a + i * h)

    return s * h / 3


"Numerical Integration using Adaptive Quadrature"
    
def quadAdaptiv(f, a, b, tol, hmin = 1.e-6):
    
    nodes = [a, b]  # List of evaluation centers
    Q = quadAdaptivRec(f, a, b, tol, hmin, nodes)
    
    #Sort the list of evaluation centers(nodes) in ascending order
    return (sorted(nodes), Q)

def quadAdaptivRec(f, a, b, tol, hmin, nodes):
    m = (a + b)/2
    nodes.append(m) # Add the midpoint to the list of evaluation centers
    h = b - a

    if h < hmin:
        raise ValueError('The minimum step size has been undercut.')
        
    #Calculate Error    
    qTrapez  = h/2 * (f(a) + f(b))
    qSimpson = h/6 * (f(a) + 4*f(m) + f(b))
    err = np.abs(qSimpson - qTrapez)
    
# (else) tol > err , keep recursing the process this will produce more nodes after each step

    if err <= tol:
        return qSimpson
    else:
        left  = quadAdaptivRec(f, a, m, tol/2, hmin, nodes)
        right = quadAdaptivRec(f, m, b, tol/2, hmin, nodes) 
        return left + right
