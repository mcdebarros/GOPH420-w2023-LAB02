import numpy as np

def newton_raphson(zeta_0, g, dgdx, eps_s):

    """Iterates newton-raphson method until error spec reached.
    Parameters
    ----------
    x0 : int or float
        Initial guess of the root
    f : int or float
        Function for root finding
    fp : int or float
        First derivative of function
    eps_s : int or float
        Stopping criteria for desired error

    Returns
    -------
    eps_a_array : array-like
        Error at each iteration
    
    """
    eps_a = 1
    zeta_r = zeta_0
    it = 0
    eps_a_arr = []


    while eps_a > eps_s:

        dx = float((-g(zeta_r))/(dgdx(zeta_r)))
        zeta_r += dx

        eps_a = np.abs(dx/zeta_r)
        eps_a_arr = np.array(eps_a_arr)
        eps_a_arr.append(eps_a)
        it = it + 1

    return (zeta_r, it, eps_a_arr)