import numpy as np

def newton_raphson(x0, fu, fp, eps_s = 1e-8):

    """Iterates newton-raphson method until error spec reached.
    Parameters
    ----------
    xo : int or float
        Initial guess of the root
    f : int or float
        Function for root finding
    fp : int or float
        First derivative of function
    eps_s : int or float
        Stopping criteria for desired error

        

    Returns
    -------
    int or float or array_like
        The sum of x and y.
    """
    eps_a = 1
    xi = x0
    it = 0
    

    while eps_a > eps_s:

        dx = float(-fu(xi)/fp(xi))
        xi += dx

        eps_a = np.abs(dx/xi)
        it = it + 1

    return (eps_a, xi, it)