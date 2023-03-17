import numpy as np

import matplotlib.pyplot as plt

from newton_raphson.rootfinder import newton_raphson

def main():

    beta_1 = 1900
    beta_2 = 3200
    rho_1 = 1800
    rho_2 = 2500
    h = 4000
    eps_s = 1e-8
    zeta_max = (h**2) * ((beta_1**-2) - (beta_2**-2))
    f = []
    zeta_r = []
    f_mode = []
    zeta_mode = []

    for n in range (6):

        f_new = (((2*n) + 1)*(0.25/zeta_max))
        f.append(f_new)

        zeta_r_list = []

        for k in range (n+1):

            zeta_0 = (((2*k) + 1) * (0.25/f_new)) - (zeta_max * 1e-3)

            def g(zeta):
        
                return ((rho_2/rho_1)*((np.sqrt((zeta_max) - (zeta**2)))/zeta)) - (np.tan(2*np.pi*f[n]*zeta))
    
            def dgdx(zeta):

                return ((rho_2/rho_1)*zeta_max/(beta_1**2 * beta_2**2 * zeta**2 * np.sqrt((zeta_max)) - zeta**2)) - ((2*np.pi*f[n])*((1/np.cos(2*np.pi*f[n]*zeta))**2))
            
            zeta_r_new = newton_raphson(zeta_0, g, dgdx, eps_s)

            zeta_r_list.append(zeta_r_new)
        
        zeta_r.append(zeta_r_list)

    def c(zeta):
        
        return np.sqrt(1/((1/(beta_1**2))-((zeta**2)/(h**2))))
    
    print(zeta_mode)
    print(f_mode)
    print(f)
    print(zeta_r)

    for m in range (4):

        zeta_mode_new = []
        f_mode_new = []

        for k in range (6):

            zeta_mode_new.append(zeta_r[k][m])
            f_mode_new.append(f[k])

        zeta_mode.append((zeta_mode_new))
        f_mode.append((f_mode_new))

if __name__ == '__main__':
    main()