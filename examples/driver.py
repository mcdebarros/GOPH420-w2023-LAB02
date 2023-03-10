import numpy as np

import matplotlib.pyplot as plt

from newton_raphson.rootfinder import newton_raphson

def main():

    beta_1 = 1900
    beta_2 = 3200
    rho_1 = 1800
    rho_2 = 2500
    h = 4000
    f = 1
    zeta_0 = 1
    zeta_max = (h**2) * ((beta_1**-2) - (beta_2**-2))

    def g(zeta):
        
        return (np.tan(2*np.pi*f*zeta)) - ((rho_2/rho_1)*((np.sqrt((h**2 * (beta_1**-2 - beta_2**-2)) - zeta**2))/zeta))
    
    def dgdx(zeta):

        return ((2*np.pi*f)*((1/np.cos(2*np.pi*f*zeta))**2)) - ((rho_2/rho_1)*zeta_max/(beta_1**2 * beta_2**2 * zeta**2 * np.sqrt((h**2 * (beta_1**-2 - beta_2**-2)) - zeta**2)))
    
    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    print("Value of the root: " , zeta_r)
    print("Final error: " , eps_a_arr)
    print("Iterations taken: " , it)

    zeta_plot = np.linspace(-zeta_max,(zeta_max),100)
    plt.plot(zeta_plot,g(zeta_plot),'r')
    plt.plot(zeta_r,g(zeta_r),'k')
    plt.plot(zeta_0,g(zeta_0),'c')
    plt.plot()
    plt.grid()
    plt.show()

if __name__ == '__main__':
    main()