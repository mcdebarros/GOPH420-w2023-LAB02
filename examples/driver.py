import numpy as np

import matplotlib.pyplot as plt

from newton_raphson.rootfinder import newton_raphson

def main():

    beta_1 = 1900
    beta_2 = 3200
    rho_1 = 1800
    rho_2 = 2500
    h = 4000
    f = 0.2
    f_list_1 = np.array([0.2,0.5,1,2,5,10,20,50,100])
    f_list_2 = np.array([0.5,1,2,5,10,20,50,100])
    f_list_3 = np.array([1,2,5,10,20,50,100])
    zeta_0 = 0.9
    mode_1_c = np.zeros(0)
    mode_2_c = np.zeros(0)
    mode_3_c = np.zeros(0)
    zeta_max = (h**2) * ((beta_1**-2) - (beta_2**-2))
    zeta_plot = np.linspace(0.005,zeta_max,1000)

    def g(zeta):
        
        return (np.tan(2*np.pi*f*zeta)) - ((rho_2/rho_1)*((np.sqrt((zeta_max) - (zeta**2)))/zeta))
    
    def dgdx(zeta):

        return ((2*np.pi*f)*((1/np.cos(2*np.pi*f*zeta))**2)) - ((rho_2/rho_1)*zeta_max/(beta_1**2 * beta_2**2 * zeta**2 * np.sqrt((zeta_max)) - zeta**2))
    
    def c(zeta):
        
        return np.sqrt(1/((1/(beta_1**2))-((zeta**2)/(h**2))))

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)
    zeta_r1 = zeta_r

    f = 0.5
    zeta_0 = 0.43
    zeta_01 = zeta_0

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)
    zeta_r1 = zeta_r

    zeta_0 = 1.27
    zeta_02 = zeta_0

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_2_c = np.append(mode_2_c,cl)
    zeta_r2 = zeta_r

    f = 1
    zeta_0 = 0.23
    zeta_01 = zeta_0

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)
    zeta_r1 = zeta_r

    zeta_0 = 0.69
    zeta_02 = zeta_0

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_2_c = np.append(mode_2_c,cl)
    zeta_r1 = zeta_r

    zeta_0 = 1.15
    zeta_03 = zeta_0

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_3_c = np.append(mode_3_c,cl)
    zeta_r3 = zeta_r

    f = 2
    zeta_0 = 0.12
    zeta_01 = zeta_0

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)
    zeta_r1 = zeta_r

    zeta_0 = 0.36
    zeta_02 = zeta_0

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_2_c = np.append(mode_2_c,cl)
    zeta_r2 = zeta_r

    zeta_0 = 0.60
    zeta_03 = zeta_0

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_3_c = np.append(mode_3_c,cl)
    zeta_r3 = zeta_r

    f = 5
    zeta_0 = 0.049

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.148

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 0.246

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_3_c = np.append(mode_3_c,cl)

    f = 10
    zeta_0 = 0.0248

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.0745

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 0.124

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_3_c = np.append(mode_3_c,cl)

    f = 20
    zeta_0 = 0.0124

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.037

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 0.06

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_3_c = np.append(mode_3_c,cl)

    f = 50
    zeta_0 = 0.005

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.015

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 0.025

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_3_c = np.append(mode_3_c,cl)

    f = 100
    zeta_0 = 0.002

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.007

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 0.012

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = c(zeta_r)
    mode_3_c = np.append(mode_3_c,cl)

    print("Mode 1 Cs: ", mode_1_c)
    print("mode 2 Cs: ", mode_2_c)
    print("mode 3 Cs: ", mode_3_c)

    plt.plot(f_list_1,mode_1_c,'k') 
    plt.plot(f_list_2,mode_2_c, 'r')
    plt.plot(f_list_3,mode_3_c, 'c')
    plt.grid()
    plt.show()

if __name__ == '__main__':
    main()