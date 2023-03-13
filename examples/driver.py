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
    f_list_1 = np.array([0.1,0.2,0.5,1,2,5,10])
    f_list_2 = np.array([0.5,1,2,5,10])
    f_list_3 = np.array([1,2,5,10])
    zeta_0 = 1.3
    mode_1_c = np.zeros(0)
    mode_2_c = np.zeros(0)
    mode_3_c = np.zeros(0)
    zeta_max = (h**2) * ((beta_1**-2) - (beta_2**-2))

    def g(zeta):
        
        return (np.tan(2*np.pi*f*zeta)) - ((rho_2/rho_1)*((np.sqrt((zeta_max) - zeta**2))/zeta))
    
    def dgdx(zeta):

        return ((2*np.pi*f)*((1/np.cos(2*np.pi*f*zeta))**2)) - ((rho_2/rho_1)*zeta_max/(beta_1**2 * beta_2**2 * zeta**2 * np.sqrt((zeta_max)) - zeta**2))
    
    def c(f):
        
        return np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    
    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    f = 0.1

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_1_c = np.append(mode_1_c,cl)

    f = 0.2
    zeta_0 = 0.9

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_1_c = np.append(mode_1_c,cl)

    f = 0.5
    zeta_0 = 0.43

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 1.27

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_2_c = np.append(mode_2_c,cl)

    f = 1
    zeta_0 = 0.23

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.69

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 1.15

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_3_c = np.append(mode_3_c,cl)

    f = 2
    zeta_0 = 0.36

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.60

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 0.84

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_3_c = np.append(mode_3_c,cl)

    f = 5
    zeta_0 = 0.049

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.148

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 0.246
    print("hello!")

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_3_c = np.append(mode_3_c,cl)

    f = 10
    zeta_0 = 0.0248

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_1_c = np.append(mode_1_c,cl)

    zeta_0 = 0.0745

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_2_c = np.append(mode_2_c,cl)

    zeta_0 = 0.124

    eps_a_arr, zeta_r, it = newton_raphson(zeta_0, g, dgdx, eps_s = 1e-8)

    cl = np.sqrt(1/((1/(beta_1**2))-((zeta_r**2)/(h**2))))
    mode_3_c = np.append(mode_3_c,cl)

    print("Mode 1 Cs: ", mode_1_c)
    print("mode 2 Cs: ", mode_2_c)
    print("mode 3 Cs: ", mode_3_c)
    print(c(0.1))
    print(c(0.2))
    print(c(0.5))

    plt.plot(f_list_1,mode_1_c,'k') 
    plt.plot(f_list_2,mode_2_c, 'r')
    plt.plot(f_list_3,mode_3_c, 'c')
    plt.grid()
    plt.show()

if __name__ == '__main__':
    main()