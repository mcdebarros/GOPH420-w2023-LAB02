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
    zeta_max = np.sqrt((h**2) * ((beta_1**-2) - (beta_2**-2)))
    f = []
    zeta_r = []
    f_mode = []
    zeta_mode = []
    cl_mode = []
    lambda_mode = []
    eps_n_list = []

    for n in range (6):

        f_new = (((2*n) + 1)/(4*zeta_max))
        f.append(f_new)
        zeta_r_list = []
        eps_it_list = []

        for k in range (n+1):

            zeta_0 = ((((2*k) + 1) / (4*f_new))) - (zeta_max * 1e-5)

            def g(zeta):
        
                return ((rho_2/rho_1)*((np.sqrt((zeta_max**2)-(zeta**2)))/zeta)) - (np.tan(2*(np.pi)*f_new*zeta))
    
            def dgdx(zeta):

                return -((rho_2*(zeta_max**2))/(rho_1*(zeta**2)*np.sqrt((zeta_max**2)-(zeta**2)))) - ((2*np.pi*f_new)*((1/np.cos(2*np.pi*f_new*zeta))**2))
            
            zeta_r_new, it, eps_a_arr, eps_a = newton_raphson(zeta_0, g, dgdx, eps_s)
            zeta_r_list.append(zeta_r_new)
            eps_it_list.append(eps_a_arr)
        
        zeta_r.append(zeta_r_list)
        eps_n_list.append(eps_it_list)

    for m in range (4):

        zeta_mode_new = []
        f_mode_new = []

        for k in range (m,6):

            zeta_mode_new.append(zeta_r[k][m])
            f_mode_new.append(f[k])

        zeta_mode_new = np.array(zeta_mode_new)
        f_mode_new = np.array(f_mode_new)
        cl_mode_new = ((beta_1**-2) - ((zeta_mode_new**2)/(h**2))) ** -0.5
        lambda_mode_new = cl_mode_new/f_mode_new
        zeta_mode.append(zeta_mode_new)
        f_mode.append(f_mode_new)
        cl_mode.append(cl_mode_new)
        lambda_mode.append(lambda_mode_new)

    plt.figure(figsize = (6,8))

    plt.subplot(2,1,1)
    for f,cl in zip(f_mode,cl_mode):
        plt.plot(f,cl)
    plt.xlabel('f [Hz]')
    plt.ylabel('c_L [m/s]')
    plt.title('figure 1: c_L as a function of f', fontsize = 7)
    plt.legend([f'mode {k}' for k,_ in enumerate (f_mode)])

    plt.subplot(2,1,2)
    for f,lam in zip(f_mode,lambda_mode):
        plt.plot(f,lam)
    plt.xlabel('f [Hz]')
    plt.ylabel('lambda [m]')
    plt.title('figure 2: lambda_L as a function of f', fontsize = 7)
    plt.legend([f'mode {k}' for k,_ in enumerate (f_mode)])

    plt.savefig('examples/modes.png')

if __name__ == '__main__':
    main()