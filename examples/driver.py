import numpy as np

from my_python_package.operators import (

        newton_raphson,
        )

def main():

    i = input("Input a function of x: ")
    p = input("Input the first derivative of the function: ")
    x0 = float(input("Input an initial guess of the root: "))
    fu = eval("lambda x: " + i)
    fp = eval("lambda x: " + p)
    

    eps_a, xi, it = newton_raphson(x0, fu, fp, eps_s=1e-8)

    print("Value of the root: " , xi)
    print("Final error: " , eps_a)
    print("Iterations taken: " , it)

if __name__ == '__main__':
    main()