"""Simple MD simulation of triatomic model.

Usage:
======

python triatomic_model.py file [-h] [-l] [-n NB_ITER] [-g] [-s]


positional arguments:
  file                  File name of MD result

optional arguments:
  -h, --help            show this help message and exit
  -l, --launch          Launch MD
  -n NB_ITER, --nb_iter NB_ITER
                        Number of iteration for MD (default: 10000)
  -g, --graph           Show graph
  -s, --save            Save graph
"""

__authors__ = ("William Amory", "Lucas Rouaud")
__date__ = "2022-04-03"

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os, sys
# Module for get arguments, options and path file.
import argparse, pathlib


# Define initial position value of atoms
X_B0 = 1.23 # in A
X_C0 = X_B0/2# in A
Y_C0 = (X_B0**2 + X_C0**2)**0.5 #  in A

# Define initial euclidean distance.
# initial euclidean distance of A-B = X_B0
l_AC0 = (X_C0**2 + Y_C0**2)**0.5
l_BC0 = ((X_C0 - X_B0)**2 + Y_C0**2)**0.5

# Define constants for energy calculation.
K = 620.0 # in kcal/mol/A^2
M_C = 12.011e-3 # in kg/mol (1 amu = 1 g/mol)
M_O = 15.999e-3 # in kg/mol (1 amu = 1 g/mol)
# Define constants for performing MD.
BOLTZMANN_CST = 1.9872041e-3 # in kcal/mol/K
N_DOF = 3 # nb of deg of freedom (only 1 -> x)
######################################
# Define constants you can change ;-).
######################################
DELTA_T_FORCE = 0.0001 # 
DELTA_t = 0.0005 # time step in ps
INIT_T = 300 # T in K
######## No Change below!!!###########


def get_arg_and_option():
    """Get arguments and options.

    Returns
    -------
    class 'argparse.Namespace' : args
        arguments and options used.
    """
    parser = argparse.ArgumentParser(
        description='Simple MD simulation of triatomic model.')
    parser.add_argument('file', type=pathlib.Path,
                        help="File name of MD result")
    parser.add_argument("-l", "--launch", help="Launch MD", action="store_true")
    parser.add_argument("-n", "--nb_iter", help="Number of iteration for MD (default: %(default)s)", type=int, default=10000,)
    parser.add_argument("-g", "--graph", help="Show graph", action="store_true")
    parser.add_argument("-s", "--save", help="Save graph\n", action="store_true")
    return parser.parse_args()


def calc_ene(x_b, x_c, y_c):
    """Calculate potential enregy.

    Parameters
    ----------
    x_b : float
        x position of atom B.
    x_c : float
        x position of atom C.
    y_c : float
        y position of atom C.

    Returns
    -------
    float : ene
        potential energy of model.
    """
    # Define current euclidean distance.
    # Current euclidean distance of A-B = x_b
    l_AC = (x_c**2 + y_c**2)**0.5
    l_BC = ((x_c  - x_b)**2  + y_c**2)**0.5
    
    ene = K *  ((x_b - X_B0)**2 + (l_AC - l_AC0)**2 + (l_BC - l_BC0)**2)
    return ene


def calc_force(x_b, x_c, y_c):
    """Calculate force.

    Parameters
    ----------
    x_b : float
        x position of atom B.
    x_c : float
        x position of atom C.
    y_c : float
        y position of atom C.

    Returns
    -------
    float : f_xb
        Force applied on atom B's x position.
    float : f_xc
        Force applied on atom C's x position.
    float : f_yc
        Force applied on atom C's y position.
    """
    # Force is the derivative of Epot, unit of force is kcal/mol/A.
    # Here we have an numerical derivative of Epot.
    f_xb = -( (calc_ene(x_b + DELTA_T_FORCE, x_c, y_c)
            -calc_ene(x_b - DELTA_T_FORCE, x_c, y_c))
            /(2*DELTA_T_FORCE) )
    f_xc = -( (calc_ene(x_b, x_c + DELTA_T_FORCE, y_c)
            -calc_ene(x_b, x_c - DELTA_T_FORCE, y_c))
            /(2*DELTA_T_FORCE) )
    f_yc = -( (calc_ene(x_b, x_c, y_c + DELTA_T_FORCE)
            -calc_ene(x_b, x_c, y_c - DELTA_T_FORCE))
            /(2*DELTA_T_FORCE) )
    return f_xb, f_xc, f_yc


def calc_acc(f_xb, f_xc, f_yc):
    """Calculate acceleration.

    Parameters
    ----------
    f_xb : float
        Force applied on atom B's x position.
    f_xc : float
        Force applied on atom c's x position.
    f_yc : float
        Force applied on atom C's y position.

    Returns
    -------
    float : acc_xb
        Acceleration applied on atom B's x position.
    float : acc_xc
        Acceleration applied on atom c's x position.
    float : acc_yc
        Acceleration applied on atom C's y position.
    """
    acc_xb = (f_xb / M_C) * 0.418
    acc_xc = (f_xc / M_C) * 0.418
    acc_yc = (f_yc / M_C) * 0.418
    return acc_xb, acc_xc, acc_yc


def start_Verlet(x_b, x_c, y_c):
    """Start MD with Verlet algorithm.

    Parameters
    ----------
    x_b : float
        x position of atom B.
    x_c : float
        x position of atom C.
    y_c : float
        y position of atom C.

    Returns
    -------
    float : xb_prev
        Previous x position of atom B.
    float : xc_prev
        Previous x position of atom C.
    float : yc_prev
        Previous y position of atom C.
    """
    # Choose a velocity (for our C atom) which corresponds to T=300K.
    # v calculated from: v^2 = kT/m.
    # Factor 0.418 to convert m/s to A/ps.
    v2 = ((INIT_T*BOLTZMANN_CST)/M_C)*0.418
    v = v2**.5
    # Compute inital acc.
    # (should be 0 if spring is at rest ;-) )
    f_xb, f_xc, f_yc = calc_force(x_b, x_c, y_c)
    acc_xb, acc_xc, acc_yc = calc_acc(f_xb, f_xc, f_yc)
    # Initialize Verlet by determining previous position (at t0-dt).
    # Use Taylor expansion truncated after third term --> error O(dt**3):
    # r(t0-dt) = r(t0) - v(t0)*dt + 0.5*a(t0)*dt**2
    xb_prev = x_b - v*DELTA_t + 0.5*acc_xb*DELTA_t**2
    xc_prev = x_c - v*DELTA_t + 0.5*acc_xc*DELTA_t**2
    yc_prev = y_c - v*DELTA_t + 0.5*acc_yc*DELTA_t**2
    return xb_prev, xc_prev, yc_prev


def do_MD(x_b, x_c, y_c, nsteps):
    """Launch MD of triatomic model.

    Parameters
    ----------
    x_b : float
        x position of atom B.
    x_c : float
        x position of atom C.
    y_c : float
        y position of atom C.
    """
    with args.file.open("w") as f_out:
        out = ("step t x_b x_c y_c ene_pot ene_kin ene_tot v_xb v_xc v_yc "+
               "f_xb f_xc f_yc T acc_xb acc_xc acc_yc\n")
        f_out.write(out)
        for step in range(nsteps):
            t = step * DELTA_t
            if step == 0:
                xb_prev, xc_prev, yc_prev = start_Verlet(x_b, x_c, y_c)
            # Calc force.
            f_xb, f_xc, f_yc = calc_force(x_b, x_c, y_c)
            # Compute acceleration (f in kcal/mol/A, M_C in kg/mol -> use conversion
            #                       factor to get acceleration in A/ps^2).
            acc_xb, acc_xc, acc_yc = calc_acc(f_xb, f_xc, f_yc)
            # Compute new pos with Verlet algorithm.
            xb_new = (2 * x_b) - xb_prev + (DELTA_t**2 * acc_xb)
            xc_new = (2 * x_c) - xc_prev + (DELTA_t**2 * acc_xc)
            yc_new = (2 * y_c) - yc_prev + (DELTA_t**2 * acc_yc)
            # Calc v (unit is A/ps).
            v_xb = (xb_new - xb_prev) / (2 * DELTA_t)
            v_xc = (xc_new - xc_prev) / (2 * DELTA_t)
            v_yc = (yc_new - yc_prev) / (2 * DELTA_t)
            # Calc Ekin.
            # M in kg/mol, v in m/s (1 A/ps = 100 m/s) ===> Ekin is in J/mol
            ene_kin = ( (.5 * M_C * (v_xb*100.0)**2) +
                       (.5 * M_C * (v_xc*100.0)**2) +
                       (.5 * M_C * (v_yc*100.0)**2) )
            # Convert Ekin to kcal/mol.
            ene_kin /= (1000.0 * 4.18)
            # calc Temperature (Ekin in kcal/mol, BOLTZMANN_CST in kcal/mol/K)
            T = (2 * ene_kin) / (N_DOF * BOLTZMANN_CST)
            # Calc Epot & Etot.
            ene_pot = calc_ene(x_b, x_c, y_c)
            ene_tot = ene_pot + ene_kin
            # Monitor properties here!
            out = (f"{step} {t:.4f} {x_b:.8f} {x_c:.8f} {y_c:.8f} {ene_pot:.8f} "+
                   f"{ene_kin:.8f} {ene_tot:.8f} {v_xb:.8f} {v_xc:.8f} "+
                   f"{v_yc:.8f} {f_xb:.8f} {f_xc:.8f} {f_yc:.8f} {T:.8f} "+
                   f"{acc_xb:.8f} {acc_xc:.8f} {acc_yc:.8f}\n")
            f_out.write(out)
           	# Update new and prev pos for the next iteration.
            xb_prev = x_b
            xc_prev = x_c
            yc_prev = y_c
            
            x_b = xb_new
            x_c = xc_new
            y_c = yc_new


def print_graph(filin):
    if not os.path.exists(filin):
        sys.exit("\nError !!!!\n\nRES.dat don't exist\n\n ---\n "+
                 "lauch MD with : python3 triatomic_model.py -l\n ---\n\n")
    data = pd.read_csv(filin, sep=" ")

    # reponse = input("entrer vos options : ")
    # plt.plot(data[reponse])
    # plt.show()

    reponse = ""
    while reponse != 'q' and reponse != 'quit':
        reponse = input("Choose variable to plot\n\n"+
                        "step : \n"+
                        "t : \n"+
                        "x_b : \n"+
                        "x_c : \n"+
                        "y_c : \n"+
                        "ene_pot : potential energy\n"+
                        "ene_kin : kinetic energy\n"+
                        "ene_tot : total energy\n"+
                        "v_xb : \n"+
                        "v_xc : \n"+
                        "v_yc : \n"+
                        "f_xb : \n"+
                        "f_xc : \n"+
                        "f_yc : \n"+
                        "T : Temperature\n"+
                        "acc_xb : \n"+
                        "acc_xc : \n"+
                        "acc_yc : \n\n"
                        +"entrer vos options : ")
        if reponse == 'q' or reponse == 'quit':
            break
        reponse = reponse.split()
        for rep in reponse:
            plt.plot(data[rep])
        plt.xlim(0, 1000)
        plt.show() 


########
# MAIN #
########
if __name__ == "__main__":
    
    # Get file name and option.
    args = get_arg_and_option()
    
    # Quit script and print message if no options are specified.
    if args.launch == False and args.graph == False and args.save == False:
        sys.exit(("\nError !!!!\n\nno option specified (used one or more of "+
                  "the following options: -l -g -s)\n\n ---\n view option "+
                  "with : python3 triatomic_model.py -h\n ---\n\n"))
    
    # Lauch MD if launch option is specified.
    if args.launch == True:
        # Define starting coordinates.
        # xa = 0.0 ==> centered at origin (static point).
        # xb is x coor of B (moves along x axis = mobile point).
        # xc is x coor of C (moves along x axis = mobile point).
        # yc is x coor of C (moves along y axis = mobile point).
        #     ==> at starting position, the spring is at rest.
        x_b = X_B0
        x_c = X_C0
        y_c = Y_C0
        # Launch MD!
        do_MD(x_b, x_c, y_c, args.nb_iter)

    if args.graph == True or args.save == True:
        print_graph(args.file)
    
