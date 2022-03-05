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
  -g, --generate_graph  Generate graph of MD analysis
  -s, --show            Show graph
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
L_AC0 = (X_C0**2 + Y_C0**2)**0.5
L_BC0 = ((X_C0 - X_B0)**2 + Y_C0**2)**0.5

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
    parser.add_argument("-g", "--generate_graph", action="store_true",
                        help="Generate graph of MD")
    parser.add_argument("-s", "--show", help="Show graph\n", action="store_true")
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
    
    ene = K *  ((x_b - X_B0)**2 + (l_AC - L_AC0)**2 + (l_BC - L_BC0)**2)
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


def calculate_area(x_b, y_c):
    """_summary_

    Parameters
    ----------
    x_b : float
        Length side of the triangle.
    y_c : float
        Length of the triangle height.

    Returns
    -------
    float
        Area of the triangle
    """
    return ((x_b * y_c) / 2)


def calculate_barycenter(x_b, x_c, y_c):
    g_x = (x_b + x_c) / 3
    g_y = (y_c) / 3
    return g_x, g_y

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
        out = ("step time x_b x_c y_c ene_pot ene_kin ene_tot v_xb v_xc v_yc "+
               "f_xb f_xc f_yc T acc_xb acc_xc acc_yc area g_x g_y\n")
        f_out.write(out)
        for step in range(nsteps):
            time = step * DELTA_t
            if step == 0:
                xb_prev, xc_prev, yc_prev = start_Verlet(x_b, x_c, y_c)
            # Calculate force.
            f_xb, f_xc, f_yc = calc_force(x_b, x_c, y_c)
            # Compute acceleration (f in kcal/mol/A, M_C in kg/mol -> use
            # conversion factor to get acceleration in A/ps^2).
            acc_xb, acc_xc, acc_yc = calc_acc(f_xb, f_xc, f_yc)
            # Compute new pos with Verlet algorithm.
            xb_new = (2 * x_b) - xb_prev + (DELTA_t**2 * acc_xb)
            xc_new = (2 * x_c) - xc_prev + (DELTA_t**2 * acc_xc)
            yc_new = (2 * y_c) - yc_prev + (DELTA_t**2 * acc_yc)
            # Calc v (unit is A/ps).
            v_xb = (xb_new - xb_prev) / (2 * DELTA_t)
            v_xc = (xc_new - xc_prev) / (2 * DELTA_t)
            v_yc = (yc_new - yc_prev) / (2 * DELTA_t)
            # Calculate Ekin.
            # M in kg/mol, v in m/s (1 A/ps = 100 m/s) ===> Ekin is in J/mol.
            ene_kin = ( (.5 * M_C * (v_xb*100.0)**2) +
                       (.5 * M_C * (v_xc*100.0)**2) +
                       (.5 * M_C * (v_yc*100.0)**2) )
            # Convert Ekin to kcal/mol.
            ene_kin /= (1000.0 * 4.18)
            # calculate Temperature (Ekin in kcal/mol, BOLTZMANN_CST
            # in kcal/mol/K).
            T = (2 * ene_kin) / (N_DOF * BOLTZMANN_CST)
            # Calculate Epot & Etot.
            ene_pot = calc_ene(x_b, x_c, y_c)
            ene_tot = ene_pot + ene_kin
            # Calculate area of molecule in A^2.
            area = calculate_area(x_b, y_c)
            # Calculate barycenter
            g_x, g_y = calculate_barycenter(x_b, x_c, y_c)
            # Monitor properties here !
            out = (f"{step} {time:.4f} {x_b:.8f} {x_c:.8f} {y_c:.8f} "
                   f"{ene_pot:.8f} {ene_kin:.8f} {ene_tot:.8f} {v_xb:.8f} "+
                   f"{v_xc:.8f} {v_yc:.8f} {f_xb:.8f} {f_xc:.8f} {f_yc:.8f} "+
                   f"{T:.8f} {acc_xb:.8f} {acc_xc:.8f} {acc_yc:.8f} "+
                   f"{area:.8f} {g_x:.8f} {g_y:.8f}\n")
            f_out.write(out)
           	# Update new and prev pos for the next iteration.
            xb_prev = x_b
            xc_prev = x_c
            yc_prev = y_c
            
            x_b = xb_new
            x_c = xc_new
            y_c = yc_new

def user_input_graph():
    """Print user choices variables to plot.

    Returns
    -------
    str : answer
        User choices of variables to plot.
    """
    answer = input("\nChoose variable to plot\n\n"+
                        "Variable  Description\n"+
                        "––––––––  –––––––––––\n"+
                        "step :\t  Number of itération\n"+
                        "t :\t  Time of MD in ps\n"+
                        "x_b :\t  x position of atom B\n"+
                        "x_c :\t  x position of atom C\n"+
                        "y_c :\t  y position of atom C\n"+
                        "ene_pot : Potential energy\n"+
                        "ene_kin : Kinetic energy\n"+
                        "ene_tot : Total energy\n"+
                        "v_xb :\t  Velocity on atom B's x position \n"+
                        "v_xc :\t  Velocity on atom C's x position\n"+
                        "v_yc :\t  Velocity on atom C's y position\n"+
                        "f_xb :\t  Force on atom B's x position\n"+
                        "f_xc :\t  Force on atom C's x position\n"+
                        "f_yc :\t  Force on atom C's y position\n"+
                        "T :\t  Temperature\n"+
                        "acc_xb :  Acceleration on atom B's x position\n"+
                        "acc_xc :  Acceleration on atom C's x position\n"+
                        "acc_yc :  Acceleration on atom C's y position\n\n"+
                        "To quit enter q or quit\n\n"+
                        "Enter the variables to represent : ")
    return answer


def user_graph(filin):
    """User generate graph function of MD résults.

    Parameters
    ----------
    filin : str
        file name.
    """
    # Exit and print error message if filin don't exist.
    if not os.path.exists(filin):
        sys.exit(f"\nError !!!!\n\n{filin} don't exist\n\n ---\n "+
                 "Launch MD with : python3 triatomic_model.py -l\n"+
                 "Or\n View help with : python3 triatomic_model.py -h\n"+
                 "---\n\n")
    data = pd.read_csv(filin, sep=" ")
    answer = user_input_graph()
    while answer.lower() != 'q' and answer.lower() != 'quit':
        answer = answer.split()
        sns.set_style("whitegrid")
        plt.figure()
        for rep in answer:
            plt.plot(data["time"], data[rep])
        #plt.xlim(1000, 2000)
        plt.xlabel("Time in ps")
        plt.show()
        answer2 = input("new graph ? (y/n) : ")
        if answer2 == "n":
            break
        else:
            answer = input("Enter the variables to represent : ")
        

def generate_graph(filin):
    """User generate graph function of MD résults.

    Parameters
    ----------
    filin : str
        file name.
    """
    # Exit and print error message if filin don't exist.
    if not os.path.exists(filin):
        sys.exit(f"\nError !!!!\n\n{filin} don't exist\n\n ---\n "+
                 "Launch MD with : python3 triatomic_model.py -l\n"+
                 "Or\n View help with : python3 triatomic_model.py -h\n"+
                 "---\n\n")
    data = pd.read_csv(filin, sep=" ")
    sns.set_style("whitegrid")
    
    # Plot energy of molecule.
    plt.figure(figsize=(20, 10))
    plt.rc("font", size=15)
    plt.plot(data["time"], data["ene_tot"])
    plt.plot(data["time"], data["ene_pot"])
    plt.plot(data["time"], data["ene_kin"])
    plt.legend(["Total energy", "Potential energy", "Kinetic energy"])
    plt.title("Energy of the molecule")
    plt.xlabel("Time of MD (in ps)")
    plt.ylabel("Energy (in kcal/mol)")
    plt.xlim(0, 0.5)
    plt.savefig("energy.png")
    
    # Plot temperature of molecule.
    plt.figure(figsize=(20, 10))
    plt.rc("font", size=15)
    plt.plot(data["time"], data["T"])
    plt.title("Temperature of the molecule")
    plt.xlabel("Time of MD (in ps)")
    plt.ylabel("Temperature (in K)")
    plt.xlim(0, 0.5)
    plt.savefig("temperature.png")
    
    # Plot force applied of atoms.
    plt.figure(figsize=(20, 10))
    plt.rc("font", size=15)
    plt.plot(data["time"], data["f_xb"])
    plt.plot(data["time"], data["f_xc"])
    plt.plot(data["time"], data["f_yc"])
    plt.title("Force applied to atoms of the molecule")
    plt.xlabel("Time of MD (in ps)")
    plt.ylabel("Force (in kcal/mol/A)")
    plt.legend(["Force applied on atom B's x position",
                "Force applied on atom C's x position",
                "Force applied on atom C's y position"])
    plt.xlim(0, 0.5)
    plt.savefig("force.png")
    
    # Plot acceleration of atoms.
    plt.figure(figsize=(20, 10))
    plt.rc("font", size=15)
    plt.plot(data["time"], data["acc_xb"])
    plt.plot(data["time"], data["acc_xc"])
    plt.plot(data["time"], data["acc_yc"])
    plt.title("acceleration applied to atoms of the molecule")
    plt.xlabel("Time of MD (in ps)")
    plt.ylabel("Acceleration (in A/ps^2)")
    plt.legend(["Acceleration applied on atom B's x position",
                "Acceleration applied on atom C's x position",
                "Acceleration applied on atom C's y position"])
    plt.xlim(0, 0.5)
    plt.savefig("acceleration.png")
    
    # Plot Velocity of atoms.
    plt.figure(figsize=(20, 10))
    plt.rc("font", size=15)
    plt.plot(data["time"], data["v_xb"])
    plt.plot(data["time"], data["v_xc"])
    plt.plot(data["time"], data["v_yc"])
    plt.title("Vitesse applied to atoms of the molecule")
    plt.xlabel("Time of MD (in ps)")
    plt.ylabel("Velocity (in m/s)")
    plt.legend(["Velocity applied on atom B's x position",
                "Velocity applied on atom C's x position",
                "Velocity applied on atom C's y position"])
    plt.xlim(0, 0.5)
    plt.savefig("velocity.png")
    
    # Plot position of atoms.
    plt.figure(figsize=(20, 10))
    plt.rc("font", size=15)
    plt.plot(data["time"], data["x_b"])
    plt.plot(data["time"], data["x_c"])
    plt.plot(data["time"], data["y_c"])
    plt.title("x position of atom B of atoms of the molecule")
    plt.xlabel("Time of MD (in ps)")
    plt.ylabel("Position")
    plt.legend(["x position of atom B", "x position of atom C",
                "y position of atom C"])
    plt.xlim(0, 0.5)
    plt.savefig("position.png")
    
    # Plot area of molecule
    plt.figure(figsize=(20, 10))
    plt.rc("font", size=15)
    plt.plot(data["time"], data["area"])
    plt.title("Area the molecule")
    plt.xlabel("Time of MD (in ps)")
    plt.ylabel("Area (in A^2)")
    plt.xlim(0, 0.5)
    plt.savefig("area.png")


########
# MAIN #
########
if __name__ == "__main__":
    
    # Get file name and option.
    args = get_arg_and_option()
    
    # Quit script and print message if no options are specified.
    if (args.launch == False and args.generate_graph == False and
        args.show == False):
        sys.exit(("\nError !!!!\n\nNo option specified (used one or more of "+
                  "the following options: -l -g -s)\n\n ---\n View option "+
                  "with : python3 triatomic_model.py -h\n ---\n\n"))
    
    # Lauch MD if launch option is specified.
    if args.launch == True:
        print("MD progress...")
        # Define starting coordinates.
        # xa = 0.0 ==> centered at origin (static point).
        # xb is x coor of B (moves along x axis = mobile point).
        # xc is x coor of C (moves along x axis = mobile point).
        # yc is y coor of C (moves along y axis = mobile point).
        #     ==> at starting position, the spring is at rest.
        x_b = X_B0
        x_c = X_C0
        y_c = Y_C0
        # Launch MD!
        do_MD(x_b, x_c, y_c, args.nb_iter)
        print("MD completed")

    if args.show == True:
        print("MD user generate graph progress...")
        user_graph(args.file)
        print("MD user generate graph completed")
    
    if args.generate_graph == True:
        print("MD generate graph progress...")
        generate_graph(args.file)
        print("MD generate graph completed")
    
    print("\n\nThank's see you soon ;)\n\u00A9 William Amory, "+
          "Lucas Rouaud - 2022\n\n")