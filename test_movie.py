from tkinter import Tk, Canvas
import numpy as np
import pandas as pd


FACTOR_DIST = 200
DELTA_POS = 10

# import dat of MD
data = pd.read_csv("RES.dat", sep=" ")
# Augmentation de la distance entre les atomes.
# ––––––––––––––––––––––––––––––––––––––––––––
x_b = np.array(data["x_b"])*FACTOR_DIST
x_c = np.array(data["x_c"])*FACTOR_DIST
y_c = np.array(data["y_c"])*FACTOR_DIST
# Barycenter
g_x = np.array(data["g_x"])*FACTOR_DIST
g_y = np.array(data["g_y"])*FACTOR_DIST

# Augmentation du delta de la différences entre la position
# d'origine et la position actuelle.
diff_xb = np.array([x_b[0] - x_b[i] for i in range(len(x_b))])
xb_modif = (x_b - diff_xb) + (diff_xb * DELTA_POS)
diff_xc = np.array([x_c[0] - x_c[i] for i in range(len(x_c))])
xc_modif = (x_c - diff_xc) + (diff_xc * DELTA_POS)
diff_yc = np.array([y_c[0] - y_c[i] for i in range(len(y_c))])
yc_modif = (y_c - diff_yc) + (diff_yc * DELTA_POS)
# Barycenter
diff_gx = np.array([g_x[0] - g_x[i] for i in range(len(g_x))])
gx_modif = (g_x - diff_gx) + (diff_gx * DELTA_POS)
diff_gy = np.array([g_y[0] - g_y[i] for i in range(len(g_y))])
gy_modif = (g_y - diff_gy) + (diff_gy * DELTA_POS)

# création de l'objet tkinter
root = Tk()
# création du canevas
cnv = Canvas(root, width=500, height=500)
cnv.pack()

# Initialisation des objets
dist_ab = cnv.create_line(75, 75, 100, 75, fill='green')
dist_ac = cnv.create_line(75, 75, 100, 100, fill='green')
dist_bc = cnv.create_line(300, 75, 200, 350, fill='green')
atom_a = cnv.create_oval(50, 50, 100, 100, fill='blue')
atom_b = cnv.create_oval(50, 50, 100, 100, fill='blue')
atom_c = cnv.create_oval(50, 50, 100, 100, fill="blue")
barycenter = cnv.create_oval(50, 50, 100, 100, fill="red")

def move(x_b, x_c, y_c, g_x, g_y):
    cnv.coords(dist_ab, 75, 75, 100+x_b, 75)
    cnv.coords(dist_ac, 75, 75, 75+x_c, 100+y_c)
    cnv.coords(dist_bc, 75+x_b, 75, 75+x_c, 100+y_c)
    
    cnv.coords(atom_b, 50+x_b, 50, 100+x_b, 100)
    cnv.coords(atom_c, 50+x_c, 50+y_c, 100+x_c, 100+y_c)
    cnv.coords(barycenter, 70+g_x, 70+g_y, 80+g_x, 80+g_y)

for i in range(len(x_b)):
    root.after(50*i, move, xb_modif[i], xc_modif[i],
               yc_modif[i], gx_modif[i], gy_modif[i])

root.mainloop()