from tkinter import Tk, Canvas

root = Tk()
cnv = Canvas(root, width=400, height=400)
cnv.pack()


position_x = [10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10]


cible = cnv.create_oval(100, 100, 200, 200, fill='blue')
cible2 = cnv.create_oval(300, 300, 400, 400, fill='blue')

def move(x, y):
    cnv.coords(cible, 100+x, 100, 200+x, 200)
    cnv.coords(cible2, 300+x, 300, 400+x, 400)

for i, pos in enumerate(position_x):
    root.after(500*i, move, pos, 0)

root.mainloop()