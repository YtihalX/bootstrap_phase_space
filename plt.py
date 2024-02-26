import matplotlib.pyplot as plt
import pandas as pd

data8 = pd.read_csv("data8.csv")
data9 = pd.read_csv("data9.csv")
data15 = pd.read_csv("data15.csv")
data25 = pd.read_csv("data25.csv")
plt.figure(dpi = 300)
plt.rc('text', usetex = True)
plt.scatter(data8['E'], data8['xsq'], s = 0.8, marker = '.', edgecolor = None, label = r'$k = 8$')
plt.scatter(data9['E'], data9['xsq'], s = 0.8, marker = '.', edgecolor = None, label = r'$k = 9$')
plt.scatter(data15['E'], data15['xsq'], s = 0.8, marker = '.', edgecolor = None, label = r'$k = 15$')
plt.scatter(data25['E'], data25['xsq'], s = 0.8, marker = '.', edgecolor = None, label = r'$k = 25$')
plt.xlabel(r'$E$')
plt.ylabel(r'$x^2$')
plt.legend()
# Get the y-axis object
# y_axis = plt.gca().axes.get_yaxis()

# Hide the y-axis
# y_axis.set_visible(False)

plt.savefig("double_well_xp.png")
