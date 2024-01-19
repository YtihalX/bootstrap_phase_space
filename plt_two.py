import matplotlib.pyplot as plt
import pandas as pd

data0 = pd.read_csv("data0.csv")
data1 = pd.read_csv("data1.csv")
plt.figure(dpi = 300)
plt.rc('text', usetex = True)
plt.scatter(data0['E'], data0['xsq'], s = 0.2, marker = '.', edgecolor = None, label = r'$e^{mx}(k = 9)$')
plt.scatter(data1['E'], data1['xsq'], s = 0.2, marker = '.', edgecolor = None, label = r'$e^{mx} p^n(k = 4)$')
plt.xlabel(r'$E$')
plt.ylabel(r'$e^x$')
plt.legend()
# Get the y-axis object
# y_axis = plt.gca().axes.get_yaxis()

# Hide the y-axis
# y_axis.set_visible(False)

plt.savefig("plot_two.png")
