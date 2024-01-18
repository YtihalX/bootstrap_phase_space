import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("data.csv")
plt.figure(dpi = 300)
plt.rc('text', usetex = True)
plt.scatter(data['E'], data['xsq'], s = 0.2, marker = '.', edgecolor = None)
plt.xlabel(r'$E$')
plt.ylabel(r'$x^2$')
# Get the y-axis object
# y_axis = plt.gca().axes.get_yaxis()

# Hide the y-axis
# y_axis.set_visible(False)

plt.savefig("plot.png")
