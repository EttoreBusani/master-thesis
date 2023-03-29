import numpy as np
import matplotlib.pyplot as plt

# Load data from file
x, y = np.loadtxt('data.txt', unpack=True)

# Plot the data
plt.plot(x, y)
plt.show()
