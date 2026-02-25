import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data_oscilador")

x = data[:,0]

for i in range(1,7):
    plt.plot(x, data[:,i], label=f"n={i-1}")

plt.xlabel("x")
plt.ylabel(r"$|\psi_n(x)|^2$")
plt.legend()
plt.grid()
plt.show()