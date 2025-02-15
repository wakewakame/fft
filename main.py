import numpy as np
import matplotlib.pyplot as plt

N = 8
x = np.sin(np.linspace(0, 2*np.pi, N+1)[:-1])
X = np.fft.fft(x)
np.set_printoptions(suppress=True, formatter={'float': '{:.03f}'.format})
print(np.abs(X))
plt.plot(np.linspace(0, N-1, N), np.abs(X))
plt.show()
