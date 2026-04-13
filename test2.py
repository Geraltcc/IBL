import numpy as np

pts = np.loadtxt("airfoil/NACA0012_processed.gnu")

pts[:, 0] += 1.0

np.savetxt("airfoil/NACA0012_processed.gnu", pts, fmt="%.6f")
