import numpy as np

pts = np.loadtxt("NACA0012_processed.gnu")

pts[:, 0] += 1.0

np.savetxt("NACA0012_blunt_new.gnu", pts, fmt="%.6f")
