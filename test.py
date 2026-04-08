import numpy as np

pts = np.loadtxt("NACA0012_processed.gnu")
x_cut = -0.520

upper = pts[:499]
lower = pts[499:]

idx_upper = np.argmax(upper[:, 0] < x_cut)
idx_lower = np.argmax(lower[:, 0] > x_cut)

upper_cut = upper[idx_upper:]    # 以 blunt_upper 开头
lower_cut = lower[:idx_lower]    # 以 blunt_lower 结尾

# 直接合并，不添加任何额外点
# Basilisk 自动将 lower_cut[-1] 连回 upper_cut[0]，这就是钝尾缘竖线段
new_pts = np.vstack([upper_cut, lower_cut])

new_pts[:, 0] += 1.0

np.savetxt("NACA0012_blunt_new.gnu", new_pts, fmt="%.6f")

blunt_upper = upper_cut[0]
print(f"TE thickness: {2*blunt_upper[1]:.5f} ({2*blunt_upper[1]/0.00488:.1f} Delta)")
print(f"Points: {len(new_pts)}")
print(f"Start: {new_pts[0]}")
print(f"End:   {new_pts[-1]}")