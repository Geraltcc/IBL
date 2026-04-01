原始待求解变量 $(\theta, H^*, C_\tau)$ ，转换为 $(\theta, H, C_\tau)$
已知关系：
$$
H^* = f(H_k, Re_\theta), \quad H_k = f(H, M_e), \quad Re_\theta = f(\theta, \rho_e, u_e, \mu_e)
$$
得
$$
\frac{dH^*}{d\xi} = \frac{\partial H^*}{\partial H_k}\left(\frac{\partial H_k}{\partial H}\frac{dH}{d\xi} + \frac{\partial H_k}{\partial M_e}\frac{d M_e}{d \xi}\right) + \frac{\partial H^*}{\partial Re_\theta}\left(\frac{d Re_\theta}{d\xi}\right)
$$
其中 
$$
\frac{d Re_\theta}{d\xi} = \theta \frac{d}{d\xi}\left(\frac{\rho_e u_e}{\mu_e}\right)+\frac{Re_\theta}{\theta}\frac{d\theta}{d\xi}
$$代入方程(2)得
$$
\theta\frac{\partial H^*}{\partial H_k}\left(\frac{\partial H_k}{\partial H}\frac{dH}{d\xi} + \frac{\partial H_k}{\partial M_e}\frac{dM_e}{d\xi}\right) + \theta \frac{\partial H^*}{\partial Re_\theta}\left(\theta\frac{d}{d\xi}\left(\frac{\rho_e u_e}{\mu_e}\right) + \frac{Re_\theta}{\theta}\frac{d\theta}{d\xi}\right)
=
2C_D - H^* \frac{C_f}{2} - [2H^{**} + H^*(1 - H)]\frac{\theta}{u_e}\frac{du_e}{d\xi}
$$
整理得
$$
\theta \frac{\partial H^*}{\partial H_k}\frac{\partial H_k}{\partial H}\frac{dH}{d\xi} + \frac{\partial H^*}{\partial Re_\theta}Re_\theta\frac{d\theta}{d\xi} = -\theta \frac{\partial H^*}{\partial H_k}\frac{\partial H_k}{\partial M_e}\frac{d M_e}{d\xi} - \theta^2\frac{\partial H^*}{\partial Re_\theta}\frac{d}{d\xi}\left(\frac{\rho_e u_e}{\mu_e}\right) + 2C_D -  H^* \frac{C_f}{2} -  [2H^{**} + H^*(1 - H)]\frac{\theta}{u_e}\frac{du_e}{d\xi}
$$
简化方程得 
$$f_1 \cdot \frac{dH}{d\xi} + f_2 \cdot \frac{d\theta}{d\xi} = RHS$$
代入方程(1)消去$\frac{d\theta}{d\xi}$得
$$
f_1\cdot \frac{dH}{d\xi} + f_2\cdot \left[\frac{C_f}{2} - (2 + H - M_e^2)\frac{\theta}{u_e}\frac{d u_e}{d\xi}\right] = RHS
$$
得到关于 $\frac{d H}{d\xi}$ 的方程
$$
f_1\cdot \frac{dH}{d\xi} = RHS - f_2 \cdot \left[\frac{C_f}{2} - (2 + H - M_e^2)\frac{\theta}{u_e}\frac{d u_e}{d\xi}\right]
$$
计算 $f_1$ 和 $f_2$ 的偏导
(1) $\frac{\partial H_k}{\partial H}$ 和 $\frac{\partial H_k}{\partial M_e}$
$$
\frac{\partial H_k}{\partial H} = \frac{1}{1 + 0.113M_e^2}
$$
$$
\frac{\partial H_k}{\partial M_e} = -\frac{M_e(0.590+0.226H)}{(1+0.113M_e^2)^2}
$$
(2) $\frac{\partial H^*}{\partial H_k}$ 和 $\frac{\partial H^*}{\partial Re_\theta}$
$$
\frac{\partial H^*}{\partial H_k} = 
\begin{cases}
-\left(0.165-\frac{1.6}{\sqrt{Re_\theta}}\right)\frac{(H_0-H_k)^{0.6}(H_0 + 0.6H_k)}{H_k^2}, & H_k < H_0 \\
2(H_k - H_0)\left[\frac{0.04}{H_k} + \frac{0.007 log(Re_\theta)}{(H_k-H_0+\frac{4}{log(Re_\theta)})^2}\right] + (H_k-H_0)^2\left[-\frac{0.04}{H_k^2}-\frac{0.014log(Re_\theta)}{(H_k-H_0+\frac{4}{log(Re_\theta)})^3}\right], & H_k > H_0
\end{cases}
$$

$$
\frac{d H_0}{d Re_\theta} = 
\begin{cases}
0, & Re_\theta < 400 \\
-\frac{400}{Re_\theta^2}, & Re_\theta > 400
\end{cases}
$$
$$
\frac{\partial H^*}{\partial Re_\theta} = 
\begin{cases}
-\frac{4}{Re_\theta^2} + 0.8Re_\theta^{-\frac{3}{2}}\frac{(H_0 - H_k)^{1.6}}{H_k}
, & H_k < H_0 \\
-\frac{4}{Re_\theta^2}-2(H_k-H_0)\frac{dH_0}{dRe_\theta}\left(\frac{0.04}{H_k} +\frac{0.007L}{(H_k-H_0+\frac{4}{log(Re_\theta)})^2}\right) + 0.007(H_k-H_0)^2\left[\frac{1}{Re_\theta(H_k-H_0+\frac{4}{log(Re_\theta)})^2} + 2log(Re_\theta)\frac{1}{(H_k-H_0)^3}\left(\frac{dH_0}{dRe_\theta} + \frac{4}{(log(Re_\theta))^2 Re_\theta}\right)\right]
, & H_k > H_0
\end{cases}
$$
