import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Cargar campo temporal guardado por Fortran
# ============================================================
data = np.loadtxt("nlse_normalizada.dat")

xi = data[:, 0]
tau = data[:, 1]
u_re = data[:, 2]
u_im = data[:, 3]

u = u_re + 1j*u_im

xi_vals = np.unique(xi)
tau_vals = np.unique(tau)

nxi = len(xi_vals)
nt = len(tau_vals)

U = u.reshape(nxi, nt)

dtau = tau_vals[1] - tau_vals[0]

# ============================================================
# Calcular espectro de intensidad para cada xi
# ============================================================

U_spec = np.fft.fft(U, axis=1)

# Normalización tipo transformada continua
U_spec = U_spec * dtau / np.sqrt(2*np.pi)

# Centrar frecuencias
U_spec = np.fft.fftshift(U_spec, axes=1)

omega = 2*np.pi*np.fft.fftfreq(nt, d=dtau)
omega = np.fft.fftshift(omega)

# Intensidad espectral
I_spec = np.abs(U_spec)**2

# Normalización global
I_spec = I_spec / np.max(I_spec)

# ============================================================
# Rango de la figura
# ============================================================

xi_max_plot = 1.6
omega_min_plot = -5
omega_max_plot = 5

mask_xi = xi_vals <= xi_max_plot
mask_omega = (omega >= omega_min_plot) & (omega <= omega_max_plot)

xi_plot = xi_vals[mask_xi]
omega_plot = omega[mask_omega]
I_plot = I_spec[mask_xi, :][:, mask_omega]

num_curvas = 14
indices = np.linspace(0, len(xi_plot)-1, num_curvas).astype(int)

# ============================================================
# Gráfico 3D tipo waterfall del espectro de intensidad
# ============================================================

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111, projection="3d")

for idx in indices:
    x = omega_plot
    y = np.full_like(omega_plot, xi_plot[idx])
    z = I_plot[idx, :]

    ax.plot(x, y, z, color="black", linewidth=0.9)

ax.set_xlabel(r"$(\omega-\omega_0)T_0$")
ax.set_ylabel(r"$z/L_D$")
ax.set_zlabel("Intensity")

ax.set_xlim(omega_min_plot, omega_max_plot)
ax.set_ylim(0, xi_max_plot)
ax.set_zlim(0, 1.05)

ax.view_init(elev=22, azim=-60)

ax.set_title(r"Spectral intensity evolution, $N=3$")

plt.tight_layout()
plt.savefig("espectro_intensidad_waterfall.png", dpi=300)
plt.show()
