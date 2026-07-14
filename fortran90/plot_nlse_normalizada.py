import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Cargar archivo generado por Fortran
data = np.loadtxt("nlse_normalizada.dat")

# Columnas:
# 0 = xi = z/LD
# 1 = tau = T/T0
# 2 = Re(u)
# 3 = Im(u)
# 4 = |u|^2
xi_all = data[:, 0]
tau_all = data[:, 1]
intensity_all = data[:, 4]

# Reconstruir grillas
tau_vals = np.unique(tau_all)
nt = len(tau_vals)

nxi = len(intensity_all) // nt
xi_vals = xi_all[::nt]

intensity = intensity_all.reshape((nxi, nt))

# Opcional: recortar la ventana temporal para verla mejor
tau_min_plot = -5.0
tau_max_plot = 5.0

mask_tau = (tau_vals >= tau_min_plot) & (tau_vals <= tau_max_plot)

tau_plot = tau_vals[mask_tau]
intensity_plot = intensity[:, mask_tau]

# Crear gráfico 3D tipo waterfall
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")

# Dibujar una curva cada cierto número de pasos en xi
step = max(1, nxi // 35)

for j in range(0, nxi, step):
    x = tau_plot
    y = np.full_like(tau_plot, xi_vals[j])
    z = intensity_plot[j, :]

    ax.plot(x, y, z, linewidth=1.0)

# Etiquetas de los ejes
ax.set_xlabel(r"$T/T_0$")
ax.set_ylabel(r"$z/L_D$")
ax.set_zlabel("Intensity")

ax.set_title("Evolución de la intensidad - NLSE normalizada")

# Ajustar vista para parecerse más a la figura del libro
ax.view_init(elev=22, azim=-60)

# Límites sugeridos
ax.set_xlim(tau_min_plot, tau_max_plot)
ax.set_ylim(xi_vals[0], xi_vals[-1])

plt.tight_layout()
plt.savefig("waterfall_nlse.png", dpi=300)
plt.show()
