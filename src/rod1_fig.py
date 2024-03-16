import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rc('xtick', labelsize=10)

axs_style = {
    'color': 'black',
    'weight': 'normal',
    'size': 12,
}

def plot() :
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
  # cilindro
  ax.plot_surface(s_x, s_y, s_z, color='k', alpha=0.4, cmap=cm.coolwarm, linewidth=0, antialiased=False, zorder=0)
  # linee
  for i in range(0, len(c_x)) :
    ax.plot([q_x[i], p_x[i]], [q_y[i], p_y[i]], [q_z[i], p_z[i]], c='r', linewidth=2, alpha=.2, zorder=5)
  # punti
  ax.plot(p_x, p_y, p_z, c='g', linewidth=3, marker="o", markersize=4, alpha=.6, label='P chain', zorder=10)
  ax.plot(q_x, q_y, q_z, c='g', linewidth=3, marker="o", markersize=4, alpha=.6, label='Q chain', zorder=10)
  # ax.plot(c_x, c_y, c_z, c='r', linewidth=3, marker="o", markersize=4, alpha=.2, label='Backbones')
  ax.set_xlabel('X (nm)', fontdict=axs_style)
  ax.set_ylabel('Y (nm)', fontdict=axs_style)
  ax.set_zlabel('Z (nm)', fontdict=axs_style)
  # ax.legend()
  ax.axis('equal')
  ax.set_xlim(-7, 7)
  ax.set_ylim(-7, 7)
  ax.set_zlim(0, 23)
  # make the panes transparent
  ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  # ax.set_axis_off()
  # reduce white space
  fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
  ax.view_init(azim=0)
  plt.savefig(f"../fig/r1_{sgm[r]}_{s}.png", bbox_inches='tight')
  plt.savefig(f"../fig/r1_{sgm[r]}_{s}.pdf", bbox_inches='tight')
  # plt.savefig(f"../fig/r1_{sgm[r]}_2000000_{s}.png", bbox_inches='tight')
  # plt.savefig(f"../fig/r1_{sgm[r]}_2000000_{s}.pdf", bbox_inches='tight')
  # plt.show()
  plt.close()

def plot_p() :
  fig = plt.figure(layout='constrained')
  ax = fig.add_subplot()
  ax.plot(e, v, c='r', linewidth=2)
  ax.set_xlabel('Steps', fontdict=axs_style)
  if p == "energy" :
    ax.set_ylabel('Energy (pNnm)', fontdict=axs_style)
  elif p == "wrapping" :
    ax.set_ylabel('Wrapping', fontdict=axs_style)
  # ax.legend()
  plt.savefig(f"../fig/r1_{sgm[r]}_{p}.png", bbox_inches='tight')
  plt.savefig(f"../fig/r1_{sgm[r]}_{p}.pdf", bbox_inches='tight')
  # plt.show()
  plt.close()

# File name
PATH = "../data/rod_1_0/"
STEPS = [0, 500000, 750000, 1000000, 1250000, 1500000, 1750000, 2000000]
# STEPS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
PARMS = ["energy", "wrapping"]
RADIUSES = [1., 1.5, 2.]
# Make data
z = np.linspace(0, 34, 200)
theta = np.linspace(0, 2*np.pi, 50)
theta_grid, z_grid=np.meshgrid(theta, z)

sgm = {1. : 'A_', 1.5 : 'B_', 2. : 'C_'}

for r in RADIUSES :
  for s in STEPS :
    s_x = r * np.cos(theta_grid)
    s_y = r * np.sin(theta_grid)
    s_z = z_grid
    
    FILE = PATH + sgm[r] + str(s)
    # FILE = PATH[:14] + str(s) + '/' + sgm[r] + str(2000000)
    # Read coordinates from three different files
    c_file = FILE + "_0_0_c.txt"
    p_file = FILE + "_0_0_p.txt"
    q_file = FILE + "_0_0_q.txt"

    # Read coordinates from file1
    c_data = np.loadtxt(c_file)
    c_x = c_data[:, 0]
    c_y = c_data[:, 1]
    c_z = c_data[:, 2]

    # Read coordinates from file2
    p_data = np.loadtxt(p_file)
    p_x = p_data[:, 0]
    p_y = p_data[:, 1]
    p_z = p_data[:, 2]

    # Read coordinates from file3
    q_data = np.loadtxt(q_file)
    q_x = q_data[:, 0]
    q_y = q_data[:, 1]
    q_z = q_data[:, 2]
    plot()

  for p in PARMS :
    # Read coordinates from file
    file = PATH + sgm[r] + p

    # Read coordinates from file1
    e_data = np.loadtxt(file + ".txt", delimiter='\t', skiprows=1)
    e = e_data[:, 0]
    v = e_data[:, 1]
    plot_p()
