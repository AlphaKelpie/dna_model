import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize=10)

axs_style = {
    'color': 'black',
    'weight': 'normal',
    'size': 12,
}

def plot() :
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
  # linee
  for i in range(0, len(c_x)) :
    ax.plot([q_x[i], p_x[i]], [q_y[i], p_y[i]], [q_z[i], p_z[i]], c='r', linewidth=2, alpha=.2)
  # punti
  ax.plot(p_x, p_y, p_z, c='g', linewidth=3, marker="o", markersize=4, alpha=.6, label='P chain')
  ax.plot(q_x, q_y, q_z, c='g', linewidth=3, marker="o", markersize=4, alpha=.6, label='Q chain')
  # ax.plot(c_x, c_y, c_z, c='r', linewidth=3, marker="o", markersize=4, alpha=.2, label='Backbones')
  ax.set_xlabel('X (nm)', fontdict=axs_style)
  ax.set_ylabel('Y (nm)', fontdict=axs_style)
  ax.set_zlabel('Z (nm)', fontdict=axs_style)
  # ax.legend()
  ax.axis('equal')
  ax.set_xlim(-12, 12)
  ax.set_ylim(-12, 12)
  ax.set_zlim(0, 27)
  # make the panes transparent
  ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  # reduce white space
  fig.subplots_adjust(left=0, right=1, bottom=0, top=1) 
  fig.tight_layout()
  plt.savefig(f"../fig/bw_{p}_{t}.png", bbox_inches='tight', pad_inches=.3)
  plt.savefig(f"../fig/bw_{p}_{t}.pdf", bbox_inches='tight', pad_inches=.3)
  # plt.show()
  plt.close()

# File name
PATH = "./bend_writhe/"
THETA = [0, 1, 3, 5]
PHI = [-4, 4]
for p in PHI :
  for t in THETA :
    FILE = PATH + str(p) + "_" + str(t)
    # Read coordinates from three different files
    c_file = FILE + "_c.txt"
    p_file = FILE + "_p.txt"
    q_file = FILE + "_q.txt"

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

FILE = PATH + "energy.txt"
data = np.loadtxt(FILE, delimiter='\t', skiprows=1)
x = data[:, 0]
t_0 = data[:, 1]
t_1 = data[:, 2]
t_3 = data[:, 3]
t_5 = data[:, 4]
fig = plt.figure(layout='constrained')
ax = fig.add_subplot()
ax.plot(x, t_1, c='b', linewidth=1, label='Theta 1', linestyle='solid')
ax.plot(x, t_3, c='g', linewidth=1, label='Theta 3', linestyle='dashed')
ax.plot(x, t_5, c='r', linewidth=1, label='Theta 5', linestyle='dashdot')
ax.set_xlabel(r'Dihedral angle $\Phi$ (deg)', fontdict=axs_style)
ax.set_ylabel('Energy (pNnm)', fontdict=axs_style)
ax.legend()
plt.savefig("../fig/bw_energy.png", bbox_inches='tight')
plt.savefig("../fig/bw_energy.pdf", bbox_inches='tight')
# plt.show()
