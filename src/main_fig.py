import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

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
    ax.plot([q_x[i], p_x[i]], [q_y[i], p_y[i]], [q_z[i], p_z[i]], c='r', linewidth=2, alpha=.8, zorder=5)
  # punti
  ax.plot(p_x, p_y, p_z, c='g', linewidth=5, marker="o", markersize=7, alpha=.8, label='P chain', zorder=10)
  ax.plot(q_x, q_y, q_z, c='b', linewidth=5, marker="o", markersize=7, alpha=.8, label='Q chain', zorder=0)
  ax.plot(c_x, c_y, c_z, c='r', linewidth=5, marker="o", markersize=7, alpha=.8, label='Backbones', zorder=5)
  ax.set_xlabel('X (nm)', fontdict=axs_style)
  ax.set_ylabel('Y (nm)', fontdict=axs_style)
  ax.set_zlabel('Z (nm)', fontdict=axs_style)
  # ax.legend()
  ax.axis('equal')
  # ax.set_xlim(-1.5, 1.5)
  # ax.set_ylim(-1.5, 1.5)
  # ax.set_zlim(1.8, 5.2)
  # make the panes transparent
  ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  # reduce white space
  fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
  plt.axis('off')
  ax.view_init(elev=6, azim=39, roll=0)
  # ax.view_init(elev=14, azim=-50, roll=0)
  plt.savefig(f"../fig/main_5.png", bbox_inches='tight')
  plt.savefig(f"../fig/main_5.pdf", bbox_inches='tight')
  plt.show()
  plt.close()

FILE = "./5"
# Read coordinates from three different files
c_file = FILE + "0_0_c.txt"
p_file = FILE + "0_0_p.txt"
q_file = FILE + "0_0_q.txt"

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
