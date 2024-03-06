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
  ax.plot(p_x, p_y, p_z, c='g', linewidth=3, marker="o", markersize=4, alpha=.4, label='P chain')
  ax.plot(q_x, q_y, q_z, c='g', linewidth=3, marker="o", markersize=4, alpha=.4, label='Q chain')
  # ax.plot(c_x, c_y, c_z, c='r', linewidth=3, marker="o", markersize=4, alpha=.2, label='Backbones')
  for i in range(0, len(c_x)) :
    ax.plot([q_x[i], p_x[i]], [q_y[i], p_y[i]], [q_z[i], p_z[i]], c='r', linewidth=2, alpha=.2)
  ax.set_xlabel('X (nm)', fontdict=axs_style)
  ax.set_ylabel('Y (nm)', fontdict=axs_style)
  ax.set_zlabel('Z (nm)', fontdict=axs_style)
  # ax.legend()
  ax.axis('equal')
  # ax.set_xlim(-5, 5)
  # ax.set_ylim(-5, 5)
  # make the panes transparent
  ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  # make the grid lines transparent
  # ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
  # ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
  # ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
  # ax.set_axis_off()
  # reduce white space
  fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
  # plt.savefig(FILE + ".png")
  # plt.savefig(FILE + ".pdf")
  plt.show()
  plt.close()

def plot_p() :
  fig = plt.figure(layout='constrained')
  ax = fig.add_subplot()
  ax.plot(e, v, c='r', linewidth=2)
  ax.set_xlabel('Steps', fontdict=axs_style)
  ax.set_ylabel(f'{p}', fontdict=axs_style)
  # ax.legend()
  # plt.savefig(FILE + ".png")
  # plt.savefig(FILE + ".pdf")
  plt.show()
  plt.close()

# File name
PATH = "./histone/"
STEPS = [0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000]
for s in STEPS :
  FILE = PATH + str(s)
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

PARMS = ["energy", "wrapping", "chirality"]
for p in PARMS :
  # Read coordinates from file
  file = PATH + p + ".txt"

  # Read coordinates from file1
  e_data = np.loadtxt(file, delimiter='\t', skiprows=1)
  e = e_data[:, 0]
  v = e_data[:, 1]
  plot_p()
