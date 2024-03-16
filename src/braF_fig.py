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
  # linee A
  for i in range(0, len(A_c_x)) :
    ax.plot([A_q_x[i], A_p_x[i]], [A_q_y[i], A_p_y[i]], [A_q_z[i], A_p_z[i]], c='r', linewidth=2, alpha=.2, zorder=5)
  # punti A
  ax.plot(A_p_x, A_p_y, A_p_z, c='g', linewidth=3, marker="o", markersize=4, alpha=.6, label='P chain A', zorder=10)
  ax.plot(A_q_x, A_q_y, A_q_z, c='g', linewidth=3, marker="o", markersize=4, alpha=.6, label='Q chain A', zorder=10)
  # ax.plot(c_x, c_y, c_z, c='r', linewidth=3, marker="o", markersize=4, alpha=.2, label='Backbones')
  # linee B
  for i in range(0, len(B_c_x)) :
    ax.plot([B_q_x[i], B_p_x[i]], [B_q_y[i], B_p_y[i]], [B_q_z[i], B_p_z[i]], c='r', linewidth=2, alpha=.2, zorder=5)
  # punti B
  ax.plot(B_p_x, B_p_y, B_p_z, c='b', linewidth=3, marker="o", markersize=4, alpha=.6, label='P chain B', zorder=10)
  ax.plot(B_q_x, B_q_y, B_q_z, c='b', linewidth=3, marker="o", markersize=4, alpha=.6, label='Q chain B', zorder=10)
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
  # save = A_FILE.replace("A_", "")
  # plt.savefig(save + ".png")
  # plt.savefig(save + ".pdf")
  plt.show()
  plt.close()

def plot_p() :
  fig = plt.figure(layout='constrained')
  ax = fig.add_subplot()
  ax.plot(e, v, c='r', linewidth=2)
  ax.set_xlabel('Steps', fontdict=axs_style)
  ax.set_ylabel(f'{p}', fontdict=axs_style)
  # ax.legend()
  # plt.savefig(file + ".png")
  # plt.savefig(file + ".pdf")
  plt.show()
  plt.close()

# File name
PATH = "../data/brand/"
DIST = 70
STEPS = [500000, 750000, 1000000, 1250000, 1500000, 1750000, 2000000]
PARMS = ["energy", "wrapping"]
F = [0, 20, 40, 60, 80, 100]
for f in F :
  for s in STEPS :
    A_FILE = PATH + str(DIST) + '_' + str(f) + '_A_' + str(s)
    # Read coordinates from three different files
    A_c_file = A_FILE + "_0_0_c.txt"
    A_p_file = A_FILE + "_0_0_p.txt"
    A_q_file = A_FILE + "_0_0_q.txt"

    # Read coordinates from file1
    A_c_data = np.loadtxt(A_c_file)
    A_c_x = A_c_data[:, 0]
    A_c_y = A_c_data[:, 1]
    A_c_z = A_c_data[:, 2]

    # Read coordinates from file2
    A_p_data = np.loadtxt(A_p_file)
    A_p_x = A_p_data[:, 0]
    A_p_y = A_p_data[:, 1]
    A_p_z = A_p_data[:, 2]

    # Read coordinates from file3
    A_q_data = np.loadtxt(A_q_file)
    A_q_x = A_q_data[:, 0]
    A_q_y = A_q_data[:, 1]
    A_q_z = A_q_data[:, 2]

    B_FILE = PATH + str(DIST) + '_' + str(f) + '_B_' + str(s)
    # Read coordinates from three different files
    B_c_file = B_FILE + "_0_0_c.txt"
    B_p_file = B_FILE + "_0_0_p.txt"
    B_q_file = B_FILE + "_0_0_q.txt"

    # Read coordinates from file1
    B_c_data = np.loadtxt(B_c_file)
    B_c_x = B_c_data[:, 0]
    B_c_y = B_c_data[:, 1]
    B_c_z = B_c_data[:, 2]

    # Read coordinates from file2
    B_p_data = np.loadtxt(B_p_file)
    B_p_x = B_p_data[:, 0]
    B_p_y = B_p_data[:, 1]
    B_p_z = B_p_data[:, 2]

    # Read coordinates from file3
    B_q_data = np.loadtxt(B_q_file)
    B_q_x = B_q_data[:, 0]
    B_q_y = B_q_data[:, 1]
    B_q_z = B_q_data[:, 2]

    plot()

  for p in PARMS :
    file = PATH + str(DIST) + '_' + str(f) + '_' + p
    # Read coordinates from file1
    e_data = np.loadtxt(file + ".txt", delimiter='\t', skiprows=1)
    e = e_data[:, 0]
    v = e_data[:, 1]
    plot_p()
