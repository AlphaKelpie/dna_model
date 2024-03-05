import numpy as np
import matplotlib.pyplot as plt

# File name
FILE = "./histone/1000000_0_0"
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

# Plot the coordinates
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(p_x, p_y, p_z, c='g', linewidth=2, marker="o", markersize=5, alpha=.2, label='P chain')
ax.plot(q_x, q_y, q_z, c='b', linewidth=2, marker="o", markersize=5, alpha=.2, label='Q chain')
# ax.plot(c_x, c_y, c_z, c='r', linewidth=2, marker="o", markersize=5, alpha=.2, label='Backbones')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
ax.axis('equal')
plt.show()
