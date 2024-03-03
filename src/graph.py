import numpy as np
import matplotlib.pyplot as plt

# Read coordinates from three different files
c_file = "./file_c.txt"
p_file = "./file_p.txt"
q_file = "./file_q.txt"

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
ax = fig.add_subplot(111, projection='3d')
ax.plot(c_x, c_y, c_z, c='r', linewidth=2, marker="o", markersize=5, label='Backbones')
ax.plot(p_x, p_y, p_z, c='g', linewidth=2, marker="o", markersize=5, label='P chain')
ax.plot(q_x, q_y, q_z, c='b', linewidth=2, marker="o", markersize=5, label='Q chain')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()
