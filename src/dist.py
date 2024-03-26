import matplotlib.pyplot as plt
import numpy as np

plt.rc('xtick', labelsize=10)

axs_style = {
    'color': 'black',
    'weight': 'normal',
    'size': 12,
}

PATH = "../data/brand_D_"
VAR = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
TYPE = "wrapping.txt"
# data = np.zeros(shape=(0, 1))
# for v in VAR :
#   FILE = PATH[:20] + v + "/"
#   # Load data
#   d = np.loadtxt(FILE + TYPE, delimiter='\t', skiprows=1)
#   np.nan_to_num(d, copy=False)
#   d = d[300000:-1, :]
#   data = np.append(data, d[:, 1])
# # all_data = len(data)
# # print("All data: ", all_data)
# # pos_data = len(data[data > 0])
# # print("Positive data: ", pos_data)
# # neg_data = len(data[data < 0])
# # print("Negative data: ", neg_data)
# # all_data = pos_data + neg_data
# hist, bins = np.histogram(data, bins=80, density=True)
# # for i in range(len(hist)) :
# #   hist[i] = hist[i] / len(data)
# fig, ax = plt.subplots(layout='constrained')
# ax.plot(bins[:-1], hist, "r-")
# ax.set_xlabel('Chirality C', fontdict=axs_style)
# ax.set_ylabel('Frequency', fontdict=axs_style)
# ax.set_xlim(-1.05, .01)
# # ax.set_ylim(0., 6.)
# # ax.text(-.4, 5, "{:.2f}".format(neg_data/all_data*100)+"%", fontsize=12, ha='center', va='center')
# # ax.text(.4, 5, "{:.2f}".format(pos_data/all_data*100)+"%", fontsize=12, ha='center', va='center')
# plt.savefig("../fig/hist_ch_pr.png")
# plt.savefig("../fig/hist_ch_pr.pdf")
# plt.show()

def perc(d) :
  data = np.zeros(shape=(0, 1))
  for v in VAR :
    FILE = PATH[:16] + v + "/" + d + "_40_"
    # FILE = PATH[:16] + v + "/70_" + f + "_"
    # Load data
    dd = np.loadtxt(FILE + TYPE, delimiter='\t', skiprows=1)
    np.nan_to_num(dd, copy=False)
    dd = dd[200000:-1, :]
    data = np.append(data, dd[:, 1])
  return data

pos = []
pos_err = []
neg = []
neg_err = []
D = ["5", "75", "100", "125", "150"]
# F = ["0", "20", "40", "60", "80", "100"]
for d in D :
  data = perc(d)
  pos_data = data[data > 0]
  neg_data = data[data < 0]
  all_data = len(pos_data) + len(neg_data)
  pos.append(len(pos_data)/all_data*100)
  pos_err.append(np.std(pos_data)*100/pos[-1]**.5)
  neg.append(len(neg_data)/all_data*100)
  neg_err.append(np.std(neg_data)*100/neg[-1]**.5)

# print(pos_err)
# print(neg_err)
fig, ax = plt.subplots(layout='constrained')
D_label = [0.5, 0.75, 1.0, 1.25, 1.5]
ax.errorbar(D_label, neg, neg_err, marker='o', color='b', linestyle='-', label="Left-handed")
ax.errorbar(D_label, pos, pos_err, marker='*', color='r', linestyle='--', label="Right-handed")
ax.set_xlabel(r'$D_{DNA}$ (pNnm)', fontdict=axs_style)
# ax.set_xlabel(r'$F$ (pN)', fontdict=axs_style)
ax.set_ylabel('Percentage %', fontdict=axs_style)
plt.xticks(D_label)
ax.legend()
plt.show()
