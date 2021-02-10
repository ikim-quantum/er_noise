from D1_v2 import width_pcc_dmera_1d, depth_pcc_dmera_1d
import numpy as np


n=8
D=10
widths = []
depths = []
vols = []

for x in range(2**n):
    print("x={}/{}".format(x+1, 2**n))
    supp = [x, (x+1)%(2**n)]
    width = width_pcc_dmera_1d(n,D, supp)
    depth = depth_pcc_dmera_1d(n,D, supp)
    widths.append(width)
    depths.append(depth)
    vols.append(width*depth)

print("max width={}".format(max(widths)))
print("average width={}".format(np.mean(widths)))
print("Standard deviation={}".format(np.std(widths)))

print("max depth={}".format(max(depths)))
print("average depth={}".format(np.mean(depths)))
print("Standard deviation={}".format(np.std(depths)))

print("max vol={}".format(max(vols)))
print("average vol={}".format(np.mean(vols)))
print("Standard deviation={}".format(np.std(vols)))
