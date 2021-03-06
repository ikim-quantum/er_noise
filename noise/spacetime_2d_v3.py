from D2_v3 import width_pcc_dmera_2d, depth_pcc_dmera_2d
import numpy as np


n=3
D=4
widths = []
depths = []
vols = []


for x in range(3**n):
    for y in range(3**n):
        print("x={}/{}, y={}/{}".format(x+1, 3**n, y+1, 3**n))
        supp_x = [(x,y), ((x+1)%(3**n),y)]
        supp_y = [(x,y), (x,(y+1)%(3**n))]

        width = width_pcc_dmera_2d(n,D, supp_x)
        depth = depth_pcc_dmera_2d(n,D, supp_x)
        widths.append(width)
        depths.append(depth)
        vols.append(width*depth)

        width = width_pcc_dmera_2d(n,D, supp_y)
        depth = depth_pcc_dmera_2d(n,D, supp_y)
        widths.append(width)
        depths.append(depth)
        vols.append(width*depth)

for w in widths:
    print(w)
        
print("max width={}".format(max(widths)))
print("average width={}".format(np.mean(widths)))
print("Standard deviation={}".format(np.std(widths)))

print("max depth={}".format(max(depths)))
print("average depth={}".format(np.mean(depths)))
print("Standard deviation={}".format(np.std(depths)))

print("max vol={}".format(max(vols)))
print("average vol={}".format(np.mean(vols)))
print("Standard deviation={}".format(np.std(vols)))
