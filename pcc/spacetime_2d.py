from D2 import width_pcc_dmera_2d, depth_pcc_dmera_2d, width_pcc_dmera_2d_nocompression
import numpy as np


n=5
D=10
widths = []
widths_nc = []
depths = []
vols = []


for x in range(2**n):
    for y in range(2**n):
        print("x={}/{}, y={}/{}".format(x+1, 2**n, y+1, 2**n))
        supp_x = [(x,y), ((x+1)%(2**n),y)]
        supp_y = [(x,y), (x,(y+1)%(2**n))]

#        width = width_pcc_dmera_2d(n,D, supp_x)
        width_nc = width_pcc_dmera_2d_nocompression(n, D, supp_x)
#        depth = depth_pcc_dmera_2d(n,D, supp_x)
#        widths.append(width)
        widths_nc.append(width_nc)
#        depths.append(depth)
#        vols.append(width*depth)

#        width = width_pcc_dmera_2d(n,D, supp_y)
        width_nc = width_pcc_dmera_2d_nocompression(n, D, supp_y)
#        depth = depth_pcc_dmera_2d(n,D, supp_y)
#        widths.append(width)
        widths_nc.append(width_nc)
#        depths.append(depth)
#        vols.append(width*depth)


print("max width (uncompressed)={}".format(max(widths_nc)))
print("average width (uncompresed)={}".format(np.mean(widths_nc)))
print("Standard deviation={}".format(np.std(widths)))
        
print("max width={}".format(max(widths)))
print("average width={}".format(np.mean(widths)))
print("Standard deviation={}".format(np.std(widths)))

print("max depth={}".format(max(depths)))
print("average depth={}".format(np.mean(depths)))
print("Standard deviation={}".format(np.std(depths)))

print("max vol={}".format(max(vols)))
print("average vol={}".format(np.mean(vols)))
print("Standard deviation={}".format(np.std(vols)))
