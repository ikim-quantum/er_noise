from D2 import width_pcc_dmera_2d, depth_pcc_dmera_2d, width_pcc_dmera_2d_nocompression
import numpy as np


n=5
D=4
ls= [l+2 for l in range(9)]
widths_nc = []

for l in ls:
    for x in range(2**n):
        for y in range(2**n):
            #print("x={}/{}, y={}/{}".format(x+1, 2**n, y+1, 2**n))
            supp = [((xp+x)%(2**n),(yp+y)%(2**n)) for xp in range(l) for yp in range(l)]
            width_nc = width_pcc_dmera_2d_nocompression(n, D, supp)
            widths_nc.append(width_nc)



    print("l={}".format(l))
    print("|Supp|={}".format(len(supp)))
    print("max width (uncompressed)={}".format(max(widths_nc)))
    print("average width (uncompresed)={}".format(np.mean(widths_nc)))
    print("Standard deviation={}".format(np.std(widths_nc)))
