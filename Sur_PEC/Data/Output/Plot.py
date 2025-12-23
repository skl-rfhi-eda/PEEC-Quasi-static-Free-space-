import matplotlib.pyplot as plt
import numpy as np

f=[]
S11 = []
S12 = []

with open("map.txt", 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        arr = np.array(line.split(','), dtype='float64')
        f.append(arr[0])
        S11.append(arr[1])
        S12.append(arr[3])

plt.plot(f, S11, label='S11')
plt.plot(f, S12, label='S12')

plt.xlabel("f/Hz")
plt.ylabel("S/dB")
plt.legend()
plt.savefig("1.png",dpi=1000)
plt.show()
