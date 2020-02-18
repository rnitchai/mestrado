import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

x = []
while 'True':
    try:
        x.append(input().strip().split())
    except EOFError:
        break


x=np.array(x, dtype=np.float32)


# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

o=[0,0]

plt.ylim(-2, 2)
plt.xlim(-2, 2)
# for i in range(len(x)):
#     if x[i,4] == 1.0:
#         plt.scatter(x[i,2], x[i,3], s=0.2, c='blue')
# plt.scatter(x[0,0], x[0,1], s=5, c='red')
plt.scatter(x[:,0], x[:,1], s=0.01, marker='.', c='blue')
plt.ylabel('Ptheta')
plt.xlabel('theta')


plt.autoscale()
plt.show()