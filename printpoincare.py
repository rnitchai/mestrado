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
print(x)

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)



plt.ylim(-0.00061, 0.00061)
plt.xlim(1.57066, 1.57095)

for i in range(len(x)):
	if x[i,4] == 1.0:
		plt.scatter(x[i,2], x[i,3], s=0.2, c='blue')
plt.scatter(x[0,2], x[0,3], s=5, c='red')
#plt.scatter(x[:,2], x[:,3], s=0.2, c='blue')
plt.ylabel('Ptheta')
plt.xlabel('theta')
plt.legend()

# plt.autoscale()
plt.show()