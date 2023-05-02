import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

M=1

def N(r):
    return 1-2*M/r


n_data = 2000000
data = np.zeros([n_data,2])
der_data = np.zeros([n_data,2])

X = np.asarray(np.linspace(2*M,1e3*M,n_data))

data[:,0] = X
data[:,1]= N(data[:,0])

inter_data = interp1d(data[:,0], data[:,1])

der_data[:,0] = X
dr = 1e-8
for i in range(1,len(X)-1):
    der_data[i,1] = (inter_data(der_data[i,0]+dr) - inter_data(der_data[i,0]-dr))/(2*dr)

der_data[0,1] = (inter_data(der_data[0,0]+dr) - inter_data(der_data[0,0]))/(dr)
der_data[-1,1] = (inter_data(der_data[-1,0]) - inter_data(der_data[-1,0]-dr))/(dr)



np.savetxt('N.txt', data)
np.savetxt('derN.txt', der_data)
plt.figure()
plt.plot(data[:,0], data[:,1])
plt.plot(der_data[:,0], der_data[:,1])
plt.show()