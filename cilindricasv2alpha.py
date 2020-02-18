import numpy as np
import time 
import scipy.integrate as ode 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab
from matplotlib import rc
from mpl_toolkits import mplot3d

#Parametros
e=0.43
eaux=e
l=1.05
Energy = -0.15

def f_x(t, x):
    f = np.empty((4,))
    f[:] = np.nan
    f[0] = x[1]
    f[1] = -(x[0]/(pow(x[0]**2 + x[2]**2, 3/2))-(e**2)*x[0]/(5*pow(x[0]**2 + x[2]**2, 5/2))+(e**2)*x[0]*(x[0]**2 - 2*(x[2]**2))/(2*pow(x[0]**2 + x[2]**2, 7/2))- l**2/pow(x[0],3))
    f[2] = x[3]
    f[3] = -(x[2]/(pow(x[0]**2 + x[2]**2, 3/2))+2*(e**2)*x[2]/(5*pow(x[0]**2 + x[2]**2, 5/2))+(e**2)*x[2]*(x[0]**2 - 2*(x[2]**2))/(2*pow(x[0]**2 + x[2]**2, 7/2)))
    return f
def pot(x):
    auxpot = -1/((x[0]**2 + x[2]**2)**(1/2)) - (e**2)*(x[0]**2 - 2*(x[2]**2))/(10*(x[0]**2 + x[2]**2)**(5/2)) + (l**2)/(2*(x[0]**2))
    return auxpot


################## Condicoes Iniciais ########################
r0=2.5
z0=3.37749/2
x=[r0,0,z0,0]
pr0=z0*np.sqrt(2*(Energy-pot(x)))/(np.sqrt(r0**2 + z0**2))
pz0=-r0*np.sqrt(2*(Energy-pot(x)))/(np.sqrt(r0**2 + z0**2))
x=[r0,pr0,z0,pz0]

#Tempo final para integracao
t_final = 7000
#Intervalo de tempo que ficara guardado
t_intervalo = 0.1

#Inicia clock RK
begin1 = time.perf_counter()

#Solucao das EDOs atraves do RK45
solution = ode.solve_ivp(f_x, [0, t_final], x, method = "RK45", t_eval = np.arange(0, t_final, t_intervalo), rtol = 0.0000000001, atol = 0.0000000001)

#imprime tempo operacional do RK
print('Tempo Rk:', time.perf_counter()-begin1)

##########################  PLOT ALPHA  ######################
#Inicia clock alpha
begin4 = time.perf_counter()
#print alphaxt
plt.plot(solution.t[:], np.arccos(np.sqrt(1/(1+(solution.y[2, :]/solution.y[0, :])**2 + ((solution.y[0, :]**2)*((solution.y[0, :]*solution.y[3, :] - solution.y[1, :]*solution.y[2, :])/(solution.y[0, :]**2))/l)**2
))), c='blue', linewidth=0.5)
#imprime tempo operacional do RK
print('Tempo alpha:', time.perf_counter()-begin4)
#diagrama

######################### Plot alpha e=0 #####################
e=0

#Inicia clock RK2
begin6 = time.perf_counter()

#Solucao das EDOs atraves do RK45
solution = ode.solve_ivp(f_x, [0, t_final], x, method = "RK45", t_eval = np.arange(0, t_final, t_intervalo), rtol = 0.0000000001, atol = 0.0000000001)

#imprime tempo operacional do RK2
print('Tempo Rk:', time.perf_counter()-begin6)

#Inicia clock alpha2
begin5 = time.perf_counter()
#print alphaxt
plt.plot(solution.t[:], np.arccos(np.sqrt(1/(1+(solution.y[2, :]/solution.y[0, :])**2 + ((solution.y[0, :]**2)*((solution.y[0, :]*solution.y[3, :] - solution.y[1, :]*solution.y[2, :])/(solution.y[0, :]**2))/l)**2
))), c='red', linewidth=0.5)
#imprime tempo operacional do RK
print('Tempo alpha p/ e=0:', time.perf_counter()-begin5)
#diagrama
##############################################################
plt.title('E=%.2f; e=%.2f; lz=%.2f;  t_final=%i \n x0=(%.2f, %.2f, %.2f, %.2f);'  %(Energy, eaux, l, t_final, r0, pr0, z0, pz0))
plt.ylabel('alpha')
plt.xlabel('t')
plt.autoscale()
plt.show()

##############################################################

# #################### Plota r+ e x0############################
# rmais = (l**2 + np.sqrt(l**4 - 6*e**2/5))/2
# plt.scatter(rmais, 0, c='red', s=10)
# plt.scatter(r0,z0, c='blue', s=10)
# ##############################################################

# #Inicia clock Plot
# begin2 = time.perf_counter()
# # #########   IMPRIME ÓRBITA EM (r, z)   ############
# plt.plot(solution.y[0, :], solution.y[2, :], c='blue', linewidth=0.1)
# #imprime tempo operacional do Plot
# print('Tempo Plot:', time.perf_counter()-begin2)

# #Inicia clock CN
# begin3 = time.perf_counter()

# ################### CURVA EQUIPOTENCIAL ######################
# r_vals = np.linspace(0.1, 7, 1000)
# z_vals = np.linspace(-6, 6, 1000)
# R, Z = np.meshgrid(r_vals, z_vals)
# P = -1/((R**2 + Z**2)**(1/2)) - (e**2)*(R**2 - 2*(Z**2))/(10*(R**2 + Z**2)**(5/2)) + (l**2)/(2*(R**2))
# cp = plt.contour(R, Z, P, [Energy])
# ###############################################################


# #imprime tempo operacional do CN
# print('Tempo Curva de nivel:', time.perf_counter()-begin3)


# #Legenda
# plt.title('x0=(%.2f, %.2f, %.2f, %.2f); Energy=%.2f; e=%.2f; lz=%.2f'  %(r0, pr0, z0, pz0, Energy, e, l))
# plt.ylabel('z')
# plt.xlabel('r')
# plt.autoscale()
# plt.show()
