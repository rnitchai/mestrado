import numpy as np
import time 
import scipy.integrate as ode 
import matplotlib.pyplot as plt

#Parâmetros
mu = 10
e=0.9
beta=mu*e/10
c=15.8
Lz=pow(c*beta*mu, 0.25)

def f0(x, mu, beta, Lz):
    return x[1]
def f1(x, mu, beta, Lz):
    phi1=(pow(x[3],2))/(pow(x[0],3)) + (pow(Lz,2))/((pow(x[0],3))*pow((np.sin(x[2])),2)) - (mu)/(pow(x[0],2)) - (3*beta*(1-3*pow((np.cos(x[2])),2)))/(pow(x[0],4))
    return phi1
def f2(x, mu, beta, Lz):
    return x[3]/(pow(x[0],2))
def f3(x, mu, beta, Lz):
    phi3=(pow(Lz,2)*np.cos(x[2]))/((pow(x[0],2))*pow((np.sin(x[2])),3)) + (6*beta*np.cos(x[2])*np.sin(x[2]))/(pow(x[0],3))
    return phi3
    

def f_x(t, x):
    f = np.empty((4,))
    f[:] = np.nan
    f[0] = x[1]
    f[1] = (pow(x[3],2))/(pow(x[0],3)) + (pow(Lz,2))/((pow(x[0],3))*pow((np.sin(x[2])),2)) - (mu)/(pow(x[0],2)) - (3*beta*(1-3*pow((np.cos(x[2])),2)))/(pow(x[0],4))
    f[2] = x[3]/(pow(x[0],2))
    f[3] = (pow(Lz,2)*np.cos(x[2]))/((pow(x[0],2))*pow((np.sin(x[2])),3)) + (6*beta*np.cos(x[2])*np.sin(x[2]))/(pow(x[0],3))
    return f

#Retorna o ponto de equilíbrio estável
def Rmais(mu, beta, Lz): 
    return (pow(Lz,2) + np.sqrt(pow(Lz, 4) - 12*beta*mu))/(2*mu)

#Retorna o ponto de equilíbrio instável
def Rmenos(mu, beta, Lz):
    return (pow(Lz,2) - np.sqrt(pow(Lz, 4) - 12*beta*mu))/(2*mu)

#Retorna energia dado um ponto e Lz
def E(x, mu, beta, Lz):
    auxE = (pow(x[1],2) + pow((x[3]/x[0]),2))/2 + (pow(Lz,2))/(2*pow(x[0]*np.sin(x[2]), 2)) - (mu)/(x[0]) - (beta*(1-3*(pow(np.cos(x[2]),2))))/(pow(x[0], 3))
    return auxE

#Dado r, theta, ptheta, retorna pr no nível de energia E
def fpr(E, x, mu, beta, Lz):
    auxxx = 2*(E - (pow(Lz,2))/(2*pow(x[0]*np.sin(x[2]),2)) + (mu)/(x[0]) + (beta*(1-3*(pow(np.cos(x[2]),2))))/(pow(x[0], 3))) - (pow(x[3], 2))/(pow(x[0],2))
    return np.sqrt(auxxx)



#Passo, Rk
hs=0.01

#Equilíbrios
r=Rmenos(mu, beta, Lz)
R=Rmais(mu, beta, Lz)
eq_inst=[r, 0, np.pi/2, 0]
eq_est=[R, 0, np.pi/2, 0]

#Definindo energia crítica (E(eq_instável))
Ec=E(eq_inst, mu, beta, Lz)
print("Energia crítica eq. instável:", Ec)
#Definindo energia crítica (E(eq_estável))
Ec_est=E(eq_est, mu, beta, Lz)
print("Energia crítica eq. estável:", Ec_est)
#Definindo energia
Energy = Ec-0.005
print("Energia:", Energy)

#Definindo Condições iniciais que estejam na sec. de Poincare
#Primeira - Vermelha
theta0= np.pi*180/360
ptheta0 = 0.001
x0=np.array([R,0.0,theta0,ptheta0])
#Segunda - Azul
theta1= np.pi*180/360
ptheta1 = 0.002
x1=np.array([R,0.0,theta1,ptheta1])
#Terceira - Verde
theta2= np.pi*180/360
ptheta2 = 0.003
x2=np.array([R,0.0,theta2,ptheta2])
#Quarta - Amarelo
theta3= np.pi*180/360
ptheta3 = 0.004
x3=np.array([R,0.0,theta3,ptheta3])

#pr que está na seção
pr0= fpr(Energy, x0, mu, beta, Lz)
pr1= fpr(Energy, x1, mu, beta, Lz)
pr2= fpr(Energy, x2, mu, beta, Lz)
pr3= fpr(Energy, x3, mu, beta, Lz)

#condição inicial e auxiliar para comparaçao com passo anterior
x0=np.array([R,pr0,theta0,ptheta0])
x1=np.array([R,pr1,theta1,ptheta1])
x2=np.array([R,pr2,theta2,ptheta2])
x3=np.array([R,pr3,theta3,ptheta3])

#Tempo final para integração
t_final = 16000
#Intervalo de tempo que ficará guardado
t_intervalo = 0.01

#Inicia clock RK
begin1 = time.perf_counter()

#Solução das EDOs através do RK45
solution0 = ode.solve_ivp(f_x, [0, t_final], x0, method = "RK45", t_eval = np.arange(0, t_final, t_intervalo), rtol = 0.000000001, atol = 0.0000000000001)
solution1 = ode.solve_ivp(f_x, [0, t_final], x1, method = "RK45", t_eval = np.arange(0, t_final, t_intervalo), rtol = 0.000000001, atol = 0.0000000000001)
solution2 = ode.solve_ivp(f_x, [0, t_final], x2, method = "RK45", t_eval = np.arange(0, t_final, t_intervalo), rtol = 0.000000001, atol = 0.0000000000001)
solution3 = ode.solve_ivp(f_x, [0, t_final], x3, method = "RK45", t_eval = np.arange(0, t_final, t_intervalo), rtol = 0.000000001, atol = 0.0000000000001)

##########   IMPRIME ÓRBITA EM (x, z)   ############
# plt.plot(solution.y[0, :]*np.sin(solution.y[2, :]), solution.y[0, :]*np.cos(solution.y[2, :]), c='red', linewidth=0.1)
# plt.autoscale()
# plt.show()

# ##########   IMPRIME ÓRBITA EM (R, pR)   ############
# plt.scatter(solution.y[0, :], solution.y[1, :], c='red', s=0.1)
# plt.autoscale()
# plt.show()

#imprime tempo operacional do RK
print('Tempo Rk:', time.perf_counter()-begin1)

#Definindo matriz que irá receber os pontos que cruzam a seção
poincare0 = np.empty((4,))
poincare1 = np.empty((4,))
poincare2 = np.empty((4,))
poincare3 = np.empty((4,))
# poincare = np.copy(x)
poincarelist0 = [list(poincare0)]
poincarelist1 = [list(poincare1)]
poincarelist2 = [list(poincare2)]
poincarelist3 = [list(poincare3)]


#Imprime tempo (da órbita) no ponto final da integração
if(t_intervalo*len(solution0.y[0,:]) < t_final):
    print("Integração da primeira CI interrompida prematuramente em t(final)=", t_intervalo*len(solution0.y[0,:]))
if(t_intervalo*len(solution0.y[0,:]) == t_final):
    print("Integração da primeira CI atingiu t(final)=", t_intervalo*len(solution0.y[0,:]))
if(t_intervalo*len(solution1.y[0,:]) < t_final):
    print("Integração da segunda CI interrompida prematuramente em t(final)=", t_intervalo*len(solution1.y[0,:]))
if(t_intervalo*len(solution1.y[0,:]) == t_final):
    print("Integração da segunda CI atingiu t(final)=", t_intervalo*len(solution1.y[0,:]))
if(t_intervalo*len(solution2.y[0,:]) < t_final):
    print("Integração da terceira CI interrompida prematuramente em t(final)=", t_intervalo*len(solution2.y[0,:]))
if(t_intervalo*len(solution2.y[0,:]) == t_final):
    print("Integração da terceira CI atingiu t(final)=", t_intervalo*len(solution2.y[0,:]))
if(t_intervalo*len(solution3.y[0,:]) < t_final):
    print("Integração da quarta CI interrompida prematuramente em t(final)=", t_intervalo*len(solution3.y[0,:]))
if(t_intervalo*len(solution3.y[0,:]) == t_final):
    print("Integração da quarta CI atingiu t(final)=", t_intervalo*len(solution3.y[0,:]))

#Inicia clock Plot
begin2 = time.perf_counter()

#Seleciona as linhas de "solution" que cruzam a seção e copia para "poincarelist"
i=1
for i in range(len(solution0.y[0,:])):
        #Verifica se R>0
        if(solution0.y[0, i]>0): 
            if(solution0.y[0, i-1] < R and solution0.y[0, i] >= R and solution0.y[1, i]>0): 
                poincarelist0.append(list(solution0.y[:,i]))

i=1
for i in range(len(solution1.y[0,:])):
        #Verifica se R>0
        if(solution1.y[0, i]>0): 
            if(solution1.y[0, i-1] < R and solution1.y[0, i] >= R and solution1.y[1, i]>0): 
                poincarelist1.append(list(solution1.y[:,i]))

i=1
for i in range(len(solution2.y[0,:])):
        #Verifica se R>0
        if(solution2.y[0, i]>0): 
            if(solution2.y[0, i-1] < R and solution2.y[0, i] >= R and solution2.y[1, i]>0): 
                poincarelist2.append(list(solution2.y[:,i]))

i=1
for i in range(len(solution3.y[0,:])):
        #Verifica se R>0
        if(solution3.y[0, i]>0): 
            if(solution3.y[0, i-1] < R and solution3.y[0, i] >= R and solution3.y[1, i]>0): 
                poincarelist3.append(list(solution3.y[:,i]))


arr0 = np.array(poincarelist0)
arr1 = np.array(poincarelist1)
arr2 = np.array(poincarelist2)
arr3 = np.array(poincarelist3)

#Plotando pontos da seção
##### PARA ALPHAS DECRESCENTES #####
# alphas = np.linspace(0.1, 1, len(arr[:,0]))
# rgba_colors = np.zeros((len(arr[:,0]),4))
# # for red the first column needs to be one
# rgba_colors[:,0] = 1.0
# # the fourth column needs to be your alphas
# rgba_colors[:, 3] = 1-alphas
# plt.scatter(arr[:,2], arr[:,3], s=0.1, c=rgba_colors)

##### PARA ALPHA Fixo #####
plt.scatter(arr0[:,2], arr0[:,3], s=0.1, c='red')
plt.scatter(arr1[:,2], arr1[:,3], s=0.1, c='blue')
plt.scatter(arr2[:,2], arr2[:,3], s=0.1, c='green')
plt.scatter(arr3[:,2], arr3[:,3], s=0.1, c='yellow')



print("Tempo Poincaré:", time.perf_counter()-begin2)
print("Número de pontos da primeira CI que cruzam a seção :", len(arr0[:,0]))
print("Número de pontos da segunda CI que cruzam a seção:", len(arr1[:,0]))
print("Número de pontos da terceira CI que cruzam a seção:", len(arr2[:,0]))
print("Número de pontos da quarta CI que cruzam a seção:", len(arr3[:,0]))
plt.legend(['Energy = %E'  %Energy])
plt.autoscale()
plt.show()              
