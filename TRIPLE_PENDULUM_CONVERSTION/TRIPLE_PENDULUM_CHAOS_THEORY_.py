#Imports
import numpy as np # Biblioteca para cáculos numéricos
import matplotlib.pyplot as plt # Biblioteca para criação de gráficos
from scipy.integrate  import solve_ivp # Função para integração de equações diferenciais

#Funções
def deriv(t, y):
    θ1, θ1p, θ2, θ2p, θ3, θ3p = y
    dθ1dt = θ1p
    dθ2dt = θ2p
    dθ3dt = θ3p
    dθ1pdt = (-g*(2*m1+m2)*np.sin(θ1) - m2*g*np.sin(θ1-2*θ2) - 2*np.sin(θ1-θ2)*m2*(θ2p**2*l2 + θ1p**2*l1*np.cos(θ1-θ2))) / (l1*(2*m1+m2-m2*np.cos(2*θ1-2*θ2)))
    dθ2pdt = (2*np.sin(θ1-θ2)*((θ1p**2)*l1*(m1+m2) + g*(m1+m2)*np.cos(θ1) + (θ2p**2)*l2*m2*np.cos(θ1-θ2))) / (l2*(2*m1+m2-m2*np.cos(2*θ1-2*θ2)))
    dθ3pdt = -g/l3*np.sin(θ3)
    
    return dθ1dt, dθ1pdt, dθ2dt, dθ2pdt, dθ3dt, dθ3pdt

# Definindo as propriedades
m1 = 1 # Massa do pêndulo 1 (kg)
m2 = 1 # Massa do pêndulo 2 (kg)
m3 = 1 # Massa do pêndulo 3 (kg)
l1 = 1 # Comprimento da haste 1 (m)
l2 = 10 # Comprimento da haste 2 (m)
l3 = 1 # Comprimento da haste 3 (m)
g = 10 # Aceleração da gravidade (m/s^2)

# Definindo as condições iniciais do sistema
θ1_0 = 30 * np.pi/180 # Ângulo inicial do pêndulo superior em radianos
θ1p_0 = 0.0           # Velocidade angular inicial do pêndulo superior em rad/s
θ2_0 = 30 * np.pi/180 # Ângulo inicial do meio superior em radianos
θ2p_0 = 0.0           # Velocidade angular inicial do pêndulo do meio em rad/s
θ3_0 = 30 * np.pi/180 # Ângulo inicial do pêndulo inferior em radianos
θ3p_0 = 0.0           # Velocidade angular inicial do pêndulo do inferior em rad/s

y0 = [θ1_0, θ1p_0, θ2_0, θ2p_0, θ3_0, θ3p_0] # Junta os parâmetros iniciais

t_span = (0.0, 40.0) # Define o tempo de intervalo
t = np.linspace(0, 50,3) # Cria um array de mil pontos entre 0 até 40 segundos
sol = solve_ivp(deriv, t_span, y0, t_eval=np.linspace(t_span[0], t_span[1], 500)) # Faz o cálculo diferencial a partir da função 'deriv'


#Gráfico 1
plt.figure(figsize=(20, 12))
plt.plot(sol.t, sol.y[0], label='θ1')
plt.plot(sol.t, sol.y[2], label='θ2')
plt.plot(sol.t, sol.y[4], label='θ3')
plt.xlabel('Tempo (s)')
plt.ylabel('Ângulo (rad)')
plt.legend()
plt.grid()
plt.show()
