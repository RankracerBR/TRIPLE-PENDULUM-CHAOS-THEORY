#Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate  import solve_ivp

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

def kinetic_energy(vx1, vy1, vx2, vy2, vx3, vy3, m1, m2, m3):
    return 0.5 * m1 * (vx1**2 + vy1**2) + 0.5 * m2 * (vx2**2 + vy2**2) + 0.5 * m3 * (vx3**2 + vy3**2)

def gravitational_potential_energy(y1, y2, y3, m1, m2, m3, g):
    return m1 * g * y1 + m2 * g * y2 + m3 * g * y3

def mechanical_energy(vx1, vy1, vx2, vy2, vx3, vy3, y1, y2, y3, m1, m2, m3, g):
    return kinetic_energy(vx1, vy1, vx2, vy2, vx3, vy3, m1, m2, m3) + gravitational_potential_energy(y1, y2, y3, m1, m2, m3, g) + 0


#
m1 = 1.0 # Massa da partícula 1 (kg)
m2 = 1.0 # Massa da partícula 2 (kg)
m3 = 1.0 # Massa da partícula 3 (kg)
l1 = 1.0 # Comprimento da haste 1 (m)
l2 = 1.0 # Comprimento da haste 2 (m)
l3 = 1.0 # Comprimento da haste 3 (m)
g = 10# Aceleração da gravidade (m/s^2)


θ1_0 = 30 * np.pi/180
θ1p_0 = 0.0
θ2_0 = 30 * np.pi/180
θ2p_0 = 0.0
θ3_0 = 30 * np.pi/180
θ3p_0 = 0.0

y0 = [θ1_0, θ1p_0, θ2_0, θ2p_0, θ3_0, θ3p_0]

t_span = (0.0, 10.0)
t = np.linspace(0, 40, 1000)
sol = solve_ivp(deriv, t_span, y0, t_eval=np.linspace(t_span[0], t_span[1], 500))

kinetic_energies = []
potential_energies = []
mechanical_energies = []

for i in range(len(sol.t)):
    y = sol.y[:, i]
    θ1, θ1p, θ2, θ2p, θ3, θ3p = y
    
    # Calcular velocidades
    vx1 = l1*θ1p*np.cos(θ1)
    vy1 = l1*θ1p*np.sin(θ1)
    vx2 = vx1 + l2*θ2p*np.cos(θ2)
    vy2 = vy1 + l2*θ2p*np.sin(θ2)
    vx3 = vx2 + l3*θ3p*np.cos(θ3)
    vy3 = vy2 + l3*θ3p*np.sin(θ3)
    
    # Calcular energias
    kinetic = kinetic_energy(vx1, vy1, vx2, vy2, vx3, vy3, m1, m2, m3)
    potential = gravitational_potential_energy(-l1*np.cos(θ1), -l2*np.cos(θ2), -l3*np.cos(θ3), m1, m2, m3, g)
    mechanical = mechanical_energy(vx1, vy1, vx2, vy2, vx3, vy3, -l1*np.cos(θ1), -l2*np.cos(θ2), -l3*np.cos(θ3), m1, m2, m3, g)
    
    kinetic_energies.append(kinetic)
    potential_energies.append(potential)
    mechanical_energies.append([mechanical])
#
plt.figure(figsize=(20, 12))
plt.plot(sol.t, sol.y[0], label='θ1')
plt.plot(sol.t, sol.y[2], label='θ2')
plt.plot(sol.t, sol.y[4], label='θ3')
plt.xlabel('Tempo (s)')
plt.ylabel('Ângulo (rad)')
plt.legend()
plt.grid()
plt.show()

#
fig, ax = plt.subplots()

ax.plot(sol.y[0], sol.y[1], label='Pêndulo 1')
ax.plot(sol.y[2], sol.y[3], label='Pêndulo 2')
ax.plot(sol.y[4], sol.y[5], label='Pêndulo 3')

ax.set_xlabel('Ângulo θ')
ax.set_ylabel('Velocidade angular dθ/dt')
ax.set_title('Trajetórias dos pêndulos')
ax.legend()
plt.grid()
plt.show()

# Plotar o gráfico
plt.figure(figsize=(12, 8))
plt.plot(sol.t, kinetic_energies, label='Energia cinética')
plt.plot(sol.t, potential_energies, label='Energia potencial')
#plt.plot(sol.t, mechanical_energies, label='Energia mecânica') #não conservativa pois está sendo influenciada pela gravidade
plt.xlabel('Tempo (s)')
plt.ylabel('Energia (J)')
plt.legend()
plt.grid()
plt.show()
