import matplotlib.pyplot as plt
from math import floor, pi
from numpy import log as ln
from numpy import interp
from json import loads

from Vector import Vector2, Vector3 # github.com/Thales625/Vector

# get data
with open("./Data/altitude.txt", "r") as f:
    altitude_data = loads(f.read())

with open("./Data/pressure.txt", "r") as f:
    pressure_data = loads(f.read())

with open("./Data/density.txt", "r") as f:
    density_data = loads(f.read())

with open("./Data/temperature.txt", "r") as f:
    temperature_data = loads(f.read())

# define
half_pi = pi*.5

# body
g0 = 9.806
GM = 3531600035840
rb = 600000

# PARAMS
isp_sea = 240
isp_vac = 270
thrust = 89028.625
mdry = 1000
mfuel = 5000

cd = 0.139
radius = 2

# CONSTS
mrate = thrust / (isp_sea*g0)
tburning = mfuel / mrate
atm_depth = altitude_data[-1]

area = pi * radius**2
c = 0.5 * cd * area

# auto pilot
start_turn = 100
end_turn = 40000
delta_turn = end_turn - start_turn

# graph
list_x = []
list_y = []
list_v = []

# initializing
s = Vector2(0, rb)
v = Vector2(0, 0)
d = s.normalize() # direction

t = 0
dt = 0.1

s_mag = s.magnitude()
dv = (isp_sea*g0) * ln((mdry+mfuel) / mdry)
twr = (thrust / (mdry+mfuel)) / (GM / s_mag**2)

index_burning = floor(tburning/dt)
max_alt = s_mag

print("dV:", dv)
print("twr:", twr)
print("t_burning:", tburning)

while s_mag >= rb and t < 20000:
    t += dt
    m = mdry + max(0, mfuel)
    a = -(s/s_mag) * (GM / s_mag**2)
    alt = s_mag - rb

    v_mag = v.magnitude()
    v_norm = v / (v_mag if v_mag != 0 else 1)

    # direction
    if alt < start_turn: pitch = 0
    if start_turn <= alt and alt <= end_turn: pitch = (alt-start_turn)/delta_turn
    if alt > end_turn: 1
    d = Vector2(thetta=s.normalize().get_angle() - pitch * half_pi)
    
    # drag force
    if alt < atm_depth:
        # pho = 1.113 * exp(-1.24 / 10000 * alt) # atm density
        pho = interp(alt, altitude_data, density_data)
        k = c * pho
        a -= (k * v_mag**2 / m) * v_norm

    # thrust
    if mfuel > 0:
        isp = isp_vac + (isp_sea - isp_vac) * (2*pho/3)
        ve = isp * g0
        thrust_ = ve * mrate # considerando o isp variando
        a += (thrust_ / m) * d
        mfuel -= mrate*dt

    v += a * dt
    s += v * dt 

    list_x.append(s.x)
    list_y.append(s.y)
    list_v.append(v_mag)

    s_mag = s.magnitude()
    if s_mag > max_alt: max_alt = s_mag

# display info
print(f"Altitude Max: {max_alt-rb}")
print(f"Altitude MECO: {Vector2(list_x[index_burning], list_y[index_burning]).magnitude()-rb}")
print(f"Velocity MECO: {list_v[index_burning]}")

# graph
plt.style.use("dark_background")

fig, ax = plt.subplots(1, 1)


max_pho = interp(rb, altitude_data, density_data)
h = 0
while True:
    radius = rb + h
    pho = interp(radius, altitude_data, density_data)
    color = Vector3(255 * float(pho / max_pho), 0, 0)
    ax.add_patch(plt.Circle((0, 0), radius, color=tuple(color.normalize()))) # Atm 60km
    h+=10
    if pho <= 0.01: break
    

# ax.add_patch(plt.Circle((0, 0), rb+60000, color=(0, 0, .2))) # Atm 60km
# ax.add_patch(plt.Circle((0, 0), rb+20000, color=(0, 0, .25))) # Atm 20km
# ax.add_patch(plt.Circle((0, 0), rb+8000, color=(0, 0, .5))) # Atm 8km
ax.add_patch(plt.Circle((0, 0), rb+20, color="w")) # Body Border
ax.add_patch(plt.Circle((0, 0), rb, color="g")) # Body

ax.scatter(list_x[:index_burning], list_y[:index_burning], s=1, color="red")
ax.scatter(list_x[index_burning:], list_y[index_burning:], s=1, color="yellow")
ax.axis("equal")
# plt.plot(list_x, list_y, "blue")

ax.ticklabel_format(useOffset=False, style='plain', axis='both') # remove cientific notation

ax.set_xlabel('X')
ax.set_ylabel('Y')

plt.show()
