import matplotlib.pyplot as plt
from math import ceil, pi
from json import loads
from numpy import interp

from Vector import Vector2

# get data
with open("./Data/altitude.txt", "r") as f:
    altitude_data = loads(f.read())

with open("./Data/pressure.txt", "r") as f:
    pressure_data = loads(f.read())

with open("./Data/density.txt", "r") as f:
    density_data = loads(f.read())

with open("./Data/temperature.txt", "r") as f:
    temperature_data = loads(f.read())

def a(x): # dv/dt
    dist = x.magnitude()
    alt = dist-rb

    accel = -(GM/dist**3)*x # gravity

    # drag force
    if alt < atm_depth:
        pho = float(interp(alt, altitude_data, density_data))
        v_mag = v[-1].magnitude()
        f_drag = c*pho*v_mag**2
        accel -= (f_drag/m) * (v[-1]/v_mag)
        # accel -= (k * v_mag**2 / m) * (v[-1]/v_mag)
        # print(tuple((k * v_mag**2 / m) * (v[-1]/v_mag)))

    return accel

# params
GM = 3.5316 * (10**12)
rb = 600000

cd = 0.14
radius = 2.6
m = 1000

# consts
area = pi * radius**2
c = 0.5 * cd * area

# init
r = [Vector2(0, rb + 100000)]
v = [Vector2(2200, 0)]
atm_depth = altitude_data[-1]
RUNNING = True

# rk4
step_size = 10 # todo: dynamic step size
total_time = 1000 #24 * 60**2
steps = ceil(total_time / step_size)

for step in range(steps):

    k1v = a(r[-1])
    k1r = v[-1]

    k2v = a(r[-1]+k1r*step_size*.5)
    k2r = v[-1]+k1v*step_size*.5

    k3v = a(r[-1]+k2r*step_size*.5)
    k3r = v[-1]+k2v*step_size*.5

    k4v = a(r[-1]+k3r*step_size)
    k4r = v[-1]+k3v*step_size

    v.append(v[-1] + (step_size/6)*(k1v + 2*k2v + 2*k3v + k4v))
    r.append(r[-1] + (step_size/6)*(k1r + 2*k2r + 2*k3r + k4r))

    if r[-1].magnitude() < rb: break
    #print(tuple(r[-1] + (step_size/6)*(k1r + 2*k2r + 2*k3r + k4r)))


# graph
plt.style.use("dark_background")
fig, ax = plt.subplots(1, 1)

ax.axis("equal")
ax.set_xlabel('X')
ax.set_ylabel('Y')

# draw atm
max_pho = float(interp(0, altitude_data, density_data))
h = 0
list_color_atm = []
while True:
    pho = float(interp(h, altitude_data, density_data))
    blue = min(1, .05 + float(pho / max_pho))
    list_color_atm.append((rb+h, (0, 0, blue)))
    h+=1000
    if pho <= 0.01: break

for i in reversed(list_color_atm):
    ax.add_patch(plt.Circle((0, 0), i[0], color=i[1]))

# draw body
ax.add_patch(plt.Circle((0, 0), rb, color="green"))

# draw trajectory
ax.plot([_.x for _ in r], [_.y for _ in r])

alt_list = [_.magnitude() - rb for _ in r]

print(f"Apoapsis: {max(alt_list):0f}")
print(f"Periapsis: {min(alt_list):0f}")

plt.show()