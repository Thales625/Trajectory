import matplotlib.pyplot as plt

from Vector import Vector2

GM = 3.5316 * (10**12)
rb = 600000

def a(x): # dv/dt
    dist = x.magnitude()
    return -GM*x/dist**3

step_size = 1

r = [Vector2(0, rb + 76290.6)]
v = [Vector2(2286.7, 0)]

for step in range(1000):
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

    #print(tuple(r[-1] + (step_size/6)*(k1r + 2*k2r + 2*k3r + k4r)))


# graph
plt.style.use("dark_background")
fig, ax = plt.subplots(1, 1)

ax.axis("equal")
ax.add_patch(plt.Circle((0, 0), rb, color="green")) # Body

ax.plot([_.x for _ in r], [_.y for _ in r])
#ax.scatter([_.x for _ in r], [_.y for _ in r], s=1, color="yellow")

alt_list = [_.magnitude() - rb for _ in r]

print(f"Apoapsis: {max(alt_list):0f}")
print(f"Periapsis: {min(alt_list):0f}")

ax.set_xlabel('X')
ax.set_ylabel('Y')

plt.show()