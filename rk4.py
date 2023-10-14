import matplotlib.pyplot as plt
from math import ceil, pi, exp

from Vector import Vector2

class Body:
    def __init__(self, name, r0, v0, radius, mass, has_atmosphere=False, rho0=0, atm_h=1, color="grey") -> None:
        self.name = name
        self.r = [r0]
        self.v = [v0]
        self.radius = radius
        self.mass = mass
        self.GM = mass * G
        self.has_atmosphere = has_atmosphere

        self.rho0 = rho0
        self.atm_h = atm_h

        self.color = color

    def air_density(self, altitude):
        return self.rho0 * exp(-altitude / self.atm_h)

def body_rk4(p, body):
    accel = Vector2(0, 0)

    for other_body in bodies:
        if other_body is not body:
            delta_pos = other_body.r[-1] - p
            dist = delta_pos.magnitude()
            direction = delta_pos / dist

            accel += (other_body.GM/(dist**2)) * direction
            # accel += (body.GM/dist**3)*delta_pos

    return accel

def bodies_calc(step_size, total_time):
    steps = ceil(total_time / step_size)

    for step in range(steps):
        for body in bodies:
            k1v = body_rk4(body.r[-1], body)
            k1r = body.v[-1]

            k2v = body_rk4(body.r[-1]+k1r*step_size*.5, body)
            k2r = body.v[-1]+k1v*step_size*.5

            k3v = body_rk4(body.r[-1]+k2r*step_size*.5, body)
            k3r = body.v[-1]+k2v*step_size*.5

            k4v = body_rk4(body.r[-1]+k3r*step_size, body)
            k4r = body.v[-1]+k3v*step_size

            dv = (step_size/6)*(k1v + 2*k2v + 2*k3v + k4v)
            dr = (step_size/6)*(k1r + 2*k2r + 2*k3r + k4r)

            body.v.append(body.v[-1] + dv)
            body.r.append(body.r[-1] + dr)

def vessel_rk4(x, v):
    accel = Vector2(0, 0)

    for body in bodies:
        delta = body.r[-1] - x
        dist = delta.magnitude()

        accel += (body.GM/dist**3)*delta # gravity

        # drag force
        if body.has_atmosphere:
            rho = body.air_density(dist - body.radius)
            v_mag = v.magnitude()
            f_drag = c*rho*v_mag**2
            accel -= (f_drag/m) * (v/v_mag)
        
    return accel

def vessel_calc(r0, v0, step_size, total_time):
    r = [r0]
    v = [v0]
    for step in range(ceil(total_time/step_size)):
        k1v = vessel_rk4(r[-1], v[-1])
        k1r = v[-1]

        k2v = vessel_rk4(r[-1] + k1r*step_size*.5, v[-1])
        k2r =   v[-1] + k1v*step_size*.5

        k3v = vessel_rk4(r[-1] + k2r*step_size*.5, v[-1])
        k3r =   v[-1] + k2v*step_size*.5

        k4v = vessel_rk4(r[-1] + k3r*step_size, v[-1])
        k4r =   v[-1] + k3v*step_size

        v.append(v[-1] + (step_size/6)*(k1v + 2*k2v + 2*k3v + k4v))
        r.append(r[-1] + (step_size/6)*(k1r + 2*k2r + 2*k3r + k4r))

    return r, v

# vessel
cd = 0.14
radius = 2.6
m = 1000

# consts
area = pi * radius**2
c = 0.5 * cd * area
G = 6.674264790194568e-11

# bodies load
earth = Body("Earth", Vector2(0, 0), Vector2(0, 0), 6371e3, 5.9722e24, True, 1.225, 7400, "green")
bodies = [earth]
bodies.append(Body("Moon", Vector2(0, earth.radius+384400000), Vector2(1000, 0), 1737400, 7.35e22, False, "grey"))

# vessel position
r0 = Vector2(0, earth.radius + 400000)
v0 = Vector2(7660, 0)

# propagate orbit
r, v = vessel_calc(r0, v0, 1, 1000)
bodies_calc(10000, 10000000)

# maneuver
maneuver_index = 0
maneuver_delta_v = Vector2(1000, 0)
man_r, man_v = vessel_calc(r[maneuver_index], v[maneuver_index] + maneuver_delta_v, 50, 10000)

# graph
plt.style.use("dark_background")
fig, ax = plt.subplots(1, 1)

ax.axis("equal")
ax.set_xlabel('X')
ax.set_ylabel('Y')

# draw bodies
for body in bodies:
    if body.has_atmosphere:
        h = 0
        list_color_atm = []
        while True:
            rho = body.air_density(h)
            blue = min(1, .05 + float(rho / body.rho0))
            list_color_atm.append((body.radius+h, (0, 0, blue)))
            h+=1000
            if rho <= 1e-5: break

        for i in reversed(list_color_atm):
            ax.add_patch(plt.Circle((0, 0), i[0], color=i[1]))
        list_color_atm.clear()

    # draw body
    ax.add_patch(plt.Circle(body.r[0], body.radius, color=body.color))

    # draw body orbit
    ax.plot([_.x for _ in body.r], [_.y for _ in body.r], color="purple")

# draw vessel trajectory
ax.plot([_.x for _ in r], [_.y for _ in r], color="orange")

# draw maneuver
# ax.add_patch(plt.Circle(r[maneuver_index], 100, color="grey"))
# ax.plot([_.x for _ in man_r], [_.y for _ in man_r], color="red")

# a = v[maneuver_index].normalize()
# b = a.copy()
# alpha = a.x
# a.x = a.y
# a.y = -alpha

# ax.plot([r[maneuver_index].x, r[maneuver_index].x+a.x*1000], [r[maneuver_index].y, r[maneuver_index].y+a.y*1000], color="green")
# ax.plot([r[maneuver_index].x, r[maneuver_index].x+b.x*1000], [r[maneuver_index].y, r[maneuver_index].y+b.y*1000], color="blue")

plt.show()