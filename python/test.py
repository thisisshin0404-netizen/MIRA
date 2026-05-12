import MIRA

positions = [[1.0, 1.0], [3.0, 3.0]]
velocities = [[0.0, 0.0], [0.0, 0.0]]

dt = 0.001
box_length = 10.0

result = MIRA.step(
    positions=positions,
    velocities=velocities,
    dt=dt,
    box_length=box_length,
    epsilon=1.0,
    sigma=1.0,
)

print("positions =", result["positions"])
print("velocities =", result["velocities"])
print("PE =", result["potential_energy"])
print("KE =", result["kinetic_energy"])