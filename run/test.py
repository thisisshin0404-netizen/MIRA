import MIRA

# Simple input section
n_particles = 16
box_length = 10.0
position_rule = "square"  # "square" or "random"
seed = 42

dt = 0.001
epsilon = 1.0
sigma = 1.0

system = MIRA.initialize_system(
    n_particles=n_particles,
    box_length=box_length,
    position_rule=position_rule,
    seed=seed,
)

positions = system["positions"]
velocities = system["velocities"]

print("Initial system")
print("n_particles =", system["n_particles"])
print("box_length =", system["box_length"])
print("position_rule =", system["position_rule"])
print("positions =", positions)
print("velocities =", velocities)

result = MIRA.step(
    positions=positions,
    velocities=velocities,
    dt=dt,
    box_length=box_length,
    epsilon=epsilon,
    sigma=sigma,
)

print("\nAfter one step")
print("positions =", result["positions"])
print("velocities =", result["velocities"])
print("PE =", result["potential_energy"])
print("KE =", result["kinetic_energy"])
print("Total energy =", result["total_energy"])
