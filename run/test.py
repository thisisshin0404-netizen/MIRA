import MIRA

# User input section
n_particles = 16
box_length = 10.0
seed = 42

n_steps = 1000
integrator = "verlet"  # "verlet" or "euler"

dt = 0.001
epsilon = 1.0
sigma = 1.0

print_every = 100

system = MIRA.initialize_system(
    n_particles=n_particles,
    box_length=box_length,
    seed=seed,
)

positions = system["positions"]
velocities = system["velocities"]

print("Initial system")
print("n_particles =", system["n_particles"])
print("box_length =", system["box_length"])
print("position_rule =", system["position_rule"])
print("seed =", system["seed"])
print("n_steps =", n_steps)
print("integrator =", integrator)
print("initial positions =", positions)
print("initial velocities =", velocities)

last_result = None
energy_history = []

for step_index in range(1, n_steps + 1):
    last_result = MIRA.step(
        positions=positions,
        velocities=velocities,
        dt=dt,
        box_length=box_length,
        epsilon=epsilon,
        sigma=sigma,
        integrator=integrator,
    )

    positions = last_result["positions"]
    velocities = last_result["velocities"]

    energy_history.append(
        (
            step_index,
            last_result["potential_energy"],
            last_result["kinetic_energy"],
            last_result["total_energy"],
        )
    )

    if step_index == 1 or step_index % print_every == 0 or step_index == n_steps:
        print(
            f"step {step_index:6d} | "
            f"PE = {last_result['potential_energy']: .8e} | "
            f"KE = {last_result['kinetic_energy']: .8e} | "
            f"E = {last_result['total_energy']: .8e}"
        )

print("\nFinal system")
print("positions =", positions)
print("velocities =", velocities)
print("integrator used =", last_result["integrator"])
print("final PE =", last_result["potential_energy"])
print("final KE =", last_result["kinetic_energy"])
print("final total energy =", last_result["total_energy"])
