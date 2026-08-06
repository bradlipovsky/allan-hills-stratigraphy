"""Boundary-generated recumbent folds in steady full-Stokes flow."""

from pathlib import Path
import logging
import warnings

warnings.filterwarnings("ignore", message="Unable to import recommended hash")
warnings.filterwarnings("ignore", message="The .*Function.at.* method")
logging.getLogger("tsfc").setLevel(logging.ERROR)

import firedrake as fd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from IPython.display import display
from scipy.integrate import solve_ivp


LENGTH = 800.0
THICKNESS = 160.0
DENSITY = 917.0
GRAVITY = 9.81
SURFACE_SLOPE = 0.03
VISCOSITY = 2.5e8

CASES = [
    {
        "name": "persistent extrusion",
        "amplitude": 150.0,
        "width_left": 300.0,
        "width_right": 80.0,
        "basal_speed": 0.20,
        "source_limits": (400.0, 520.0),
        "ages": (0.0, 1000.0, 2000.0, 2500.0, 3000.0),
        "target_age": 2500.0,
        "plot_limits": (350.0, 1250.0),
    },
    {
        "name": "slow-flow amplification",
        "amplitude": 140.0,
        "width_left": 45.0,
        "width_right": 260.0,
        "basal_speed": 0.20,
        "source_limits": (140.0, 220.0),
        "ages": (0.0, 200.0, 400.0, 600.0, 800.0, 900.0),
        "target_age": 800.0,
        "plot_limits": (100.0, 650.0),
    },
]

MESHES = [(64, 14), (96, 21), (128, 28)]


def bed_profile(x, case):
    """Raised-cosine bed obstacle, repeated with the model period."""
    distance = (np.asarray(x) - 0.5 * LENGTH + 0.5 * LENGTH) % LENGTH - 0.5 * LENGTH
    bed = np.zeros_like(distance, dtype=float)
    left = (distance >= -case["width_left"]) & (distance <= 0.0)
    right = (distance > 0.0) & (distance <= case["width_right"])
    bed[left] = 0.5 * case["amplitude"] * (
        1.0 + np.cos(np.pi * distance[left] / case["width_left"])
    )
    bed[right] = 0.5 * case["amplitude"] * (
        1.0 + np.cos(np.pi * distance[right] / case["width_right"])
    )
    return bed


def represented_bed(x, case, nx):
    """Piecewise-linear bed represented by a mesh with nx divisions."""
    mesh_x = np.linspace(0.0, LENGTH, nx, endpoint=False)
    mesh_bed = bed_profile(mesh_x, case)
    return np.interp(
        np.mod(x, LENGTH),
        np.r_[mesh_x, LENGTH],
        np.r_[mesh_bed, mesh_bed[0]],
    )


def bed_slope_expression(x, case):
    """Analytical bed slope used in the no-penetration sliding condition."""
    distance = x - 0.5 * LENGTH
    left = fd.And(distance >= -case["width_left"], distance <= 0.0)
    right = fd.And(distance > 0.0, distance <= case["width_right"])
    left_slope = -0.5 * case["amplitude"] * np.pi / case["width_left"] * fd.sin(
        np.pi * distance / case["width_left"]
    )
    right_slope = -0.5 * case["amplitude"] * np.pi / case["width_right"] * fd.sin(
        np.pi * distance / case["width_right"]
    )
    return fd.conditional(left, left_slope, fd.conditional(right, right_slope, 0.0))


def solve_stokes(case, nx, nz):
    """Solve incompressible, linear-viscous full Stokes flow."""
    mesh = fd.PeriodicRectangleMesh(nx, nz, LENGTH, THICKNESS, direction="x")
    reference = mesh.coordinates.dat.data_ro.copy()
    sigma = reference[:, 1] / THICKNESS
    bed = bed_profile(reference[:, 0], case)
    mesh.coordinates.dat.data[:, 1] = bed + sigma * (THICKNESS - bed)
    mesh.clear_spatial_index()

    velocity_space = fd.VectorFunctionSpace(mesh, "CG", 2)
    pressure_space = fd.FunctionSpace(mesh, "CG", 1)
    mixed_space = velocity_space * pressure_space
    solution = fd.Function(mixed_space)
    velocity_symbol, pressure_symbol = fd.split(solution)
    test_velocity, test_pressure = fd.TestFunctions(mixed_space)
    body_force = fd.Constant(
        (
            DENSITY * GRAVITY * np.sin(SURFACE_SLOPE),
            -DENSITY * GRAVITY * np.cos(SURFACE_SLOPE),
        )
    )
    strain_rate = fd.sym(fd.grad(velocity_symbol))
    residual = (
        2.0 * VISCOSITY * fd.inner(strain_rate, fd.sym(fd.grad(test_velocity)))
        - pressure_symbol * fd.div(test_velocity)
        - test_pressure * fd.div(velocity_symbol)
        - fd.dot(body_force, test_velocity)
    ) * fd.dx

    x = fd.SpatialCoordinate(mesh)[0]
    basal_velocity = fd.as_vector(
        (
            case["basal_speed"],
            case["basal_speed"] * bed_slope_expression(x, case),
        )
    )
    bed_condition = fd.DirichletBC(mixed_space.sub(0), basal_velocity, 1)
    parameters = {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    solution.sub(0).interpolate(fd.Constant((case["basal_speed"], 0.0)))
    fd.solve(residual == 0, solution, bcs=[bed_condition], solver_parameters=parameters)
    velocity, pressure = solution.subfunctions
    pressure.dat.data[:] -= np.mean(pressure.dat.data_ro)
    gradient_space = fd.TensorFunctionSpace(mesh, "DG", 1)
    velocity_gradient = fd.Function(gradient_space).interpolate(fd.grad(velocity))
    return mesh, velocity, velocity_gradient


def make_sampler(case, nx, velocity, velocity_gradient, include_gradient=False):
    """Resample the finite-element solution for trajectory integration."""
    horizontal_samples = 201
    vertical_samples = 80
    x = np.linspace(0.0, LENGTH, horizontal_samples, endpoint=False)
    x_evaluation = x.copy()
    x_evaluation[0] = 1.0e-6
    sigma = np.linspace(0.002, 0.99999, vertical_samples)
    bed = represented_bed(x, case, nx)
    X = np.tile(x_evaluation, (len(sigma), 1))
    Z = bed[None, :] + sigma[:, None] * (THICKNESS - bed)[None, :]
    points = np.column_stack((X.ravel(), Z.ravel()))
    velocity_values = np.asarray(
        velocity.at(points, tolerance=1.0e-7)
    ).reshape(Z.shape + (2,))
    fields = {
        "u": velocity_values[..., 0],
        "w": velocity_values[..., 1],
    }
    if include_gradient:
        gradient_values = np.asarray(
            velocity_gradient.at(points, tolerance=1.0e-7)
        ).reshape(Z.shape + (2, 2))
        fields["uz"] = gradient_values[..., 0, 1]
    dx = LENGTH / len(x)
    dsigma = sigma[1] - sigma[0]

    def sample(x_position, z_position, names=("u", "w")):
        x_wrapped = np.mod(x_position, LENGTH)
        local_bed = float(represented_bed(np.array([x_wrapped]), case, nx)[0])
        local_sigma = (z_position - local_bed) / (THICKNESS - local_bed)
        if not names:
            return {"sigma": local_sigma}
        x_index = x_wrapped / dx
        ix0 = int(np.floor(x_index)) % len(x)
        ix1 = (ix0 + 1) % len(x)
        x_weight = x_index - np.floor(x_index)
        sigma_index = (local_sigma - sigma[0]) / dsigma
        iz0 = int(np.clip(np.floor(sigma_index), 0, len(sigma) - 2))
        iz1 = iz0 + 1
        sigma_weight = np.clip(sigma_index - iz0, 0.0, 1.0)
        output = {}
        for name in names:
            values = fields[name]
            lower = (1.0 - x_weight) * values[iz0, ix0] + x_weight * values[iz0, ix1]
            upper = (1.0 - x_weight) * values[iz1, ix0] + x_weight * values[iz1, ix1]
            output[name] = float((1.0 - sigma_weight) * lower + sigma_weight * upper)
        return output | {"sigma": local_sigma}

    return sample


def integrate_path(sample, case, nx, x0, maximum_age, stop_at_surface=True):
    """Trace one particle from its surface deposition position."""
    start_sigma = 0.9999
    start_bed = float(represented_bed(np.array([x0]), case, nx)[0])
    z0 = start_bed + start_sigma * (THICKNESS - start_bed)
    surface = sample(x0, z0)
    accumulation = -surface["w"]
    surface_speed = surface["u"]

    def right_hand_side(time, state):
        values = sample(state[0], state[1], names=("u", "w"))
        return [values["u"], values["w"]]

    def surface_return(time, state):
        if time < 1.0:
            return -1.0
        return sample(state[0], state[1], names=())["sigma"] - 0.99995

    surface_return.terminal = True
    surface_return.direction = 1

    def one_way_limit(time, state):
        return sample(state[0], state[1], names=("u",))["u"] - 1.0e-6 * case["basal_speed"]

    one_way_limit.terminal = True
    one_way_limit.direction = -1
    events = [surface_return, one_way_limit] if stop_at_surface else None
    path = solve_ivp(
        right_hand_side,
        (0.0, maximum_age),
        [x0, z0],
        rtol=1.0e-7,
        atol=1.0e-9,
        max_step=10.0,
        events=events,
        dense_output=True,
    )
    return path, accumulation, surface_speed


def trace_sources(sample, case, nx, source_x):
    """Build the steady material map from surface position and age."""
    maximum_age = max(case["ages"])
    traces = []
    for x0 in source_x:
        path, accumulation, surface_speed = integrate_path(
            sample, case, nx, x0, maximum_age
        )
        traces.append(
            {
                "x0": x0,
                "path": path,
                "accumulation": accumulation,
                "surface_speed": surface_speed,
            }
        )
    return traces


def material_curve(traces, age):
    """Return one boundary-generated isochrone and its source labels."""
    source_x = np.array([trace["x0"] for trace in traces])
    positions = np.full((len(traces), 2), np.nan)
    for index, trace in enumerate(traces):
        if trace["accumulation"] > 0.0 and age <= trace["path"].t[-1]:
            positions[index] = trace["path"].sol(age)
    valid = np.isfinite(positions[:, 0])
    return source_x, positions, valid


def order_from_material_map(traces, age):
    """Evaluate dx/dx0 directly from neighboring surface trajectories."""
    source_x, positions, valid = material_curve(traces, age)
    derivative = np.gradient(positions[valid, 0], source_x[valid])
    valid_indices = np.where(valid)[0]
    return source_x, positions, valid, valid_indices, derivative


def integrate_budget(sample, trace, age):
    """Integrate the positive and reverse vertical-shear budgets."""
    times = np.linspace(0.0, age, 501)
    positions = trace["path"].sol(times)
    values = [sample(x, z, names=("u", "uz")) for x, z in positions.T]
    speed = np.array([value["u"] for value in values])
    vertical_shear = np.array([value["uz"] for value in values])
    source_factor = trace["accumulation"] * trace["surface_speed"]
    positive_rate = source_factor * np.maximum(vertical_shear, 0.0) / speed**2
    reverse_rate = source_factor * np.maximum(-vertical_shear, 0.0) / speed**2
    dt = np.diff(times)
    positive = np.r_[0.0, np.cumsum(0.5 * (positive_rate[1:] + positive_rate[:-1]) * dt)]
    reverse = np.r_[0.0, np.cumsum(0.5 * (reverse_rate[1:] + reverse_rate[:-1]) * dt)]
    return {
        "time": times,
        "speed": speed,
        "positive": positive,
        "reverse": reverse,
        "order_budget": 1.0 + positive - reverse,
    }


def consecutive_groups(indices):
    if len(indices) == 0:
        return []
    return np.split(indices, np.where(np.diff(indices) > 1)[0] + 1)


def output_directory():
    directory = Path("figures") if Path("figures").exists() else Path("../figures")
    directory.mkdir(exist_ok=True)
    return directory


solutions = {}
convergence_rows = []
print("Solving the two steady full-Stokes examples.")
for case in CASES:
    source_x = np.linspace(*case["source_limits"], 81)
    for nx, nz in MESHES:
        print(f"  {case['name']}, mesh {nx} x {nz}", flush=True)
        mesh, velocity, velocity_gradient = solve_stokes(case, nx, nz)
        sample = make_sampler(case, nx, velocity, velocity_gradient)
        traces = trace_sources(sample, case, nx, source_x)
        _, _, _, _, derivative = order_from_material_map(traces, case["target_age"])
        convergence_rows.append(
            {
                "case": case["name"],
                "mesh": f"{nx} x {nz}",
                "minimum dx/dx0": derivative.min(),
                "overturned source intervals": np.count_nonzero(derivative < 0.0),
            }
        )
        solutions[(case["name"], nx, nz)] = {
            "sample": sample,
            "traces": traces,
            "velocity": velocity if (nx, nz) == MESHES[-1] else None,
            "velocity_gradient": velocity_gradient if (nx, nz) == MESHES[-1] else None,
        }

convergence = pd.DataFrame(convergence_rows)
display(convergence.round(4))


figure_directory = output_directory()
fine_nx, fine_nz = MESHES[-1]
for case in CASES:
    result = solutions[(case["name"], fine_nx, fine_nz)]
    result["sample"] = make_sampler(
        case,
        fine_nx,
        result["velocity"],
        result["velocity_gradient"],
        include_gradient=True,
    )
fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.6), constrained_layout=True)
letters = iter("abcd")
selected = {}
for row, case in enumerate(CASES):
    result = solutions[(case["name"], fine_nx, fine_nz)]
    traces = result["traces"]
    ax = axes[row, 0]
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(case["ages"])))
    curves = {}
    for age, color in zip(case["ages"], colors):
        source_x, positions, valid = material_curve(traces, age)
        curves[age] = (source_x, positions, valid)
        for group in consecutive_groups(np.where(valid)[0]):
            if len(group) > 1:
                ax.plot(
                    positions[group, 0],
                    positions[group, 1],
                    color=color,
                    linewidth=1.6,
                    label=f"{age / 1000:g} kyr" if group[0] == np.where(valid)[0][0] else None,
                )
    bed_x = np.linspace(0.0, 2.0 * LENGTH, 1001)
    ax.fill_between(bed_x, 0.0, bed_profile(np.mod(bed_x, LENGTH), case), color="0.7")
    ax.axhline(THICKNESS, color="black", linewidth=0.8)
    ax.set(
        xlim=case["plot_limits"],
        ylim=(0.0, THICKNESS),
        xlabel="Unwrapped horizontal coordinate (m)" if row == 1 else None,
        ylabel="Elevation (m)",
        title=f"({next(letters)}) {case['name']}",
    )
    ax.legend(frameon=False, ncol=3, fontsize=8, loc="lower left")

    zoom = axes[row, 1]
    source_x, positions, valid, valid_indices, derivative = order_from_material_map(
        traces, case["target_age"]
    )
    negative_local = np.where(derivative < 0.0)[0]
    if len(negative_local) == 0:
        raise RuntimeError("A refined material map did not contain the expected fold.")
    start = max(0, negative_local[0] - 7)
    stop = min(len(valid_indices), negative_local[-1] + 8)
    indices = valid_indices[start:stop]
    zoom.plot(positions[indices, 0], positions[indices, 1], color="black", linewidth=1.1)
    points = zoom.scatter(
        positions[indices, 0],
        positions[indices, 1],
        c=source_x[indices],
        cmap="viridis",
        norm=Normalize(source_x[indices].min(), source_x[indices].max()),
        s=21,
        zorder=3,
    )
    negative_indices = valid_indices[negative_local]
    zoom.scatter(
        positions[negative_indices, 0],
        positions[negative_indices, 1],
        facecolors="none",
        edgecolors="red",
        linewidths=0.9,
        s=42,
        label=r"$\partial x/\partial x_0<0$",
    )
    x_pad = max(1.0, 0.12 * np.ptp(positions[indices, 0]))
    y_pad = max(1.0, 0.10 * np.ptp(positions[indices, 1]))
    zoom.set(
        xlim=(positions[indices, 0].min() - x_pad, positions[indices, 0].max() + x_pad),
        ylim=(positions[indices, 1].min() - y_pad, positions[indices, 1].max() + y_pad),
        xlabel="Horizontal coordinate (m)",
        ylabel="Elevation (m)",
        title=f"({next(letters)}) {case['target_age'] / 1000:g} kyr fold detail",
    )
    zoom.legend(frameon=False, fontsize=8)
    fig.colorbar(points, ax=zoom, pad=0.01, label=r"Deposition coordinate $x_0$ (m)")
    selected_local = int(np.argmin(derivative))
    selected_index = valid_indices[selected_local]
    selected_trace = traces[selected_index]
    selected[case["name"]] = {
        "trace": selected_trace,
        "material_jacobian": derivative[selected_local],
        "source_x": source_x[selected_index],
        "sample": result["sample"],
    }

fold_figure = figure_directory / "steady_full_stokes_recumbent_folds.png"
fig.savefig(fold_figure, dpi=220)
plt.show()


fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.2), constrained_layout=True)
diagnostic_rows = []
for row, case in enumerate(CASES):
    choice = selected[case["name"]]
    budget = integrate_budget(choice["sample"], choice["trace"], case["target_age"])
    order = budget["order_budget"]
    crossings = np.where(np.diff(np.signbit(order)))[0]
    fold_age = budget["time"][crossings[0]] if len(crossings) else np.nan
    ax = axes[row, 0]
    ax.plot(budget["time"] / 1000.0, budget["positive"], label=r"$\mathcal{P}$")
    ax.plot(budget["time"] / 1000.0, budget["reverse"], label=r"$\mathcal{R}$")
    ax.plot(budget["time"] / 1000.0, order, color="black", label=r"$1+\mathcal{P}-\mathcal{R}$")
    ax.axhline(0.0, color="0.5", linewidth=0.8)
    if np.isfinite(fold_age):
        ax.axvline(fold_age / 1000.0, color="0.5", linestyle=":", linewidth=0.9)
    ax.set(
        xlabel="Age (kyr)",
        ylabel="Order-budget terms",
        title=f"({chr(97 + 2 * row)}) {case['name']}: shear budget",
    )
    ax.legend(frameon=False, fontsize=8)

    ax = axes[row, 1]
    speed_ratio = budget["speed"] / choice["trace"]["surface_speed"]
    ax.plot(budget["time"] / 1000.0, speed_ratio, color="black")
    ax.axhline(1.0, color="0.5", linewidth=0.8)
    if np.isfinite(fold_age):
        ax.axvline(fold_age / 1000.0, color="0.5", linestyle=":", linewidth=0.9)
    ax.set(
        xlabel="Age (kyr)",
        ylabel=r"Horizontal speed $u/u_{s0}$",
        title=f"({chr(98 + 2 * row)}) speed along the same trajectory",
    )
    diagnostic_rows.append(
        {
            "case": case["name"],
            "source x0 (m)": choice["source_x"],
            "minimum u/us0": speed_ratio.min(),
            "minimum material dx/dx0": choice["material_jacobian"],
            "final order budget": order[-1],
            "first budget crossing (yr)": fold_age,
        }
    )

mechanism_figure = figure_directory / "steady_full_stokes_fold_mechanisms.png"
fig.savefig(mechanism_figure, dpi=220)
plt.show()

diagnostics = pd.DataFrame(diagnostic_rows)
display(diagnostics.round(4))


fig, ax = plt.subplots(figsize=(5.4, 3.8), constrained_layout=True)
for case in CASES:
    group = convergence[convergence["case"] == case["name"]]
    vertical_divisions = np.array([int(mesh.split(" x ")[1]) for mesh in group["mesh"]])
    ax.plot(
        vertical_divisions,
        group["minimum dx/dx0"],
        marker="o",
        label=case["name"],
    )
ax.axhline(0.0, color="0.5", linewidth=0.8)
ax.set(
    xlabel="Vertical mesh divisions",
    ylabel=r"Minimum $\partial x/\partial x_0$",
    title="Fold sign under mesh refinement",
)
ax.legend(frameon=False, fontsize=8)
convergence_figure = figure_directory / "steady_full_stokes_fold_convergence.png"
fig.savefig(convergence_figure, dpi=220)
plt.show()

print("Saved figures:")
for path in (fold_figure, mechanism_figure, convergence_figure):
    print(f"  {path}")
