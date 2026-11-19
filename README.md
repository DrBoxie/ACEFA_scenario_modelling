# Reproducing the results of [paper title]

This repository contains the code and input data needed to reproduce the results presented in [paper title and link].

## Software requirements

The code was tested with the following versions.

- Python [version]
- R [version]
- RStudio [version]

### Python packages

- [Add required Python packages and versions]

### R packages

- [Add required R packages and versions]

## Repository structure

The repository root should contain `config.json` and the following directories.

| Path | Contents |
| --- | --- |
| `config.json` | Configuration settings used by the pipeline |
| `python/` | Python scripts |
| `R_code/` | R scripts |
| `data/` | Input data required by the pipeline |

## Configuration

Before running the scripts, edit `config.json` in the repository root as needed. Use lowercase `true` and `false` for Boolean values in JSON.

| Setting | Description | Default |
| --- | --- | --- |
| `paral` | Set to `false` to run the Python code on a single core, or `true` to run it in parallel using the number of cores specified by `nr_cores`. | `true` |
| `nr_cores` | Number of cores to use when `paral` is `true`. | [Add default] |
| `nr_runs_per_part` | Number of simulation runs for each particle in the initial approximate Bayesian computation (ABC) step. | `10` |
| `accept_perc` | Minimum fraction of those runs that must meet the acceptance criteria for a particle to be accepted. Enter this as a fraction between 0 and 1. | `0.9` |
| `vaccination_term` | Timing of the live attenuated influenza vaccine (LAIV) rollout. Use `"term_1"` for simulation days 59–100 or `"term_2"` for simulation days 110–170. | `"term_1"` |
| `K_beta` | Non-negative number that controls the shape of seasonal forcing in transmissibility. A value of `0` gives sinusoidal forcing. As `K_beta` increases, the forcing approaches a square wave. | `0` |
| `nr_particles` | Number of particles evaluated in the initial ABC step. | `70000` |
| `ROOT` | Root directory for outputs. Use `"default"` to save outputs under the repository root, or specify a directory path. Subsequent steps read intermediate outputs from this location. | `"default"` |
| `ABC_output_dir` | Name of the directory where the ABC outputs are saved. | `"ABC outputs"` |
| `nr_proj_per_part` | Number of forward projection runs for each accepted particle. | `50` |
| `Forw_proj_output_dir` | Name of the directory where the forward projection outputs are saved. | `"Forward projection"` |

For example, with `nr_runs_per_part` set to `10` and `accept_perc` set to `0.9`, a particle is accepted if at least 9 of its 10 simulation runs meet the acceptance criteria.

## Running the pipeline

Run these commands from the repository root in the order shown.

1. Run the ABC calibration.

   ```bash
   python python/abc_code.py
   ```

2. Once calibration has finished, run the forward projections.

   ```bash
   python python/forward_proj.py
   ```

Instructions for the remaining pipeline steps will be added as those steps are finalised.

# ACEFA Scenario Modelling

This repository contains the code for the ACEFA scenario-modelling analysis.
Tested with Python 3.10.12.

## Execution pipeline

The full execution pipeline begins with `abc_code.py`.

Before running the ABC analysis, review and adapt the following variables as required:

| Variable | Description |
|---|---|
| `ABC_output_results` | Activates or suppresses output of data frames and plots from the ABC runs. |
| `seed` | Seed used to generate the seed list supplied to the individual simulation runs. |
| `seed2` | Secondary random seed. |
| `parallel` | Determines whether simulations are run in parallel (`True`) or sequentially (`False`). |
| `nr_cores` | Number of CPU cores used for parallel processing. Ignored when `parallel = False`. |
| `nr_particles` | Number of ABC particles to run. |
| `peak_time_range` | Acceptance range for the timing of the epidemic peak. |
| `AR_range` | Acceptance range for the attack rate (AR). |
| `VE_range` | Acceptance range for vaccine efficacy (VE). |

## ABC parameter sampling

The sampling ranges for the ABC particles are defined in `list_to_sample_params`.

The following parameter ranges can be specified:

- `beta_0_range`
- `beta_1_range`
- `shift`
- `half-life`
- `phi_S1_range`
- `phi_V0_range`
- `phi_V1_range`
- `frac_S1_range`

These ranges determine the parameter space from which the ABC particles are sampled.

## Model parameters

Core model variables are defined in the `init_params` function in `main_sim.py`. If the underlying model structure or assumptions are changed, the corresponding parameters should be updated there.

| Variable | Description |
|---|---|
| `tau` | Size of the tau-leap time step. |
| `t_start`, `t_end` | Start and end days of the simulation. |
| `age_groups` | Age groups included in the model. |
| `compartments` | Epidemiological and vaccination compartments included in the model. |
| `tot_pop` | Total population size. |
| `frac_ages` | Fraction of the population belonging to each age group. |
| `C` | Age-stratified contact matrix. |
| `gamma` | Recovery rate for each age group. |
| `target_cover` | Target vaccination coverage for each age group. |
| `iiv_rates` | Weekly cumulative rollout percentages for inactivated influenza vaccine (IIV). |
| `laiv_rates` | Weekly cumulative rollout percentages for live attenuated influenza vaccine (LAIV). |
| `p_ext` | External force of infection. |

## Recommended workflow

1. Review the ABC settings and acceptance criteria in `abc_code.py`.
2. Review the ABC parameter sampling ranges in `list_to_sample_params`.
3. If the model structure or fixed parameters have changed, update `init_params` in `main_sim.py`.
4. Run `abc_code.py`.
5. Inspect the accepted ABC particles and generated outputs before proceeding to subsequent scenario analyses.
