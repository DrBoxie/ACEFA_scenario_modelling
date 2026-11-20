# ACEFA scenario modelling hub for LAIV school-based rollout program: University of Melbourne contribution

This repository contains the code and publicly available input data for reproducing the results presented in [Modelling the impact of a school-based influenza vaccination programme under uncertainty in seasonal transmission dynamics](https://doi.org/YOUR-DOI). This is part of the ACEFA scenario modelling hub.

## Software requirements

The code was tested with the following versions.

- Python 3.10.12
- R 4.5.3
- RStudio Server 2026.1.1.403

### Python packages

The Python package versions used in the analyses are listed in
[requirements.txt](requirements.txt).

With your Python environment activated, run the following command from
the repository root:

```bash
python -m pip install -r requirements.txt
```

### R packages

- [Add required R packages and versions]

## Repository structure

The repository root should contain the following:

| Path | Contents |
| --- | --- |
| `config.json` | Configuration settings used by the pipeline |
| `requirements.txt` | Python dependencies and the versions used |
| `python/` | Python scripts |
| `R_code/` | R scripts |
| `data/` | Publicly available input data |

## Data availability

Aggregate hospitalisation data from the Australian Institute of Health and Welfare (AIHW),
population data by jurisdiction from the Australian Bureau of Statistics (ABS), and a CSV file
containing influenza epidemic peak timings are publicly available in this repository.

The peak timings were derived from influenza case data from the Australian National Notifiable
Disease Surveillance System (NNDSS). The underlying case data cannot be shared under the terms of
our data access agreements. Requests for access should be submitted to
[NNDSS.datarequests@health.gov.au](mailto:NNDSS.datarequests@health.gov.au).
A data committee independent of the authorship group assesses these requests.

## Configuration

Before running the scripts, edit `config.json` in the repository root as needed. Use lowercase `true` and `false` for Boolean values in JSON.

| Setting | Description | Default |
| --- | --- | --- |
| `paral` | Set to `false` to run the Python code on a single core, or `true` to run it in parallel using the number of cores specified by `nr_cores`. | `true` |
| `nr_cores` | Number of cores to use when `paral` is `true`. | `12` |
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

After installing the Python dependencies and reviewing config.json, run these commands from the repository root in the order shown.

1. Run the ABC calibration.

   ```bash
   python python/abc_code.py
   ```

2. Once calibration has finished and the Python process has exited, run the forward projections as a separate command.

   ```bash
   python python/forward_proj.py
   ```

Instructions for the remaining pipeline steps will be added as those steps are finalised.

