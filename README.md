
# MetaSurface Lab

MetaSurface Lab is a Python application for numerical simulation of phase quantization in optical metasurfaces and multi-criteria evaluation of optical fields.

## Features

- Generation of five phase profiles: lens, axicon, vortex, linear tilt and random phase.
- Simulation of optical field propagation using the angular spectrum method.
- Comparison of continuous and quantized phase profiles.
- Calculation of the Multi-Criteria Optical Field Deviation Index (MOFDI / ІБВП).
- Visualization of phase and intensity distributions.
- Export of experimental results to CSV and image maps to PNG.
- Automated quantization sweep for L = 2, 4, 8, 16, 32, 64, 128, 256.

## Requirements

- Python 3.10 or newer
- NumPy
- SciPy
- Pillow
- Tkinter

## Installation

Clone the repository:

```bash
git clone <repository-url>
cd MetaSurface-Lab
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Run the application

```bash
python metasurface_lab.py
```

## Run validation tests

```bash
python metasurface_lab.py --self-test
```

## Run an experiment

```bash
python metasurface_lab.py --sweep results/ibvp_quantization_sweep.csv
```

## MOFDI formula

The proposed index combines three components:

D_B = 1 - (1 - E_i)(1 - E_phi)(1 - E_eta)

Where:

- E_i — normalized intensity discrepancy;
- E_phi — weighted phase inconsistency;
- E_eta — difference in energy concentration within the target region.

## Model limitations

The application uses scalar optical field propagation based on the angular spectrum method. It does not model individual nanostructures or solve the full vector Maxwell equations.

## Author

Milana Budzanovska  
Faculty of Information Technology  
Taras Shevchenko National University of Kyiv
