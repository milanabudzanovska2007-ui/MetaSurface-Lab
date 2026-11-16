
# META-USAGE

## Application
MetaSurface Lab — scalar optical field simulation and phase quantization evaluation tool.

## Launch
Run:

    python metasurface_lab.py

## Main workflow

1. Select a phase profile.
2. Set the wavelength, grid spacing and grid size.
3. Set the focal length and propagation distance.
4. Set the Gaussian beam radius.
5. Specify the number of quantization levels L.
6. Enable or disable the finite aperture.
7. Click "Порівняти поля та обчислити ІБВП".
8. Review the phase map, reference intensity and quantized intensity.
9. Review E_i, E_phi, E_eta and D_B.

## Batch experiment

Click "Експеримент L=2…256 → CSV" to calculate results for:

2, 4, 8, 16, 32, 64, 128, 256.

The results are saved as a CSV file.

## Image export

Click "Зберегти карти PNG" to export:

- phase_quantized.png
- intensity_reference.png
- intensity_quantized.png

## Validation

Run:

    python metasurface_lab.py --self-test

The built-in tests check identical fields, global phase invariance, zero fields and repeatability of calculations.
