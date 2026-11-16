
# Algorithm Description

## Processing stages

1. Read model parameters.
2. Generate spatial coordinate grid.
3. Generate Gaussian input beam.
4. Generate continuous phase profile.
5. Quantize the phase profile.
6. Propagate both fields using the angular spectrum method.
7. Calculate intensity distributions.
8. Calculate E_i, E_phi and E_eta.
9. Calculate the overall MOFDI index.
10. Display and export results.

## Algorithm flowchart

```mermaid
flowchart TD
    A[Start] --> B[Read parameters]
    B --> C[Generate coordinate grid]
    C --> D[Generate Gaussian input field]
    D --> E[Generate continuous phase profile]
    E --> F[Quantize phase]
    F --> G[Propagate reference and test fields]
    G --> H[Calculate intensity distributions]
    H --> I[Calculate E_i, E_phi, E_eta]
    I --> J[Calculate MOFDI]
    J --> K[Display results]
    K --> L[Export CSV or PNG]
    L --> M[End]
```

## Main algorithms

- Phase profile generation: `make_phase_profile()`
- Phase quantization: `quantize_phase()`
- Angular spectrum propagation: `angular_spectrum()`
- Index calculation: `calculate_index()`
- Simulation: `simulate_pair()`
- Batch experiment: `run_quantization_sweep()`
