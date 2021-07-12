# lift_coords

A basic Python wrapper for LiftOver.

## Usage

```python
import lift_coords
df2, failed = lift_coords.lift_over(df, 'grch38', 'grch37', keep_orig=True)
```

## Installation

```bash
pip install git+https://github.com/Townsend-Lab-Yale/lift_coords.git
```


<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.0.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
