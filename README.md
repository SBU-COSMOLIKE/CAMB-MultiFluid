# CAMB-MultifluidDE

Authors: Vivian Miranda, João Victor Rebouças, Diogo H. F. de Souza

## Description:

This is a CAMB modification that makes it easier to add multiple components to dark energy.

## Usage with Python Interface:

The basic usage is presented in `docs/Multifluid_demo.ipynb`.

To use a model with multiple fluids, select `dark_energy_model = 'MultiFluidDE'`. The parameter `num_of_components` chooses the number of dark energy components you want to use (so far only accepts 2 components, but this template is extensible!). The models of each component are set in the list parameter `models: list[int] = [<model1>, <model2>, <model3>, <model4>]` (the list needs to have 4 components, the last two ones can be arbitrary integers).

The model choices are:
- Late-time models are set in the first component of `models`: 1 - constant w; 2 - w0wa.
- Early-time models are set in the second component of `models`: 1 - fluid EDE; 2 - Axion-like scalar field.

Each model has parameters you need to input in the `set_params` function:
- `w0` and `wa` for late DE
- `zc`, `fde_zc`, `theta_i` and `wn` for fluid EDE or scalar field EDE
- `use_zc`, `initial_phi`, `m`, `f`, `V0` for scalar field EDE

## TODO:
- Nonlinear model
- Python interface