# CUDA kernels for Particles simulation

## Description
Simple simulation of particles in a gravitational field using CUDA.

## Authors:
- Pietro Portolani, IT PhD candidate @ Politecnico di Milano
- Tommaso Alfonsi, IT PhD candidate @ Politecnico di Milano

## [Slides (italian)](https://github.com/piepor/particles-simulation-parallel/blob/9aba280d10bd8cc36bb56d3a3e9b3f10f0379fd3/presentazione/presentazione.pptx)

## Code versions:
- [Original source code](https://github.com/piepor/particles-simulation-parallel/blob/9aba280d10bd8cc36bb56d3a3e9b3f10f0379fd3/original_code/07-Particles2D/07-Particles2D_c/particles_c.c)
- [GPU accelerated version (final)](https://github.com/piepor/particles-simulation-parallel/blob/9aba280d10bd8cc36bb56d3a3e9b3f10f0379fd3/particles_first_opt.cu)
- Other accelerated versions:
    * [ForceCompt_par_v2](https://github.com/piepor/particles-simulation-parallel/blob/9aba280d10bd8cc36bb56d3a3e9b3f10f0379fd3/particles_second_opt_not_enough.cu)
    * [ForceCompt_par_v3](https://github.com/piepor/particles-simulation-parallel/blob/9aba280d10bd8cc36bb56d3a3e9b3f10f0379fd3/particles_second_opt_two_streams.cu)
    * [Acceleration of ParticleScreen](https://github.com/piepor/particles-simulation-parallel/blob/9aba280d10bd8cc36bb56d3a3e9b3f10f0379fd3/particles_third_opt.cu)
