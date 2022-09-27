# Globular Dynamics Animation Project
This project demonstrates a method for modelling and animating viscous fluid using a physically-based particle simulator and soft collisions. This project was inspired by the work presented by Gavin Miller & Andrew Pearse in Globular Dynamics: A connected particle system for animating viscous fluids. Using the same fundamental theory, this project delves into the decisions for variable values and results that Miller & Pearse neglect. Results are expanded to show how each variable influences the animation. Specifically, how to make certain materials related to fluid. The simulation was constructed using C++ and rendered using OpenGL and AnimTcl. A Forward Euler approach with a small time step and elastic penalty collisions were used to determine the next position of each particle. The results of this project are customizable by the user to form material that looks like powder, liquid, or solid.

## Usage
Open this file in Visual Studio and select `cs174-2010.sln` to access the source code. Then, run `Windows Local Debugger` to start the TKCon simulation window.  Use the `s` key to start or pause the simulation.

Edit the file `start.tcl` in the /Build directory to customize parameters for different animation results. Below are the customizable variables:

Set number of globules by changing variable "num".

Set r0 using:

	system globSys r0 <x>

Set constants m, n, b1, timestep:

	simulator globSim constants <m> <n> <b1> <timestep>

Set material type:

	simulator globSim material <type>

	<type> includes: powder, liquid, solid

Set ground spring constants:

	simulator globSim ground <ks> <kd>

Set global gravity:

	simulator globSim gravity <dy>

To test different viscosity for liquid, change "alpha" in GlobuleSim::Cd.

## Contributing
Bug reports are welcome on Github at [Issues](https://github.com/llevasseur/GlobularProject/issues).

## Future Work
1. Test globules with additional effectors, including obstacles, containers, and external forces. Correct collision detection must be implemented which will slow the simulation.
2. Speed up simulation calculations. This can be done by switching to a semi-lagrangian method for calculating the next position of a particle, and by offloading computations from the CPU to the GPU.
3. Add temperature to the simulation to a time variance to the material type. Simulate freezing and melting of a material.

### Keywords:
Globular Dynamics, Particle Simulator, C++, OpenGL, AnimTcl, Eulerian collisions.
