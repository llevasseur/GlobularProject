------------------------------------------------------------------------
FINAL REPORT - CSC 473 - LEEVON LEVASSEUR V00868326 Dec 6/2021
------------------------------------------------------------------------

Edit the file start.tcl in Build directory to customize parameters for different animation results.

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