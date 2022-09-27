#include "GlobuleSim.h"

GlobuleSim::GlobuleSim(const std::string& name, Globules* target) :
	BaseSimulator(name),
	mySystem(target)
{
}	// GlobuleSim

GlobuleSim::~GlobuleSim()
{
}	// GlobuleSim::~GlobuleSim

double GlobuleSim::Cr(double r0)
{
	if (strcmp(material, "solid") == 0)
	{
		return 2 * r0;
	}
	else // infer powder or liquid (Cr the same for both)
	{
		return r0;
	}
} // GlobuleSim::Cr

double GlobuleSim::Cd(double r0)
{
	int alpha = 1.5; // Change this value to control viscosity of liquid. [1.5, 2.0]
	if (strcmp(material, "powder") == 0)
	{
		return r0;
	}
	else if (strcmp(material, "liquid") == 0)
	{
		return alpha * r0;
	}
	else // infer solid (Cd the same for both)
	{
		return 2 * r0;
	}
} // GlobuleSim::Cd

double GlobuleSim::b2(double r0)
{
	return b1 * pow(r0, n - m);
} // GlobuleSim::b2

void GlobuleSim::solveF(int index, vector<double*> globs, vector<double*> vel, vector<double> radii, double r0, double mass, double* FTotal)
{
	double* Fip = new double[3]{ 0.0f, 0.0f, 0.0f };
	double* pBar = new double[3]{ 0.0f, 0.0f, 0.0f };
	double* vBar = new double[3]{ 0.0f, 0.0f, 0.0f };
	double repulsionTerm;
	double dragTerm;
	double thisCr;
	double thisCd;
	double D;
	double secondTermSr;
	double Sr;
	double secondTermSd;
	double Sd;
	double thisB2;
	double vDotp;
	double influence;

	for (int i = 0; i < globs.size(); i++)
	{
		D = sqrt(pow(globs[i][0] - globs[index][0], 2) + pow(globs[i][1] - globs[index][1], 2) + pow(globs[i][2] - globs[index][2], 2));
		if (i != index)
		{
			Fip[0] = 0.0f;
			Fip[1] = 0.0f;
			Fip[2] = 0.0f;

			repulsionTerm = 0;
			dragTerm = 0;
			thisCr = Cr(r0);
			thisCd = Cd(r0);

			if (D == 0)
			{
				D = 0.001;
			}

			pBar[0] = (globs[index][0] - globs[i][0]);
			pBar[1] = (globs[index][1] - globs[i][1]);
			pBar[2] = (globs[index][2] - globs[i][2]);
			// Vector difference between globules

			// solve for Sr and Sd
			secondTermSr = pow(D, 2.0) / (pow(thisCr, 2.0) * pow(radii[i] + radii[index], 2.0));
			Sr = 1.0 - secondTermSr;

			secondTermSd = pow(D, 2.0) / (pow(thisCd, 2.0) * pow(radii[i] + radii[index], 2.0));
			Sd = 1.0 - secondTermSd;

			// Check for negligibal force, ie the globules are too far away to influence
			if (Sr <= 0.0)
			{
				Sr = 0.0;
			}
			else
			{
				thisB2 = b2(r0);
				repulsionTerm = Sr * ((b1 / pow(D, m)) - (thisB2 / pow(D, n)));

			}
			if (Sd <= 0.0)
			{
				Sd = 0.0;
			}
			else
			{

				vBar[0] = (vel[index][0] - vel[i][0]);
				vBar[1] = (vel[index][1] - vel[i][1]);
				vBar[2] = (vel[index][2] - vel[i][2]);
				// Vector differnece between globule velocities

				// find dot product of vBar and pBar
				vDotp = vBar[0] * pBar[0] + vBar[1] * pBar[1] + vBar[2] * pBar[2];

				dragTerm = Sd * (vDotp / pow(D, 2));
			}

			// Solve for Repulsion/Attraction - Drag
			influence = repulsionTerm - dragTerm;

			
			
			// Finally, solve for Fip
			Fip[0] = pBar[0] * influence;
			Fip[1] = pBar[1] * influence;
			Fip[2] = pBar[2] * influence;


		}
		// Sum forces
		FTotal[0] = FTotal[0] + Fip[0];
		FTotal[1] = FTotal[1] + Fip[1];
		FTotal[2] = FTotal[2] + Fip[2];
	}
	// Environmental forces, drag and gravity
	FTotal[1] = FTotal[1] + (mass * gravity);

} // GlobuleSim::solveF


int GlobuleSim::step(double time)
{
	// Determine which integration method to use to determine new velocity and new position based on set time step
	// retrive data from system
	vector<double*> newGlobs;
	vector<double*> newV;
	vector<double*> globs = mySystem->getGlobules();
	vector<double*> vel = mySystem->getVelocities();
	vector<double*> prevGlobs = mySystem->getPrevGlobules();

	vector<double> masses = mySystem->getMasses();
	vector<double> radii = mySystem->getRadii();
	double r0 = mySystem->getR0();
	// end retrive data

	// initialize constants
	double* groundP = new double[3]{ 0.0, -5.0, 0.0 }; // initialized in system
	double* groundN = new double[3]{ 0.0, 1.0, 0.0 };

	// use user requested time
	time = myTime;
	//animTcl::OutputMessage("MY TIME: %f, MATERIAL: %s", time, material);

	for (int p = 0; p < globs.size(); p++)
	{
		//animTcl::OutputMessage("inSim globs[p][1]: %f, vel[p][1]: %f", globs[p][1], vel[p][1]);
		double* FTotal = new double[3]{ 0.0f, 0.0f, 0.0f };
		solveF(p, globs, vel, radii, r0, masses[p], FTotal);
		//animTcl::OutputMessage("FORCE: %f %f %f", FTotal[0], FTotal[1], FTotal[2]);
		
		// forward euler
		double* newVel = new double[3]{ 0.0f, 0.0f, 0.0f };
		newVel[0] = vel[p][0] + time * FTotal[0] / masses[p];
		newVel[1] = vel[p][1] + time * FTotal[1] / masses[p];
		newVel[2] = vel[p][2] + time * FTotal[2] / masses[p];

		//animTcl::OutputMessage("newVel[1]: %f, globs[p][1]: %f", newVel[1], globs[p][1]);
		// 
		// Collision detection with ground plane (no friction)
		if (globs[p][1] + newVel[1] * time <= groundP[1])
		{
			// Collision with plane, solve for Fnormal
			double *xMinusP = new double[3]{ globs[p][0], globs[p][1] - groundP[1], globs[p][2] };
			double xDotN = xMinusP[1] * groundN[1]; // N is the normal of ground plane, only contributes in the y direction
			// Define -ks(dxDotN * N), where N = groundN = {0.0, 1.0, 0.0}, ks=groundKs
			xDotN = -groundKs * (xDotN * groundN[1]);
			double* term1 = new double[3]{ 0.0, xDotN, 0.0 };

			double vDotN = newVel[1] * groundN[1]; // N only contributes in the y direction. kd(vDotN * N), kd=groundKd
			vDotN = groundKd * vDotN * groundN[1];
			double* term2 = new double[3]{ 0.0, vDotN, 0.0 };

			// Finally, solve Fnormal
			double Fn = term1[1] - term2[1]; // x and z terms come out to 0.0, thus a single double is sufficient.

			FTotal[1] = FTotal[1] + Fn;
			// Solve and update velocity 
			newVel[1] = time * FTotal[1] / masses[p];
		}

		globs[p][0] = globs[p][0] + newVel[0] * time;
		globs[p][1] = globs[p][1] + newVel[1] * time;
		globs[p][2] = globs[p][2] + newVel[2] * time;
		//animTcl::OutputMessage("before pushing glob, globs[p][1]: %f", globs[p][1]);
		// Append new position to list
		newGlobs.push_back(globs[p]);
		// Append new velocity to list
		newV.push_back(newVel);	
	}
	//animTcl::OutputMessage("NEWGlobs[0]: %f %f %f", newGlobs[0][0], newGlobs[0][1], newGlobs[0][2]);
	mySystem->setState(newGlobs);
	mySystem->setVel(newV);
	return 0;

}	// GlobuleSim::step

int GlobuleSim::command(int argc, myCONST_SPEC char** argv)
{
	if (argc < 1)
	{
		animTcl::OutputMessage("simulator %s: wrong number of params.", m_name.c_str());
		return TCL_ERROR;
	}
	else if (strcmp(argv[0], "constants") == 0)
	{
		// m, n, b1
		if (argc == 4)
		{
			m = std::stof(argv[1]);
			n = std::stof(argv[2]);
			b1 = std::stof(argv[3]);
			myTime = std::stof(argv[4]);

		}
		else
		{
			animTcl::OutputMessage("simulator %s: wrong number of params for command 'constants'.", m_name.c_str());
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[0], "material") == 0)
	{
		if (argc == 2)
		{
			strcpy(material, argv[1]);
		}
		else
		{
			animTcl::OutputMessage("simulator %s: wrong number of params for command 'material'.", m_name.c_str());
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[0], "ground") == 0)
	{
		if (argc == 3)
		{
			groundKs = std::stod(argv[1]);
			groundKd = std::stod(argv[2]);
		}
		else
		{
			animTcl::OutputMessage("simulator %s: wrong number of params for command 'ground'.", m_name.c_str());
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[0], "gravity") == 0)
	{
		if (argc == 2)
		{
			gravity = std::stod(argv[1]);
		}
		else
		{
			animTcl::OutputMessage("simulator %s: wrong number of params for command 'gravity'.", m_name.c_str());
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[0], "drag") == 0)
	{
		if (argc == 2)
		{
			drag = std::stod(argv[1]);
		}
		else
		{
			animTcl::OutputMessage("simulator %s: wrong number of params for command 'drag'.", m_name.c_str());
			return TCL_ERROR;
		}
	}

	glutPostRedisplay();
	return TCL_OK;
} // GlobuleSim::command
