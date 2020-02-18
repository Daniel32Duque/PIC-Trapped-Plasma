/*
Written by: Daniel Duque
Last modified on 10 Dec 2019

Declarations for the Plasma class
This file contains a corresponding source file.
*/

#ifndef PLASMA_HPP
#define PLASMA_HPP

#include"Constants.hpp"
#include"PenningTrap.hpp"
#include<Eigen/SparseLU>
#include<array>
#include<vector>

class MacroRing
{
private:
	int posR;//Fixed at one of the grid divisions along r
	double posZ; //Position along z
	double speed; //speed along z, there is no r transfer.
	std::vector<double> historySpeed;
	std::vector<double> historyZ;
public:
	MacroRing(int r, double z, double aSpeed);
	~MacroRing();
	int getR() const;
	double getZ() const;
	double getSpeed() const;
	void setZ(double newZ);
	void setSpeed(double newSpeed);
	void saveState(); //Stores the position and speed.
	void reserve(int desired); //Saves space for how many times you will save state.
	void printPositions(std::ofstream& file) const;
	void printSpeeds(std::ofstream& file) const;
};

class Plasma
{
private:
	PenningTrap& refTrap;
	const std::string name;//A name to identify the output files e.g. Electrons, Antiprotons, etc.
	const double mass; //Mass of type of particles that form the plasma. It is NOT the total mass of the plasma.
	const double charge; //Not the net charge of the plasma. It is the charge of the particle that forms the plasma.
	double chargeMacro; //Charge of each MacroRing, all rings have the same charge regardless of r position.
	double massMacro; //Mass of each MacroRing.
	std::vector<MacroRing> rings;
	//double temperature;
	Eigen::VectorXd selfPotential{ refTrap.Nz * refTrap.Nr + refTrap.Nr };
	Eigen::VectorXd RHS{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// = - charge density / permittivity of free space
	void updateRHS();
	void solvePoisson();
	void moveRings(double deltaT); //Move all particles by deltaT seconds
	void saveState(); //Store state of each of the rings
	void reserve(int desired);
	void extractHistory(std::string preName) const;
	double getPotentialEnergy() const;
public:
	Plasma(PenningTrap& trap, std::string name, double mass, double charge);
	~Plasma();
	friend PenningTrap;
	void extractSelfPotential(std::string fileName) const;//Store the potential field in a document called filename
	void extractPlasmaParameters(std::string filename) const;//Extract the plasma parameters
	//Loading routines
	/*----------------------------------------------------------------------------------
	These routines below are the only thing that you should change/add to the code.
	Add a different loading routine depending on how you want the plasma to be loaded
	into the trap.
	The template for a loading routine should be:
	1) Check that input makes sense.
	2) Assign charge and mass of a MacroRing (same for all).
	3) Clear the rings vector (This is like cleaning the trap before loading a plasma).
	4) Push the MacroRings into the rings vector according to the loading you want.
	5) Solve Poisson's Equation. Just call the method; literally one line: solvePoisson();

	The code only requires the loading to populate the rings vector i.e. telling the code how
	and where are the particles initially.
	-------------------------------------------------------------------------------------*/
	void loadOneDUniform(int numMacro, double chargeMacro, double lengthLine, int r); //equally spaced particles at fixed r in a length centred in the center of the trap.
	void loadSingleRing(double chargeMacro, int r, double Z, double speed); //Single ring at a position with a velocity
};

#endif
