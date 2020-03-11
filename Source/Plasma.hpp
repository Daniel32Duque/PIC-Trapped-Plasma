/*
Written by: Daniel Duque
Last modified on 11 Mar 2020

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
#include<iostream>
#include<random>

class MacroRing
{
private:
	int posR;//Fixed at one of the grid divisions along r. This indicates the index of the radial grid point [0,Nr]
	double posZ; //Position along z. Position in meters along the trap [0, trapLength]
	double speed; //speed along z, there is no r transfer.
	std::vector<double> historySpeed;
	std::vector<double> historyZ;
	int getR() const;
	double getZ() const;
	double getSpeed() const;
	void setZ(double newZ);
	void setSpeed(double newSpeed);
	void saveState(); //Stores the position and speed.
	void reserve(int desired); //Saves space for how many times you will save state.
	void printPositions(std::ofstream& file) const;
	void printSpeeds(std::ofstream& file) const;
	MacroRing(int r, double z, double aSpeed);
	friend Plasma;
};

class Plasma
{
private:
	PenningTrap& refTrap;
	const std::string name;//A name to identify the output files e.g. Electrons, Antiprotons, etc.
	const double mass; //Mass of type of particles that form the plasma. It is NOT the total mass of the plasma. (In kg)
	const double charge; //Not the net charge of the plasma. It is the charge of the particle that forms the plasma. (In Coulombs)
	double chargeMacro; //Charge of macro ring at r=0. The charge of a general macro ring at a given r is given by (indexR == 0 ? macroCharge : 8 * indexR * macroCharge)
	double macroChargeDensity; //Charge increases in the same proportion as volume along r. Each macroparticle has a fixed charge density in Coulombs per meter cubed
	std::vector<MacroRing> rings; //All macro-Rings of a given plasma are stored here
	double temperature; //In Kelvin
	Eigen::VectorXd selfPotential{ refTrap.Nz * refTrap.Nr + refTrap.Nr };
	Eigen::VectorXd RHS{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// = - charge density / permittivity of free space
	void updateRHS();
	void solvePoisson();
	void moveRings(double deltaT); //Move all particles by deltaT seconds
	void saveState(); //Store state of each of the rings
	void reserve(int desired);
	void extractHistory(std::string preName) const;
	double getPotentialEnergy() const;
	Eigen::VectorXd initialDensity{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// Initial density with which the plasma was loaded
	void estimateDensityProportions(); //Estimate the densities w.r.t. central densities in local thermal equilibrium from total potential. i.e. a value of 1 for all radius at the centre. Need to determine the actual values based on the profile
	void fitDensityProportionToProfile(double shape, double scale); //Weight each radius to match the desired profile
	void normalizeDensityToTotalCharge(double totalCharge); //Weight each radius so that the total charge matches the desired total charge
public:
	Plasma(PenningTrap& trap, std::string name, double mass, double charge);
	~Plasma();
	friend PenningTrap;
	void extractSelfPotential(std::string fileName) const;//Store the potential field in a document called filename
	void extractPlasmaParameters(std::string filename) const;//Extract the plasma parameters
	void extractInitialDensity(std::string filename) const;//This is NOT the actual initial density obtained from macro-Particles. This is the EXPECTED initial density
	int getNumMacro() const;//Return the current number of macro-particles still in the trap
	//Loading routines
	/*----------------------------------------------------------------------------------
	These routines below are the only thing that you should change/add to the code.
	Add a different loading routine depending on how you want the plasma to be loaded
	into the trap.
	The template for a loading routine should be:
	1) Check that input makes sense.
	2) Assign charge and charge density of a MacroRing (same for all).
	3) Clear the rings vector (This is like cleaning the trap before loading a plasma).
	4) Push the MacroRings into the rings vector according to the loading you want.
	5) Solve Poisson's Equation. Just call the method; literally one line: solvePoisson();

	The code only requires the loading to populate the rings vector i.e. telling the code how
	and where are the particles initially.
	-------------------------------------------------------------------------------------*/
	void loadOneDUniform(int numMacro, double totalCharge, double lengthLine, int r); //equally spaced particles at fixed r in a length centred in the center of the trap.
	void loadSingleRing(double totalCharge, int r, double Z, double speed); //Single ring at a position with a velocity
	void loadProfile(double aTemperature, double totalCharge, double shape, double scale, int numMacro, double KSThreshold); //Load using a density profile (as obtained from MCP), shape and scale define the generalized normal fit. Assume local thermal equilibrium along each radius
	void loadDensityFile(std::string fileName, double aTemperature, int numMacro);

};

#endif
