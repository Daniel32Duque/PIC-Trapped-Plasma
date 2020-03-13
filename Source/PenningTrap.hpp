/*
Written by: Daniel Duque
Last modified on 13 Mar 2020

Declarations for the Electrode and PenningTrap classes
This file contains a corresponding source file.
*/

#ifndef PENNINGTRAP_HPP
#define PENNINGTRAP_HPP

#include<vector>
#include<Eigen/SparseLU>
#include<fstream>
#include<unordered_map>

class Plasma;

class Electrode
{
private:
	double length;
	double potential;

public:
	Electrode(double aLength, double aPotential);
	Electrode(const Electrode& copiable);
	~Electrode();
	double getLength() const;
	double getPotential() const;
	void setPotential(double finalPotential);
};


class PenningTrap
{
private:
	double trapRadius;
	std::vector<Electrode> electrodes;
	std::vector<double> gaps;
	int Nz; //Number of cells along z direction i.e. Nz + 1 grid points
	int Nr; //Same, but only Nr grid points, the last grid point is fixed at the electrodes (boundary)
	double hr, hz; //Distance between grid points
	double lengthTrap; //Total length of all electrodes and gaps
	std::vector<double> potentialEnergiesHistory; //The total potential energy of all particles in the trap
	std::vector<double> timesSaved; //Times at which the state of the plasmas was stored
	Eigen::SparseMatrix<double> coefficients{ Nz * Nr + Nr, Nz * Nr + Nr };//Sparse matrix formed from finite difference method
	void generateSparse();
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;//keep solver loaded to solve fast for varying potentials
	Eigen::VectorXd potentialsVector{ Nz * Nr + Nr };//Array that stores phi in row-major order starting from (z = 0, r = 0)
	Eigen::SparseVector<double> RHS{ Nz * Nr + Nr };//Right hand side where boundary conditions come from
	void updateRHS(); //Set the RHS according to appropriate boundary conditions
	void solveLaplace(); //Solves Laplace equation whenever the potential in an electrode changes
	std::vector<int> limitLeft;
	std::vector<int> limitRight;
	std::vector<std::reference_wrapper<Plasma>> plasmas; //Reference to the plasmas inside of the trap
	void addPlasma(Plasma&);
	std::unordered_map<int, double> eFields; //"Smart" way of obtaining electric field fast avoids repeating grid points
	double getEField(int r, int z); //Get E field at one of the grid points
	double getEField(int r, double z); //Get field anywhere, based on E field at neighbour grid points.
	double getTotalPhi(int r, int z) const; //Get the sum of the potential fields of the trap plus all the plasmas in a grid point
	double getTotalPhi(int r, double z) const; //Total sum of Phi anywhere
	
public:
	PenningTrap(double radius, const std::vector<Electrode>& theElectrodes, const std::vector<double>& theGaps, int NumCellsZ, int NumCellsR);
	~PenningTrap();
	friend class Plasma;
	void extractTrapPotential(std::string fileName) const;//Store the potential field in a document called filename
	void extractTrapLaplacian(std::string fileName) const; //Store in filename, this is expected to be all zeros 
	void extractPlasmasHistories(std::string pathAndPreName) const; //History of all plasmas, separate file for position, speed, etc.
	void extractTrapParameters(std::string filename) const;//Store the intrinsic and static parameters of the trap i.e. dimensions, number of grids, etc.
	void setPotential(int indexElectrode, double newPotential); //Change potential of a single electrode, it changes the field as well
	double getLength() const;
	double getRadius() const;
	void movePlasmas(double deltaT); //Move all plasmas in the trap by a time deltaT
	void saveStates(double aTime); //Store the states of all plasmas
	void reserve(int desired); //Reserve how many times are you going to save the states of all plasmas. Use it if you will store a lot of times
};
#endif
