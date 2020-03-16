/*
Written by: Daniel Duque
Last modified on 16 Mar 2020

Definitions for the Electrode and PenningTrap classes
*/

#include"PenningTrap.hpp"
#include"Plasma.hpp"

Electrode::Electrode(double aLength, double aPotential) : length(aLength), potential(aPotential)
{
}
Electrode::Electrode(const Electrode& copiable) : length(copiable.length), potential(copiable.potential)
{
}
Electrode::~Electrode()
{}
double Electrode::getLength() const
{
	return length;
}
double Electrode::getPotential() const
{
	return potential;
}
void Electrode::setPotential(double aPotential)
{
	potential = aPotential;
}

/*----------------------------------------------------------------------------------------------------------
Penning Trap class
----------------------------------------------------------------------------------------------------------*/

PenningTrap::PenningTrap(double radius, const std::vector<Electrode>& theElectrodes, const std::vector<double>& theGaps, int NumCellsZ, int NumCellsR)
	: trapRadius(radius), electrodes(theElectrodes), gaps(theGaps), Nz(NumCellsZ), Nr(NumCellsR)
{
	if (electrodes.size() != gaps.size() + 1)
	{
		throw std::logic_error("Error number of gaps and electrodes; No. electrods should match No. gaps + 1");
	}
	for (unsigned int i = 0; i < electrodes.size(); ++i)
	{
		lengthTrap += electrodes[i].getLength();
		if (i < gaps.size())
		{
			lengthTrap += gaps[i];
		}
	}
	hz = lengthTrap / Nz;
	hr = trapRadius / Nr;
	coefficients.reserve(Eigen::VectorXi::Constant(Nz * Nr + Nr, 5)); //There are maximum 5 non-zero elements per column, reserve to optimize speed
	RHS.reserve(Nz + 1);
	generateSparse();
	solver.analyzePattern(coefficients);
	solver.factorize(coefficients);
	solveLaplace();
	//Now find the limit left and right indices for the "hole" where the plasma is going to be trapped initially
	//This is needed for the initial equilibrium solver
	//Because if the ends of the trap are grounded, then the trap thinks that these are "holes" that need to be filled with charge
	//This tells you the index of the closest grid point to the left of the centre of the trap
	int indexZ{ (int)floor(lengthTrap / (2 * hz)) };
	//Now starte moving to the left, and to the right, until the potential stops increasing/decreasing, this would tell you where the limits of your "hole" is
	//Move to the right
	for (int indexR = 0; indexR < Nr; ++indexR)
	{
		int limit = indexZ + 1;
		double change{ potentialsVector.coeffRef((Nz + 1) * indexR + limit + 1) - potentialsVector.coeffRef((Nz + 1) * indexR + limit) };
		double changeTwo;
		do
		{
			++limit;
			changeTwo = potentialsVector.coeffRef((Nz + 1) * indexR + limit + 1) - potentialsVector.coeffRef((Nz + 1) * indexR + limit);
		} while (change * changeTwo > 0 && limit + 1 < Nz);
		limitRight.push_back(limit);
	}
	//Similar for the left
	for (int indexR = 0; indexR < Nr; ++indexR)
	{
		int limit = indexZ;
		double change{ potentialsVector.coeffRef((Nz + 1) * indexR + limit - 1) - potentialsVector.coeffRef((Nz + 1) * indexR + limit) };
		double changeTwo;
		do
		{
			--limit;
			changeTwo = potentialsVector.coeffRef((Nz + 1) * indexR + limit - 1) - potentialsVector.coeffRef((Nz + 1) * indexR + limit);
		} while (change * changeTwo > 0 && limit - 1 > 0);
		limitLeft.push_back(limit);
	}
}
PenningTrap::~PenningTrap()
{}
void PenningTrap::generateSparse()
{
	//Look at the PDF for an explanation about where the coefficients come from
	double hz2{ std::pow(hz, -2) };
	double hr2{ std::pow(hr, -2) };
	//Need to populate all diagonals by applying the equation in the documentation to each grid point
	//For the 0th grid point
	coefficients.insert(0,0) = -2 * (hr2 + hz2);
	coefficients.insert(0,1) = 2 * hz2;
	coefficients.insert(0,Nz + 1) = 2 * hr2;
	//All grid points at r=0 except last
	for (int i = 1; i < Nz; ++i)
	{
		coefficients.insert(i, i - 1) = hz2;
		coefficients.insert(i,i) = -2 * (hr2 + hz2);
		coefficients.insert(i,i + 1) = hz2;
		coefficients.insert(i,i + Nz + 1) = 2 * hr2;
	}
	//Last grid point at r = 0;
	coefficients.insert(Nz, Nz - 1) = 2 * hz2;
	coefficients.insert(Nz,Nz) = -2 * (hr2 + hz2);
	//coefficients.insert(Nz,Nz + 1) = 0; //Sparse, no need.
	coefficients.insert(Nz,Nz * 2 + 1) = 2 * hr2;
	//All grid points except the row with Dirchlet boundary
	int lastRow{ Nz*Nr + Nr - Nz - 1 };//index of first element in last row
	for (int i = Nz + 1; i < lastRow; ++i)
	{
		double radius = floor(i / (Nz + 1)) * hr;
		coefficients.insert(i, i - Nz - 1) = hr2 - std::pow(2 * radius * hr, -1);
		if (i % (Nz + 1) == 0) //Left boundary
		{
			//coefficients.insert(i, i - 1) = 0;//Sparse, no need.
			coefficients.insert(i,i + 1) = 2 * hz2;
		}
		else if ((i + 1) % (Nz + 1) == 0) //Right boundary
		{
			coefficients.insert(i, i - 1) = 2 * hz2;
			//coefficients.insert(i, i + 1) = 0; //Sparse, no need.
		}
		else
		{
			coefficients.insert(i, i - 1) = hz2;
			coefficients.insert(i, i + 1) = hz2;
		}
		coefficients.insert(i, i) = -2 * (hr2 + hz2);
		coefficients.insert(i, i + Nz + 1) = hr2 + std::pow(2 * radius * hr, -1);
	}
	//Last row with Dirichlet boundary conditions
	double radius{ trapRadius - hr };
	//First element
	coefficients.insert(lastRow,lastRow - Nz - 1) = hr2 - std::pow(2 * radius * hr, -1);
	//coefficients.insert(lastRow, lastRow - 1) = 0;//Sparse, no need
	coefficients.insert(lastRow, lastRow) = -2 * (hr2 + hz2);
	coefficients.insert(lastRow, lastRow + 1) = 2 * hz2;
	//All except last in the row
	int lastElement{ lastRow + Nz };
	for (int i = lastRow + 1; i < lastElement; ++i)
	{
		coefficients.insert(i,i - Nz - 1) = hr2 - std::pow(2 * radius * hr, -1);
		coefficients.insert(i, i - 1) = hz2;
		coefficients.insert(i,i) = -2 * (hr2 + hz2);
		coefficients.insert(i, i + 1) = hz2;
	}
	//Last element
	coefficients.insert(lastElement, lastElement - Nz - 1) = hr2 - std::pow(2 * radius * hr, -1);
	coefficients.insert(lastElement, lastElement - 1) = 2 * hz2;
	coefficients.insert(lastElement,lastElement) = -2 * (hr2 + hz2);
	coefficients.makeCompressed(); //Clear buffer, removed empty spaces and redundant innerNNZs array (see Eigen Sparse documentation)
}
void PenningTrap::updateRHS()
{
	double hr2{ std::pow(hr, -2) };
	double radius{ trapRadius - hr };
	double matrixFactor{ hr2 + std::pow(2 * radius * hr, -1) };
	int N{ Nz * Nr + Nr - Nz - 1 };//First index in the last row
	int point{ 0 };
	double totalLength{ 0 };
	double boundary;
	for (unsigned int i = 0; i < electrodes.size(); ++i)
	{
		boundary = electrodes[i].getPotential();
		while (point * hz <= electrodes[i].getLength() + totalLength)
		{
			RHS.coeffRef(point + N) = -1 * matrixFactor * boundary;
			point++;
		}
		if (i < gaps.size())
		{
			while (point * hz < electrodes[i].getLength() + gaps[i] + totalLength)//Linear interpolation inbetween electrodes
			{
				boundary = (point * hz - electrodes[i].getLength() - totalLength) * (electrodes[i + 1].getPotential() - electrodes[i].getPotential()) / gaps[i] + electrodes[i].getPotential();
				RHS.coeffRef(point + N) = -1 * matrixFactor * boundary;
				point++;
			}
			totalLength += electrodes[i].getLength() + gaps[i];
		}
	}
	//Sometimes, if the number of grid cells is large, rounding errors cause the last grid point to be missed
	//fill the last points missed due to rounding errors with the appropriate end electrode
	while (point < Nz + 1)
	{
		RHS.coeffRef(point + N) = -1 * matrixFactor * electrodes.back().getPotential();
		point++;
	}
}
void PenningTrap::solveLaplace()
{
	updateRHS();
	potentialsVector = solver.solve(RHS);
}
void PenningTrap::addPlasma(Plasma& aPlasma)
{
	plasmas.push_back(aPlasma);
}
double PenningTrap::getEField(int indexR, int indexZ)
{
	int index{ (Nz + 1) * indexR + indexZ };//index in 1D array
	auto it = eFields.find(index);
	if (it != eFields.end())//If it found the element
	{
		return it->second;//seconds means the value
	}
	//else
	int indexLeft, indexRight;
	if (indexZ == 0 || indexZ == Nz)
	{
		eFields[index] = 0;
		return 0;
	}
	//else
	indexLeft = index - 1;
	indexRight = index + 1;
	double valueLeft{ potentialsVector.coeff(indexLeft) };
	double valueRight{ potentialsVector.coeff(indexRight) };
	for (const Plasma& aPlasma : plasmas)
	{
		valueLeft += aPlasma.selfPotential.coeff(indexLeft);
		valueRight += aPlasma.selfPotential.coeff(indexRight);
	}
	double fieldValue{ (valueLeft - valueRight) / (2 * hz) };
	eFields[index] = fieldValue;
	return fieldValue;
}
void PenningTrap::extractTrapPotential(std::string fileName) const
{
	Eigen::IOFormat fastFullPrecision(Eigen::FullPrecision, Eigen::DontAlignCols, "", "\n", "", "", "", "");
	std::ofstream newFile;
	newFile.open(fileName);
	newFile << potentialsVector.format(fastFullPrecision);
	newFile.close();
}
void PenningTrap::extractTrapLaplacian(std::string fileName) const
{
	Eigen::IOFormat fastFullPrecision(Eigen::FullPrecision, Eigen::DontAlignCols, "", "\n", "", "", "", "");
	std::ofstream newFile;
	newFile.open(fileName);
	newFile << Eigen::VectorXd { (coefficients * potentialsVector) - RHS }.format(fastFullPrecision);
	newFile.close();
}
void PenningTrap::extractPlasmasHistories(std::string pathAndPreName) const
{
	std::ofstream newFile;
	newFile.open(pathAndPreName + "Times.csv");
	newFile << std::setprecision(std::numeric_limits<double>::digits10);
	for (unsigned int i = 0; i < timesSaved.size(); ++i)
	{
		newFile << timesSaved[i];
		if (i < timesSaved.size() - 1)
		{
			newFile << ",";
		}
	}
	newFile.close();
	//Now extract the total potential energy
	newFile.open(pathAndPreName + "PotentialEnergies.csv");
	newFile << std::setprecision(std::numeric_limits<double>::digits10);
	for (unsigned int i = 0; i < potentialEnergiesHistory.size(); ++i)
	{
		newFile << potentialEnergiesHistory[i];
		if (i < potentialEnergiesHistory.size() - 1)
		{
			newFile << ",";
		}
	}
	newFile.close();
	//Plasma phase space history
	int i{ 1 };
	for (const Plasma& aPlasma : plasmas)
	{
		aPlasma.extractHistory(pathAndPreName);
		++i;
	}
}
void PenningTrap::extractTrapParameters(std::string filename) const
{
	std::ofstream newFile;
	newFile.open(filename);
	newFile << std::setprecision(std::numeric_limits<double>::digits10);
	newFile << trapRadius << '\n';
	for (unsigned int i = 0; i < electrodes.size(); ++i)
	{
		newFile << electrodes[i].getLength();
		i < electrodes.size() - 1 ? newFile << ',' : newFile << '\n';
	}
	for (unsigned int i = 0; i < electrodes.size(); ++i)
	{
		newFile << electrodes[i].getPotential();
		i < electrodes.size() - 1 ? newFile << ',' : newFile << '\n';
	}
	for (unsigned int i = 0; i < gaps.size(); ++i)
	{
		newFile << gaps[i];
		i < gaps.size() - 1 ? newFile << ',' : newFile << '\n';
	}
	newFile << Nz << ',' << Nr << '\n';
	newFile << hz << ',' << hr << '\n';
	newFile << lengthTrap;
	newFile.close();
}
void PenningTrap::setPotential(int indexElectrode, double newPotential)
{
	electrodes[indexElectrode].setPotential(newPotential);
	solveLaplace();
}
double PenningTrap::getLength() const
{
	return lengthTrap;
}
double PenningTrap::getRadius() const
{
	return trapRadius;
}
double PenningTrap::getEField(int r, double z)//There is a problem if z == lengthTrap, there is no point to the right. Either make sure this never happens, or test every time (which could slow down the code)
{
	int indexZ{ (int)floor(z / hz) };
	double fieldLeft{ getEField(r, indexZ) };
	double fieldRight{ getEField(r, indexZ + 1) };
	double dz{ z - indexZ * hz };
	double weightFactor{ dz / hz };
	return ( (1 - weightFactor) * fieldLeft + weightFactor * fieldRight );
}
double PenningTrap::getTotalPhi(int r, int z) const
{
	int index{ (Nz + 1) * r + z };
	double phi{ potentialsVector.coeff(index) };
	for (const Plasma& aPlasma : plasmas)
	{
		phi += aPlasma.selfPotential.coeff(index);
	}
	return phi;
}
double PenningTrap::getTotalPhi(int r, double z) const
{
	int indexZ{ (int)floor(z / hz) };
	double dz{ z - indexZ * hz };
	double weightFactor{ dz / hz };
	return ((1 - weightFactor) * getTotalPhi(r, indexZ) + weightFactor * getTotalPhi(r, indexZ + 1));
}
void PenningTrap::movePlasmas(double deltaT)
{
	for (Plasma& aPlasma : plasmas)
	{
		aPlasma.moveRings(deltaT);
	}
	for (Plasma& aPlasma : plasmas)
	{
		aPlasma.solvePoisson();
	}
	eFields.clear();
}
void PenningTrap::saveStates(double aTime)
{
	timesSaved.push_back(aTime);
	double potentialEnergy{ 0 };
	for (Plasma& aPlasma : plasmas)
	{
		aPlasma.saveState();
		potentialEnergy += aPlasma.getPotentialEnergy();
	}
	potentialEnergiesHistory.push_back(potentialEnergy);
}
void PenningTrap::saveStates(double aTime, int indexR)
{
	timesSaved.push_back(aTime);
	double potentialEnergy{ 0 };
	for (Plasma& aPlasma : plasmas)
	{
		aPlasma.saveState(indexR);
		potentialEnergy += aPlasma.getPotentialEnergy();
	}
	potentialEnergiesHistory.push_back(potentialEnergy);
}
void PenningTrap::reserve(int desired)
{
	timesSaved.reserve(desired);
	for (Plasma& aPlasma : plasmas)
	{
		aPlasma.reserve(desired);
	}
	potentialEnergiesHistory.reserve(desired);
}
