/*
Written by: Daniel Duque
Last modified on 10 Dec 2019

Definitions for the Plasma class
This file contains a corresponding source file.
*/
#include"Plasma.hpp"

MacroRing::MacroRing(int r, double z, double aSpeed) : posR(r), posZ(z), speed(aSpeed)
{
}
MacroRing::~MacroRing()
{
}
int MacroRing::getR() const
{
	return posR;
}
double MacroRing::getZ() const
{
	return posZ;
}
double MacroRing::getSpeed() const
{
	return speed;
}
void MacroRing::setZ(double newZ)
{
	posZ = newZ;
}
void MacroRing::setSpeed(double newSpeed)
{
	speed = newSpeed;
}
void MacroRing::saveState()
{
	historyZ.push_back(posZ);
	historySpeed.push_back(speed);
}
void MacroRing::reserve(int desired)
{
	historyZ.reserve(desired);
	historySpeed.reserve(desired);
}
void MacroRing::printPositions(std::ofstream& file) const
{
	if (!historyZ.empty())
	{
		file << posR;
		file << ",";
		for (unsigned int i = 0; i < historyZ.size() - 1; ++i)
		{
			file << historyZ[i] << ",";
		}
		file << historyZ.back() << "\n";
	}
}
void MacroRing::printSpeeds(std::ofstream& file) const
{
	for (unsigned int i = 0; i < historySpeed.size(); ++i)
	{
		file << historySpeed[i];
		i < historySpeed.size() - 1 ? file << "," : file << "\n";
	}
}

/*----------------------------------------------------------------------------------------------------------
Plasma class
----------------------------------------------------------------------------------------------------------*/

Plasma::Plasma(PenningTrap& trap, std::string aName, double aMass, double aCharge)
	: refTrap(trap) , name(aName), mass(aMass), charge(aCharge)
{
	refTrap.addPlasma(*this);
}
Plasma::~Plasma()
{
	//Should remove itself from the plasma, look later into this
}
void Plasma::updateRHS()
{
	//The RHS is: (minus) chargeDensity/epsilon
	//Performs a First-order weighting (Area weighting)
	RHS.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
	double hr{ refTrap.hr };
	double hz{ refTrap.hz };
	for (const MacroRing& aRing : rings)
	{
		int indexR{ aRing.getR() };
		int indexZ{ (int)floor(aRing.getZ() / hz) };
		int indexRHS{ (refTrap.Nz + 1) * indexR + indexZ };
		double z{ aRing.getZ() - indexZ * hz };//Distance from left grid point to particle
		double weightFactor{ z / hz };//Weight going to grid point on the right (number between 0 and 1)
		double volume; //Volume of the MacroRing
		if (indexR == 0)
		{
			volume = PI * hz * hr * hr / 4;
		}
		else
		{
			volume = hz * hr * 2 * PI * indexR * hr;
		}
		RHS.coeffRef(indexRHS) += -chargeMacro * (1 - weightFactor) / (volume * epsilon);//Point to the left
		RHS.coeffRef(indexRHS + 1) += -chargeMacro * weightFactor / (volume * epsilon);//Point to the right
	}
}
void Plasma::solvePoisson()
{
	updateRHS();
	selfPotential = refTrap.solver.solve(RHS);
}
void Plasma::moveRings(double deltaT)
{
	for (unsigned int i = 0; i < rings.size(); )
	{
		//Leapfrog method
		double forceOld{ chargeMacro * refTrap.getEField(rings[i].getR(), rings[i].getZ()) };
		double vNew{ deltaT * forceOld / massMacro + rings[i].getSpeed() };
		double zNew{ deltaT * vNew + rings[i].getZ() };
		//remove elements than escape the trap
		if (zNew < refTrap.getLength() && zNew > 0)
		{
			rings[i].setZ(zNew);
			rings[i].setSpeed(vNew);
			++i;
		}
		else
		{
			std::swap(rings[i], rings.back());//check that this std::swap works efficiently as expected
			rings.pop_back();
		}
	}
}
void Plasma::extractSelfPotential(std::string fileName) const
{
	Eigen::IOFormat fastFullPrecision(Eigen::FullPrecision, Eigen::DontAlignCols, "", "\n", "", "", "", "");
	std::ofstream newFile;
	newFile.open(fileName);
	newFile << selfPotential.format(fastFullPrecision);
	newFile.close();
}
void Plasma::extractHistory(std::string preName, std::string postName) const
{
	std::ofstream newPositions, newSpeeds;
	newPositions.open(preName + "Positions" + postName);
	newSpeeds.open(preName + "Speeds" + postName);
	for (const MacroRing& aRing : rings)
	{
		aRing.printPositions(newPositions);
		aRing.printSpeeds(newSpeeds);
	}
	newPositions.close();
	newSpeeds.close();
}
void Plasma::saveState()
{
	for (MacroRing& aRing : rings)
	{
		aRing.saveState();
	}
}
void Plasma::reserve(int desired)
{
	for (MacroRing& aRing : rings)
	{
		aRing.reserve(desired);
	}
}
void Plasma::loadOneDUniform(int numMacro, double aChargeMacro, double lengthLine, int r)
{
	//Creates equally spaced charges at r = 0 along a centred line of length lengthLine
	if (lengthLine >= refTrap.getLength() || r >= refTrap.Nr - 1)
	{
		throw std::logic_error("Length of charge must be less than the length of the trap and r inside the trap");
	}
	if (aChargeMacro * charge < 0)
	{
		throw std::logic_error("Charge of MacroRing and the plasma type must have the same sign");
	}
	chargeMacro = aChargeMacro;
	massMacro = mass * chargeMacro / charge;
	rings.clear();
	rings.reserve(numMacro);
	//Divide lengthLine in numMacro + 1 cells, which corresponds to numMacro + 2 points
	//Put the particles along the points except for the first and last point
	double start{ (refTrap.getLength() - lengthLine) / 2 };
	double hz{ lengthLine / (numMacro + 1) };
	for (int i = 1; i <= numMacro; ++i)
	{
		rings.push_back(MacroRing(r, start + i * hz, 0));
	}
	solvePoisson();
}
void Plasma::loadSingleRing(double aChargeMacro, int r, double Z, double speed)
{
	if (Z >= refTrap.getLength() || Z <= 0 || r >= refTrap.Nr - 1)
	{
		throw std::logic_error("Input r,z is not inside the trap");
	}
	if (aChargeMacro * charge < 0)
	{
		throw std::logic_error("Charge of MacroRing and the plasma type must have the same sign");
	}
	chargeMacro = aChargeMacro;
	massMacro = mass * chargeMacro / charge;
	rings.clear();
	rings.push_back(MacroRing(r, Z, speed));
	solvePoisson();
}
