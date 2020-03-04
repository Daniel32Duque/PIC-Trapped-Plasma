/*
Written by: Daniel Duque
Last modified on 04 Mar 2020

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
void Plasma::extractPlasmaParameters(std::string fileName) const
{
	std::ofstream newFile;
	newFile.open(fileName);
	newFile << mass << '\n';
	newFile << charge << '\n';
	newFile << massMacro << '\n';
	newFile << chargeMacro;
	newFile.close();
}
void Plasma::extractInitialDensity(std::string fileName) const
{
	Eigen::IOFormat fastFullPrecision(Eigen::FullPrecision, Eigen::DontAlignCols, "", "\n", "", "", "", "");
	std::ofstream newFile;
	newFile.open(fileName);
	newFile << initialDensity.format(fastFullPrecision);
	newFile.close();
}
void Plasma::extractHistory(std::string preName) const
{
	std::ofstream newPositions, newSpeeds;
	newPositions.open(preName + "Positions" + name + ".csv");
	newSpeeds.open(preName + "Speeds" + name + ".csv");
	for (const MacroRing& aRing : rings)
	{
		aRing.printPositions(newPositions);
		aRing.printSpeeds(newSpeeds);
	}
	newPositions.close();
	newSpeeds.close();
}
double Plasma::getPotentialEnergy() const
{
	double potentialEnergy{ 0 };
	for (const MacroRing& aRing : rings)
	{
		int indexZ{ (int)floor(aRing.getZ() / refTrap.hz) };
		double fieldLeft{ refTrap.getTotalPhi(aRing.getR(), indexZ) };
		double fieldRight{ refTrap.getTotalPhi(aRing.getR(), indexZ + 1) };
		double dz{ aRing.getZ() - indexZ * refTrap.hz };
		double weightFactor{ dz / refTrap.hz };
		double phiHere{ (1 - weightFactor) * fieldLeft + weightFactor * fieldRight };
		potentialEnergy += phiHere * chargeMacro;
	}
	return potentialEnergy / 2;
}
void Plasma::estimateDensityProportions()
{
	int pointsZ = refTrap.Nz + 1;
	int pointsR = refTrap.Nr;
	//For each grid point, estimate density from Total Potential
	for (int indexR = 0; indexR < pointsR; ++indexR)
	{
		double phiCentralR{ refTrap.getTotalPhi(indexR, refTrap.lengthTrap / 2) };
		for (int indexZ = 0; indexZ < pointsZ; ++indexZ)
		{
			initialDensity.coeffRef(pointsZ * indexR + indexZ) = exp(-(charge / (KB * temperature)) * (refTrap.getTotalPhi(indexR, indexZ) - phiCentralR));
		}
	}
}
void Plasma::fitDensityProportionToProfile(double shape, double scale)//profile = exp(-(r/scale)^shape), No normalization required
{
	int pointsZ = refTrap.Nz + 1;
	int pointsR = refTrap.Nr;
	double hr{ refTrap.hr };
	double hz{ refTrap.hz };
	//For each grid point, estimate density from Total Potential
	for (int indexR = 0; indexR < pointsR; ++indexR)
	{
		//Convert from volume density into flat space density
		double profileAtR{ 0 };
		for (int indexZ = 0; indexZ < pointsZ; ++indexZ)
		{
			profileAtR += initialDensity.coeffRef(pointsZ * indexR + indexZ) * hz;
		}
		double neededFactor{ exp(-pow((indexR * hr * 1000 / scale), shape)) / profileAtR }; //This 1000 is because the profile is given in milimeters
		for (int indexZ = 0; indexZ < pointsZ; ++indexZ)
		{
			initialDensity.coeffRef(pointsZ * indexR + indexZ) *= neededFactor;
		}
	}
}
void Plasma::normalizeDensityToTotalCharge()
{
	int pointsZ = refTrap.Nz + 1;
	int pointsR = refTrap.Nr;
	double hr{ refTrap.hr };
	double hz{ refTrap.hz };
	//First calculate what the total charge is right now
	double currentCharge{ 0 };
	for (int indexR = 0; indexR < pointsR; ++indexR)
	{
		double volume; //Volume of the MacroRing
		if (indexR == 0)
		{
			volume = PI * hz * hr * hr / 4;
		}
		else
		{
			volume = hz * hr * 2 * PI * indexR * hr;
		}
		for (int indexZ = 0; indexZ < pointsZ; ++indexZ)
		{
			currentCharge += initialDensity.coeffRef(pointsZ * indexR + indexZ) * volume;
		}
	}
	//Now multiply every density by the appropriate correction to get the expected total charge
	double correction{ totalCharge / currentCharge };
	for (int indexR = 0; indexR < pointsR; ++indexR)
	{
		for (int indexZ = 0; indexZ < pointsZ; ++indexZ)
		{
			initialDensity.coeffRef(pointsZ * indexR + indexZ) *= correction;
		}
	}
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
void Plasma::loadProfile(double aTemperature, double aTotalCharge, double shape, double scale, int numMacro, double KSThreshold)
{
	//Check input makes sense
	if (aTemperature <= 0 || shape <= 0 || scale <= 0)
	{
		throw std::logic_error("Temperature, shape, and scale all need to be positive");
	}
	if (aTotalCharge * charge < 0)
	{
		throw std::logic_error("Total charge and the plasma type charge must have the same sign");
	}
	if (KSThreshold <= 0 || KSThreshold >= 1)
	{
		throw std::logic_error("The Kolmogorov Smirnov distance threshold needs to be a number between 0 and 1");
	}
	temperature = aTemperature;
	totalCharge = aTotalCharge;
	chargeMacro = totalCharge / numMacro;
	massMacro = mass * chargeMacro / charge;
	rings.clear();
	rings.reserve(numMacro);
	//Iteratively solve the equilibrium density
	initialDensity.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
	//This is our first guess of the charge distribution
	selfPotential = refTrap.solver.solve(initialDensity);
	estimateDensityProportions();
	fitDensityProportionToProfile(shape, scale);
	normalizeDensityToTotalCharge();
	//First estimate what the Pseudo Kolmogorov Smirnov distance is between the first and second guess
	//This will give you an indication of how low or high the temperature is i.e. how much of the new guess would it be good to add into the previous guess
	//At low temperatures you need to add less of the new estimate because the peaks are very narrow and everything starts to blow
	//This Pseudo KS distance I made up is twice the normal KS distance, this is because I assume symmetry wrt centre of the trap, and I want to compare the difference in half the trap
	//Solve Poisson's equation for the charge density
	//Then estimate a thermal equilibrium density from that potential (normalized to totalCharge).
	//Tune the obtained density according to the expected profile
	//Combine this new estimate with the previous one. Linear combination with coefficients depend on the Kolmogorov Smirnov distance between the first and second guess
	//Ugly for loop just to have the i let us know it is the first time inside
	double KSCorrection{ 0 };
	for (int i = 0; ; )
	{
		Eigen::VectorXd previousDensity(initialDensity);
		//Format the density as the right hand side of Poissons equation
		for (int j = 0; j < refTrap.Nz * refTrap.Nr + refTrap.Nr; ++j)
		{
			initialDensity.coeffRef(j) = -initialDensity.coeffRef(j) / epsilon;
		}
		//Obtain the next density guess
		selfPotential = refTrap.solver.solve(initialDensity);
		estimateDensityProportions();
		fitDensityProportionToProfile(shape, scale);
		normalizeDensityToTotalCharge();
		//Compute the Kolmogorov Smirnov distance between both estimates
		std::vector<double> oldDistribution, newDistribution; //Store here the 1D Distribution along z
		double newInput{ 0 };
		double oldInput{ 0 };
		for (int indexZ = 0; indexZ < refTrap.Nz + 1; ++indexZ)
		{
			double volume;
			for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
			{
				if (indexR == 0)
				{
					volume = PI * refTrap.hz *  refTrap.hr *  refTrap.hr / 4;
				}
				else
				{
					volume = refTrap.hz *  refTrap.hr * 2 * PI * indexR *  refTrap.hr;
				}
				newInput += volume * initialDensity.coeffRef((refTrap.Nz + 1) * indexR + indexZ);
				oldInput += volume * previousDensity.coeffRef((refTrap.Nz + 1) * indexR + indexZ);
			}
			newDistribution.push_back(newInput / totalCharge);
			oldDistribution.push_back(oldInput / totalCharge);
		}
		double maxDistance{ 0 };
		for (unsigned int h = 0; h < newDistribution.size(); ++h)
		{
			double test{ abs(newDistribution[h] - oldDistribution[h]) };
			if (test > maxDistance)
			{
				maxDistance = test;
			}
		}
		//If this is the first time running, set this as the correction factor.
		if (i == 0)
		{
			++i;
			//Pseudo KS distance, is just twice the amount because we assume symmetry along z wrt the centre of the trap so we only care about half the distribution
			KSCorrection = 2 * maxDistance;
		}
		if (maxDistance < KSThreshold)
		{
			break;
		}
		//Now combine both the new and old density based on the KSCorrection
		//This linear combinations retains the total integral and the profile
		initialDensity = KSCorrection * previousDensity + (1 - KSCorrection) * initialDensity;
	}
	//Ok, now I have a density grid. I need now to populate macro-particles to match this grid density.
}
