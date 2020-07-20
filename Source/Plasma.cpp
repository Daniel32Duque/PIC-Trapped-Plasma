/*
Written by: Daniel Duque
Last modified on 20 Jul 2020

Definitions for the Plasma class
*/
#include"Plasma.hpp"

MacroRing::MacroRing(int r, double z, double aSpeed) : posR(r), posZ(z), speed(aSpeed)
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
		RHS.coeffRef(indexRHS) += -macroChargeDensity * (1 - weightFactor) / epsilon;//Point to the left
		RHS.coeffRef(indexRHS + 1) += -macroChargeDensity * weightFactor / epsilon;//Point to the right
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
		double vNew{ deltaT * refTrap.getEField(rings[i].getR(), rings[i].getZ()) * charge / mass + rings[i].getSpeed() };
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
	newFile << std::setprecision(std::numeric_limits<double>::digits10);
	newFile << mass << '\n';
	newFile << charge << '\n';
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
int Plasma::getNumMacro() const
{
	return rings.size();
}
int Plasma::getNumMacroCentralWell() const
{
	int numIn{ 0 };
	for (const MacroRing& aRing : rings)
	{
		if (aRing.posZ >= refTrap.limitLeft[aRing.posR] * refTrap.hz && aRing.posZ <= refTrap.limitRight[aRing.posR] * refTrap.hz)
		{
			++numIn;
		}
	}
	return numIn;
}
double Plasma::getAverageTemperature() const
{
	double T{ 0 };
	//Number of speeds stored minus one (need to average with next one, so ignore the last one)
	int times{ (int)rings[0].historySpeed.size() - 1};
	double numParticles{ 0 };
	for (const MacroRing& aRing : rings)
	{
		double aMass{ aRing.posR == 0 ? massMacro : 8 * aRing.posR * massMacro };
		numParticles += aMass / mass;
	}
	for (int i = 0; i < times; ++i)
	{
		double KE{ 0 };
		for (const MacroRing& aRing : rings)
		{
			double aMass{ aRing.posR == 0 ? massMacro : 8 * aRing.posR * massMacro };
			double speed{ (aRing.historySpeed[i] + aRing.historySpeed[i + 1]) / 2 };
			KE += 0.5 * aMass * speed * speed;
		}
		T += 2 * KE / (KB * numParticles);
	}
	return T / times;
}
double Plasma::getstdDeviation() const
{
	double mean{ getAverageTemperature() };
	double T{ 0 };
	//Number of speeds stored minus one (need to average with next one, so ignore the last one)
	int times{ (int)rings[0].historySpeed.size() - 1 };
	double numParticles{ 0 };
	for (const MacroRing& aRing : rings)
	{
		double aMass{ aRing.posR == 0 ? massMacro : 8 * aRing.posR * massMacro };
		numParticles += aMass / mass;
	}
	for (int i = 0; i < times; ++i)
	{
		double KE{ 0 };
		for (const MacroRing& aRing : rings)
		{
			double aMass{ aRing.posR == 0 ? massMacro : 8 * aRing.posR * massMacro };
			double speed{ (aRing.historySpeed[i] + aRing.historySpeed[i + 1]) / 2 };
			KE += 0.5 * aMass * speed * speed;
		}
		T += pow((2 * KE / (KB * numParticles)) - mean, 2);
	}
	return sqrt(T / (times - 1));
}
double Plasma::getTemperature() const
{
	double numParticles{ 0 };
	for (const MacroRing& aRing : rings)
	{
		double aMass{ aRing.posR == 0 ? massMacro : 8 * aRing.posR * massMacro };
		numParticles += aMass / mass;
	}
	double KE{ 0 };
	for (const MacroRing& aRing : rings)
	{
		double aMass{ aRing.posR == 0 ? massMacro : 8 * aRing.posR * massMacro };
		double speed{ (aRing.historySpeed.end()[-2] + aRing.historySpeed.back()) / 2 };
		KE += 0.5 * aMass * speed * speed;
	}
	return 2 * KE / (KB * numParticles);
}
void Plasma::extractHistory(std::string preName) const
{
	std::ofstream newPositions, newSpeeds;
	newPositions.open(preName + "Positions" + name + ".csv");
	newPositions << std::setprecision(std::numeric_limits<double>::digits10);
	newSpeeds.open(preName + "Speeds" + name + ".csv");
	newSpeeds << std::setprecision(std::numeric_limits<double>::digits10);
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
		potentialEnergy += refTrap.getTotalPhi(aRing.getR(), aRing.getZ()) * (aRing.getR() == 0 ? chargeMacro : aRing.getR() * 8 * chargeMacro);
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
			if (indexZ < refTrap.limitLeft[indexR] || indexZ > refTrap.limitRight[indexR])
			{
				initialDensity.coeffRef(pointsZ * indexR + indexZ) = 0;
			}
			else
			{
				initialDensity.coeffRef(pointsZ * indexR + indexZ) = exp(-(charge / (KB * temperature)) * (refTrap.getTotalPhi(indexR, indexZ) - phiCentralR));
			}		
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
void Plasma::normalizeDensityToTotalCharge(double totalCharge)
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
void Plasma::saveState(int indexR)
{
	for (MacroRing& aRing : rings)
	{
		if (aRing.getR() == indexR)
		{
			aRing.saveState();
		}
	}
}
void Plasma::reserve(int desired)
{
	for (MacroRing& aRing : rings)
	{
		aRing.reserve(desired);
	}
}
double Plasma::getCentralDensity() const
{
	int pointsZ = refTrap.Nz + 1;
	if (pointsZ - 1 % 2 == 0)
	{
		return initialDensity.coeffRef((pointsZ - 1) / 2) / charge;
	}
	else
	{
		return (initialDensity.coeffRef(pointsZ / 2) + initialDensity.coeffRef((pointsZ / 2) - 1)) / (2 * charge);
	}
}
void Plasma::loadProfile(double aTemperature, double aTotalCharge, double shape, double scale, int numMacro, double KSThreshold)
{
	//Check input makes sense
	if (aTemperature <= 0 || shape <= 0 || scale <= 0)
	{
		throw std::logic_error("Temperature, shape, and scale all need to be positive");
	}
	if (numMacro <= 0)
	{
		throw std::logic_error("Number of macro-particles has to be a positive integer");
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
	//Iteratively solve the equilibrium density
	initialDensity.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
	//This is our first guess of the charge distribution
	selfPotential = refTrap.solver.solve(initialDensity);
	estimateDensityProportions();
	fitDensityProportionToProfile(shape, scale);
	normalizeDensityToTotalCharge(aTotalCharge);
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
		normalizeDensityToTotalCharge(aTotalCharge);
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
			newDistribution.push_back(newInput / aTotalCharge);
			oldDistribution.push_back(oldInput / aTotalCharge);
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
		std::cout << maxDistance << '\n';
		if (maxDistance < KSThreshold)
		{
			break;
		}
		//Now combine both the new and old density based on the KSCorrection
		//This linear combinations retains the total integral and the profile
		initialDensity = KSCorrection * previousDensity + (1 - KSCorrection) * initialDensity;
	}
	//Ok, now I have a density grid. I need now to populate macro-particles to match this grid density.
	//How many particles am I placing at each radial index
	std::vector<std::vector<double>> cumulativeAtR;
	for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
	{
		std::vector<double> cumulative;
		double volume;
		if (indexR == 0)
		{
			volume = PI * refTrap.hz *  refTrap.hr *  refTrap.hr / 4;
		}
		else
		{
			volume = refTrap.hz *  refTrap.hr * 2 * PI * indexR *  refTrap.hr;
		}
		double aChargeR{ 0 };
		for (int indexZ = 0; indexZ < refTrap.Nz + 1; ++indexZ)
		{
			aChargeR += volume * initialDensity.coeffRef((refTrap.Nz + 1) * indexR + indexZ);
			cumulative.push_back(aChargeR);
		}
		cumulativeAtR.push_back(cumulative);
	}
	double futureMacroCharge{ cumulativeAtR[0].back() };
	for (unsigned int i = 1; i < cumulativeAtR.size(); ++i)
	{
		futureMacroCharge += cumulativeAtR[i].back() / (8 * i);
	}
	chargeMacro = futureMacroCharge / numMacro;
	massMacro = chargeMacro * mass / charge;
	macroChargeDensity = 4 * chargeMacro / (PI * refTrap.hz * refTrap.hr * refTrap.hr);
	std::vector<int> numAtR;
	numAtR.push_back((int)round(cumulativeAtR[0].back() / chargeMacro));
	for (unsigned int i = 1; i < cumulativeAtR.size(); ++i)
	{
		numAtR.push_back((int)round(cumulativeAtR[i].back() / (8 * i * chargeMacro)));
	}
	//Now I know how many particles I want to place at each radius.
	//Distribute them using an inverse transform sampling
	//But without any randomness, just equally space the uniform distribution
	rings.clear();
	rings.reserve(numMacro);//this is not exactly the number of particles we will load, but it is close
	//Use std::random to produce boltzmann distributed speeds
	//this are simply gaussian with the appropriate std deviation as we are looking at 1 dimension only
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0, sqrt(KB * aTemperature / mass));
	for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
	{
		double deltaQ{ cumulativeAtR[indexR].back() / (numAtR[indexR] + 1) };
		int currentIndex{ 0 };
		for (int i = 0; i < numAtR[indexR]; ++i)
		{
			double currentInvert{ deltaQ * (i + 1) };
			while (abs(cumulativeAtR[indexR][currentIndex]) < abs(currentInvert))
			{
				currentIndex++;
			}
			//The point to invert is between currentIndex and currentIndex - 1
			//Assume linear interpolation between this two points
			double invertedPos{ (currentIndex - 1) * refTrap.hz + refTrap.hz / 2 + refTrap.hz * (currentInvert - cumulativeAtR[indexR][currentIndex - 1]) / (cumulativeAtR[indexR][currentIndex] - cumulativeAtR[indexR][currentIndex - 1]) };
			rings.push_back(MacroRing(indexR, invertedPos, distribution(generator)));
		}
	}
	std::cout << "Loading " << rings.size() << " macro-particles from which " << numAtR[0] << " are at r=0.\n";
	solvePoisson();
}
void Plasma::loadDensityFile(std::string fileName, double aTemperature, int numMacro)
{
	//The file should just be a single column, in the same way the output of initialDensity is given if you loadProfile for example
	//Check input makes sense
	if (aTemperature <= 0)
	{
		throw std::logic_error("Temperature needs to be positive");
	}
	if (numMacro <= 0)
	{
		throw std::logic_error("Number of macro-particles has to be a positive integer");
	}
	temperature = aTemperature;
	//Read the file into the initial density
	initialDensity.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
	//Assume the file is in the correct format.
	double readNumber;
	std::ifstream rFile(fileName);
	int currentIndex{ 0 };
	while (rFile >> readNumber)
	{
		initialDensity.coeffRef(currentIndex) = readNumber;
		++currentIndex;
	}
	if (currentIndex != refTrap.Nz * refTrap.Nr + refTrap.Nr)
	{
		throw std::logic_error("The number of grid points in the file do not match this trap.");
	}
	//Ok, now I have a density grid. I need now to populate macro-particles to match this grid density.
	//How many particles am I placing at each radial index
	std::vector<std::vector<double>> cumulativeAtR;
	for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
	{
		std::vector<double> cumulative;
		double volume;
		if (indexR == 0)
		{
			volume = PI * refTrap.hz *  refTrap.hr *  refTrap.hr / 4;
		}
		else
		{
			volume = refTrap.hz *  refTrap.hr * 2 * PI * indexR *  refTrap.hr;
		}
		double aChargeR{ 0 };
		for (int indexZ = 0; indexZ < refTrap.Nz + 1; ++indexZ)
		{
			aChargeR += volume * initialDensity.coeffRef((refTrap.Nz + 1) * indexR + indexZ);
			cumulative.push_back(aChargeR);
		}
		cumulativeAtR.push_back(cumulative);
	}
	double futureMacroCharge{ cumulativeAtR[0].back() };
	for (unsigned int i = 1; i < cumulativeAtR.size(); ++i)
	{
		futureMacroCharge += cumulativeAtR[i].back() / (8 * i);
	}
	chargeMacro = futureMacroCharge / numMacro;
	massMacro = chargeMacro * mass / charge;
	macroChargeDensity = 4 * chargeMacro / (PI * refTrap.hz * refTrap.hr * refTrap.hr);
	std::vector<int> numAtR;
	numAtR.push_back((int)round(cumulativeAtR[0].back() / chargeMacro));
	for (unsigned int i = 1; i < cumulativeAtR.size(); ++i)
	{
		numAtR.push_back((int)round(cumulativeAtR[i].back() / (8 * i * chargeMacro)));
	}
	//Now I know how many particles I want to place at each radius.
	//Distribute them using an inverse transform sampling
	//But without any randomness, just equally space the uniform distribution
	rings.clear();
	rings.reserve(numMacro);//this is not exactly the number of particles we will load, but it is close
	//Use std::random to produce boltzmann distributed speeds
	//this are simply gaussian with the appropriate std deviation as we are looking at 1 dimension only
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0, sqrt(KB * aTemperature / mass));
	for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
	{
		double deltaQ{ cumulativeAtR[indexR].back() / (numAtR[indexR] + 1) };
		int currentIndex{ 0 };
		for (int i = 0; i < numAtR[indexR]; ++i)
		{
			double currentInvert{ deltaQ * (i + 1) };
			while (abs(cumulativeAtR[indexR][currentIndex]) < abs(currentInvert))
			{
				currentIndex++;
			}
			//The point to invert is between currentIndex and currentIndex - 1
			//Assume linear interpolation between this two points
			double invertedPos{ (currentIndex - 1) * refTrap.hz + refTrap.hz / 2 + refTrap.hz * (currentInvert - cumulativeAtR[indexR][currentIndex - 1]) / (cumulativeAtR[indexR][currentIndex] - cumulativeAtR[indexR][currentIndex - 1]) };
			rings.push_back(MacroRing(indexR, invertedPos, distribution(generator)));
		}
	}
	std::cout << "Loading " << rings.size() << " macro-particles from which " << numAtR[0] << " are at r=0.\n";
	solvePoisson();
}
