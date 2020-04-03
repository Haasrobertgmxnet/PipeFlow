//Pipe flow calculations for 
// * pressure drop
// * mean velocity
// * friction factors
// !! This sketch currently does not build on a AVR ATmega 328P 
// 
// Sample output: see below after the loop() function
//
#include <cmath>

using decimal = double;
const decimal Pi = 3.1415926535;
const decimal ReKrit = 2300.0;

// The basic Calculation methods:
// * Hagen-Poiseuille: fully developed, stationary and laminar flow. https://de.wikipedia.org/wiki/Darcy-Weisbach-Gleichung
// * Darcy-Weisbach: laminar and turbulent flow. https://de.wikipedia.org/wiki/Darcy-Weisbach-Gleichung
// * Colebrook: friction factor for turbulent flow: https://de.wikipedia.org/wiki/Rohrreibungszahl
// * Churchill: friction factor for laminar and turbulent flow: https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Table_of_Approximations
// * Fanning: friction factor for laminar flow: https://en.wikipedia.org/wiki/Fanning_friction_factor
//
enum class CalcMethod : unsigned int { HagenPoiseuille, DarcyWeisbach, Colebrook, Churchill, Fanning };

// BEGIN Class PipeFlow

struct PipeFlow {

	/** \brief PipeFlow Constructor
	*
	*/
	/**<
	*
	* \param _diameter : Diameter in m
	* \param _length : Length in m
	* \param _roughness : Surface Roughness in m
	* \param _density : Density in kg/m3
	* \param _absoluteViscosity : Dynamic Visosity in Pa*s
	*/
	PipeFlow(decimal _diameter,
		decimal _length,
		decimal _roughness,
		decimal _density,
		decimal _absoluteViscosity) :
		diameter{ _diameter },
		length{ _length },
		roughness{ _roughness },
		density{ _density },
		absoluteViscosity{ _absoluteViscosity }
	{}

	/** \brief General Formula for Pressure Difference
	*
	*/
	/**<
	* \return Pressure Difference in Pa
	*
	* \param _velocity : Mean Velocity in m/s
	* \param _frictionFactor : Friction Factor, dimensionless
	*/
	template<CalcMethod meth>
	decimal PressureDifference(decimal _velocity, decimal _frictionFactor = 0.0) {
		return 0.0;
	}

	/** \brief General Formula for Mean Velocity
	*
	*/
	/**<
	* \return Mean Velocity in m/s
	*
	* \param _deltaP : Pressure Difference in Pa
	*/
	template<CalcMethod meth>
	decimal MeanVelocity(decimal _deltaP) {
		return 0.0;
	}

	/** \brief Genral Formula of Friction Factor
	*
	*/
	/**<
	* \return Mean Friction Factor, dimensionless
	*
	* \param _Re : Reynolds Number, dimensionless
	*/
	template<CalcMethod meth>
	decimal FrictionFactor(decimal _Re) {
		return 0.0;
	}

	/**<
	* \return Reynolds Number
	*
	* \param _velocity : Mean Velocity in m/s
	*/
	decimal ReynoldsNumber(decimal _velocity) {
		if (absoluteViscosity == 0) {
			return -1.0;
		}
		return density * _velocity * diameter / absoluteViscosity;
	}

	//Geometric Data
	decimal diameter{ 0.0 };
	decimal length{ 0.0 };
	decimal roughness{ 0.0 };

	//Fluid Data
	decimal density{ 0.0 };
	decimal absoluteViscosity{ 0.0 };

	//Method data
	struct PreferredCalcMethods {
		CalcMethod PressureDifferenceLaminar{ CalcMethod::DarcyWeisbach };
		CalcMethod PressureDifferenceTurbulent{ CalcMethod::DarcyWeisbach };
		CalcMethod MeanVelocityLaminar{ CalcMethod::HagenPoiseuille };
		CalcMethod MeanVelocityTurbulent{ CalcMethod::Colebrook };
		CalcMethod FrictionFactorLaminar{ CalcMethod::Fanning };
		CalcMethod FrictionFactorTurbulent{ CalcMethod::Churchill };
	} preferredCalcMethods;
};

// END Class PipeFlow

// BEGIN Template specializations

// Pressure Difference Formulas

/** \brief Pressure Difference for lamainar Flow Regime using Hagen-Poiseuille Formula
* no _frictionFactor needed
*/
template<>
decimal PipeFlow::PressureDifference<CalcMethod::HagenPoiseuille>(decimal _velocity, decimal _frictionFactor) {
	if (diameter == 0) {
		return 0.0;
	}
	return 8 * absoluteViscosity * length * _velocity / pow(0.5 * diameter, 2.0);
}

/** \brief Pressure Difference for lamainar und turbulent Flow Regime using Darcy-Weisbach Formula
* _frictionFactor mandatory
*/
template<>
decimal PipeFlow::PressureDifference<CalcMethod::DarcyWeisbach>(decimal _velocity, decimal _frictionFactor) {
	if (diameter == 0) {
		return 0.0;
	}
	return _frictionFactor * (length / diameter) * density * _velocity * _velocity / 2.0;
}

// Mean Velocity Formulas

	/** \brief Mean Velocity from transformed Colebrook Equation
	*
	*/
template<>
decimal PipeFlow::MeanVelocity<CalcMethod::Colebrook>(decimal _deltaP) {
	if (diameter * length * density <= 0) {
		return -1000.0;
	}
	decimal A{ _deltaP / (length / diameter * density / 2) };
	if (A <= 0) {
		return -1000.0;
	}
	A = sqrt(A);
	return -2 * A * log10(roughness / (3.7 * diameter) + 2.51 * absoluteViscosity / (A * density * diameter));
}

/** \brief Mean Velocity from transformed Hagen-Poiseuille Equation
*
*/
template<>
decimal PipeFlow::MeanVelocity<CalcMethod::HagenPoiseuille>(decimal _deltaP) {
	if (absoluteViscosity * length == 0) {
		return 0.0;
	}
	return pow(0.5 * diameter, 2.0) * _deltaP / (8.0 * absoluteViscosity * length);
}

//Friction Factor Formulas

	/** \brief Friction Factor for turbulent Flow Regime from Colebrook Equation
	*
	*/
template<>
decimal PipeFlow::FrictionFactor < CalcMethod::Colebrook>(decimal _Re) {
	decimal x1{ 0.0 };
	decimal x2{ 0.0 };
	uint8_t maxiter{ 30 };
	uint8_t iter{ 0 };
	do
	{
		x1 = x2;
		x2 = -2 * log10(roughness / (3.7 * diameter) + 2.51 * x1 / _Re);
		++iter;
	} while (abs(x2 - x1) > 1e-6 && iter < maxiter);
	return 1.0 / (x2 * x2);
}

/** \brief Friction Factor for laminar and turbulent Flow Regime from Churchill Equation
*
*/
template<>
decimal PipeFlow::FrictionFactor<CalcMethod::Churchill>(decimal _Re) {
	if (_Re * diameter <= 0.0) {
		return 0.0;
	}
	decimal A{ pow(2.457 * log(1.0 / (pow(7.0 / _Re, 0.9) + 0.27 * roughness / diameter)), 16.0) };
	decimal B{ pow(37530.0 / _Re, 16.0) };
	decimal f{ pow((pow(8.0 / _Re, 12.0) + 1.0 / pow(A + B, 1.5)), 1.0 / 12.0) };
	return 8.0 * f;
}

/** \brief Friction Factor for laminar Flow Regime from Fanning Equation
*
*/
template<>
decimal PipeFlow::FrictionFactor<CalcMethod::Fanning>(decimal _Re) {
	if (_Re <= 0.0) {
		return 0.0;
	}
	return 64.0 / _Re;
}

// END Template specializations

static size_t count = 0;

// the setup function runs once when you press reset or power the board
void setup() {
    Serial.begin(115200);
	{
		Serial.println("Example with dry air.");
		Serial.println("Density and dynamic viscosity are from");
		Serial.println("VDI Heat Atlas, 2nd Edition.");

		//Geometric Data
		decimal diameter{ 0.01 };
		decimal length{ 0.05 };
		decimal roughness{ 1e-6 };

		//Fluid Data (These are for dry air (from VDI Heat Atlas))
		decimal density{ 1.1124 };//kg / m3
		decimal absoluteViscosity{ 19.165e-6 };// Pa * s

		PipeFlow pipeFlow(diameter, length, roughness, density, absoluteViscosity);

		decimal area{ Pi * (diameter / 2.0) * (diameter / 2.0) };
		decimal w0{ 0.0 };
		//Case 1: deltaP is given
		//Volume flow calculated
		{
			Serial.println(" ");
			Serial.println("(1) Mean velocity and volume flow are calculated from pressure difference.");
			Serial.println(" ");

			decimal deltaP{ 19.9993 };
			decimal w{ pipeFlow.MeanVelocity<CalcMethod::Colebrook>(deltaP) };
			decimal Re{ pipeFlow.ReynoldsNumber(w) };
			decimal vFlow{ area * w };
			decimal f{ deltaP / (length * density * w * w / 2 / diameter) };

			Serial.print("Given Pressure Difference (Pa): ");
			Serial.println(deltaP,4);
			Serial.print("Calculated Reynolds Number (.): ");
			Serial.println(Re, 4);
			Serial.print("Calculated Mean Velocity (m/s): ");
			Serial.println(w, 4);
			Serial.print("Volume Flow (l/min) from Mean Velocity: ");
			Serial.println(6e4 * vFlow, 4);
			Serial.print("Calculated Friction factor (.): ");
			Serial.println(f, 4);

			//save for Case 2
			w0 = w;
		}
		//Case 2: Mean velocity is given
		//deltaP calculated
		{
			Serial.println(" ");
			Serial.println("(2) The pressure difference is calculated from mean velocity from (1).");
			Serial.println(" ");

			decimal w{ w0 };
			decimal Re{ pipeFlow.ReynoldsNumber(w) };
			Serial.print("Mean Velocity (m/s) as calculated (1): ");
			Serial.println(w, 4);
			Serial.print("Calculated Reynolds Number (.): ");
			Serial.println(Re, 4);
			if (Re <= ReKrit) {

				Serial.println(" ");
				Serial.println("Laminar flow regime:");
				Serial.println(" ");

				//Calculations
				decimal vFlow{ area * w };//Volume flow from velocity
				decimal f1{ pipeFlow.FrictionFactor<CalcMethod::Fanning>(Re) };
				decimal f2{ pipeFlow.FrictionFactor < CalcMethod::Churchill>(Re) };
				decimal deltaP1{ pipeFlow.PressureDifference<CalcMethod::DarcyWeisbach>(w,f1) };
				decimal deltaP2{ pipeFlow.PressureDifference<CalcMethod::DarcyWeisbach>(w,f2) };
				decimal deltaP3{ pipeFlow.PressureDifference<CalcMethod::HagenPoiseuille>(w) };

				//Output
				Serial.print("Calculated Fanning friction factor (.): ");
				Serial.println(f1, 4);
				Serial.print("Calculated Churchill friction factor (.): ");
				Serial.println(f2, 4);
				Serial.print("Calculated Delta p with Fanning friction factor (Pa): ");
				Serial.println(deltaP1, 4);
				Serial.print("Calculated Delta p with Churchill friction factor (Pa): ");
				Serial.println(deltaP2, 4);
				Serial.print("Delta p Hagen-Poiseuille (Pa): ");
				Serial.println(deltaP3, 4);
			}
			if (Re > ReKrit) {

				Serial.println(" ");
				Serial.println("Turbulent flow regime:");
				Serial.println(" ");

				//Calculations
				decimal f1{ pipeFlow.FrictionFactor<CalcMethod::Colebrook>(Re) };
				decimal f2{ pipeFlow.FrictionFactor<CalcMethod::Churchill>(Re) };
				decimal deltaP1{ pipeFlow.PressureDifference<CalcMethod::DarcyWeisbach>(w,f1) };
				decimal deltaP2{ pipeFlow.PressureDifference<CalcMethod::DarcyWeisbach>(w,f2) };

				//Output
				Serial.print("Calculated Colebrook friction factor (.): ");
				Serial.println(f1, 4);
				Serial.print("Calculated Churchill friction factor (.): ");
				Serial.println(f2, 4);
				Serial.print("Calculated Delta p with Colebrook friction factor (Pa): ");
				Serial.println(deltaP1, 4);
				Serial.print("Calculated Delta p with Churchill friction factor (Pa): ");
				Serial.println(deltaP2, 4);
			}
			Serial.println("Case 2 done: ");
		}
		Serial.println("Calculation done: ");
		Serial.println(" ");
	}
}


// the loop function runs over and over again forever
void loop() {
	Serial.begin(115200);
    Serial.print("Calculation done: ");
    Serial.print(count++);
	delay(5000);
}

//Sample output

//Example with dry air.
//Densityand dynamic viscosity are from
//VDI Heat Atlas, 2nd Edition.
//
//(1) Mean velocity and volume flow are calculated from pressure difference.
//
//Given Pressure Difference(Pa) : 19.9993
//Calculated Reynolds Number(.) : 8670.7337
//Calculated Mean Velocity(m / s) : 14.9384
//Volume Flow(l / min) from Mean Velocity : 70.3955
//Calculated Friction factor(.) : 0.0322
//
//(2) The pressure difference is calculated from mean velocity from(1).
//
//Mean Velocity(m / s) as calculated(1) : 14.9384
//Calculated Reynolds Number(.) : 8670.7337
//
//Turbulent flow regime :
//
//Calculated Colebrook friction factor(.) : 0.0322
//Calculated Churchill friction factor(.) : 0.0324
//Calculated Delta p with Colebrook friction factor(Pa) : 19.9993
