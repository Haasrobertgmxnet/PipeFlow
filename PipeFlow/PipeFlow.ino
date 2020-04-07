#include <functional>
#include <cmath>

//Pipe flow calculations for 
// * pressure drop
// * mean velocity
// * friction factors
// !! This sketch currently does not build on a AVR ATmega 328P 
// 
// Sample output: see below after the loop() function
//

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
enum class CalcMethod : unsigned int { HagenPoiseuille, DarcyWeisbach, Colebrook, Churchill, Fanning, Any };

struct PipeFlow;

// BEGIN Class PipeFlow

struct PipeFlow {

	/** \brief PipeFlow Creator
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
	void Create(decimal _diameter,
		decimal _length,
		decimal _roughness,
		decimal _density,
		decimal _absoluteViscosity);

	/**<
	* \return Reynolds Number
	*
	* \param _velocity : Mean Velocity in m/s
	*/
	decimal ReynoldsNumber(decimal _velocity);

	//Geometric Data
	decimal diameter{ 0.0 };
	decimal length{ 0.0 };
	decimal roughness{ 0.0 };

	//Fluid Data
	decimal density{ 0.0 };
	decimal absoluteViscosity{ 0.0 };

	//Method data std::function<decimal(PipeFlow*, decimal, decimal)>
	struct PreferredCalcMethods {
		std::function<decimal(decimal, decimal)> PressureDifferenceLaminar;
		std::function<decimal(decimal, decimal)> PressureDifferenceTurbulent;
		std::function<decimal(decimal)> MeanVelocityLaminar;
		std::function<decimal(decimal)> MeanVelocityTurbulent;
		std::function<decimal(decimal)> FrictionFactorLaminar;
		std::function<decimal(decimal)> FrictionFactorTurbulent;
	} preferredCalcMethods;
};

// END Class PipeFlow Declarations

// BEGIN Namespace PipeFlowCalculation

namespace PipeFlowCalculation {
	//Friction Factor Formulas

	/** \brief Genral Formula of Friction Factor
	*
	*/
	/**<
	* \return Mean Friction Factor, dimensionless
	*
	* \param _Re : Reynolds Number, dimensionless
	*/
	template<CalcMethod meth>
	decimal FrictionFactor(PipeFlow* _pipeFlow, decimal _Re) {
		return 0.0;
	}

	/** \brief Friction Factor for turbulent Flow Regime from Colebrook Equation
*
*/
	template<>
	decimal FrictionFactor < CalcMethod::Colebrook>(PipeFlow* _pipeFlow, decimal _Re) {
		decimal x1{ 0.0 };
		decimal x2{ 0.0 };
		uint8_t maxiter{ 30 };
		uint8_t iter{ 0 };
		do
		{
			x1 = x2;
			x2 = -2 * log10(_pipeFlow->roughness / (3.7 * _pipeFlow->diameter) + 2.51 * x1 / _Re);
			++iter;
		} while (abs(x2 - x1) > 1e-6 && iter < maxiter);
		return 1.0 / (x2 * x2);
	}

	/** \brief Friction Factor for laminar and turbulent Flow Regime from Churchill Equation
	*
	*/
	template<>
	decimal FrictionFactor<CalcMethod::Churchill>(PipeFlow* _pipeFlow, decimal _Re) {
		if (_Re * _pipeFlow->diameter <= 0.0) {
			return 0.0;
		}
		decimal A{ pow(2.457 * log(1.0 / (pow(7.0 / _Re, 0.9) + 0.27 * _pipeFlow->roughness / _pipeFlow->diameter)), 16.0) };
		decimal B{ pow(37530.0 / _Re, 16.0) };
		decimal f{ pow((pow(8.0 / _Re, 12.0) + 1.0 / pow(A + B, 1.5)), 1.0 / 12.0) };
		return 8.0 * f;
	}

	/** \brief Friction Factor for laminar Flow Regime from Fanning Equation
	*
	*/
	template<>
	decimal FrictionFactor<CalcMethod::Fanning>(PipeFlow* _pipeFlow, decimal _Re) {
		if (_Re <= 0.0) {
			return 0.0;
		}
		return 64.0 / _Re;
	}

	// Mean Velocity Formulas

	/** \brief General Formula for Mean Velocity
	*
	*/
	/**<
	* \return Mean Velocity in m/s
	*
	* \param _deltaP : Pressure Difference in Pa
	*/
	template<CalcMethod meth>
	decimal MeanVelocity(PipeFlow* _pipeFlow, decimal _deltaP) {
		return 0.0;
	}

	/** \brief Mean Velocity from transformed Colebrook Equation
	*
	*/
	template<>
	decimal MeanVelocity<CalcMethod::Colebrook>(PipeFlow* _pipeFlow, decimal _deltaP) {
		if (_pipeFlow->diameter * _pipeFlow->length * _pipeFlow->density <= 0) {
			return -1000.0;
		}
		decimal A{ _deltaP / (_pipeFlow->length / _pipeFlow->diameter * _pipeFlow->density / 2) };
		if (A <= 0) {
			return -1000.0;
		}
		A = sqrt(A);
		auto w = -2 * A * log10(_pipeFlow->roughness / (3.7 * _pipeFlow->diameter) + 2.51 * _pipeFlow->absoluteViscosity / (A * _pipeFlow->density * _pipeFlow->diameter));
		return w;// -2 * A * log10(roughness / (3.7 * diameter) + 2.51 * absoluteViscosity / (A * density * diameter));
	}

	/** \brief Mean Velocity from transformed Hagen-Poiseuille Equation
	*
	*/
	template<>
	decimal MeanVelocity<CalcMethod::HagenPoiseuille>(PipeFlow* _pipeFlow, decimal _deltaP) {
		if (_pipeFlow->absoluteViscosity * _pipeFlow->length == 0) {
			return 0.0;
		}
		return pow(0.5 * _pipeFlow->diameter, 2.0) * _deltaP / (8.0 * _pipeFlow->absoluteViscosity * _pipeFlow->length);
	}

	template<>
	decimal MeanVelocity<CalcMethod::Any>(PipeFlow* _pipeFlow, decimal _deltaP) {
		if (_deltaP <= 0.0) {
			return 0.0;
		}
		decimal threshold{ 0.25 };//Is this a good value?
		//Try with turbulent flow Regime
		auto MeanVelocity = _pipeFlow->preferredCalcMethods.MeanVelocityTurbulent;
		decimal velocity = _pipeFlow->preferredCalcMethods.MeanVelocityTurbulent(_deltaP);
		decimal Re = _pipeFlow->ReynoldsNumber(velocity);
		if (Re > ReKrit) {
			return velocity;
		}
		//Try with laminar flow Regime
		MeanVelocity = _pipeFlow->preferredCalcMethods.MeanVelocityLaminar;
		velocity = MeanVelocity(_deltaP);
		Re = _pipeFlow->ReynoldsNumber(velocity);
		if (Re <= ReKrit) {
			return velocity;
		}
		return 0.0;
	}

	// Pressure Difference Formulas

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
	decimal PressureDifference(PipeFlow* _pipeFlow, decimal _velocity) {
		return 0.0;
	}

	/** \brief Pressure Difference for lamainar Flow Regime using Hagen-Poiseuille Formula
	* no _frictionFactor needed
	*/
	template<>
	decimal PressureDifference<CalcMethod::HagenPoiseuille>(PipeFlow* _pipeFlow, decimal _velocity) {
		if (_pipeFlow->diameter == 0) {
			return 0.0;
		}
		return 8 * _pipeFlow->absoluteViscosity * _pipeFlow->length * _velocity / pow(0.5 * _pipeFlow->diameter, 2.0);
	}

	template<CalcMethod meth>
	decimal PressureDifference(PipeFlow* _pipeFlow, decimal _velocity, decimal _frictionFactor) {
		return 0.0;
	}

	/** \brief Pressure Difference for lamainar und turbulent Flow Regime using Darcy-Weisbach Formula
	* _frictionFactor mandatory
	*/
	template<>
	decimal PressureDifference<CalcMethod::DarcyWeisbach>(PipeFlow* _pipeFlow, decimal _velocity, decimal _frictionFactor) {
		if (_pipeFlow->diameter == 0) {
			return 0.0;
		}
		return _frictionFactor * (_pipeFlow->length / _pipeFlow->diameter) * _pipeFlow->density * _velocity * _velocity / 2.0;
	}

	template<>
	decimal PressureDifference<CalcMethod::Any>(PipeFlow* _pipeFlow, decimal _velocity) {
		decimal Re{ _pipeFlow->ReynoldsNumber(_velocity) };
		if (Re <= ReKrit) {
			//Calculations
			decimal f{ FrictionFactor< CalcMethod::Fanning >(_pipeFlow, Re) };

			return PressureDifference<CalcMethod::DarcyWeisbach>(_pipeFlow, _velocity, f);
		}
		if (Re > ReKrit) {
			decimal f{ FrictionFactor< CalcMethod::Churchill >(_pipeFlow, Re) };
			return PressureDifference<CalcMethod::DarcyWeisbach>(_pipeFlow, _velocity, f);
		}
		return 0.0;
	}
}

// END Namespace PipeFlowCalculation

// BEGIN Class PipeFlow Definitions

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
void PipeFlow::Create(decimal _diameter,
	decimal _length,
	decimal _roughness,
	decimal _density,
	decimal _absoluteViscosity)
{
	diameter = _diameter;
	length = _length;
	roughness = _roughness;
	density = _density;
	absoluteViscosity = _absoluteViscosity;

	preferredCalcMethods.PressureDifferenceLaminar = [this](decimal _velocity, decimal _frictionFactor) {
		return PipeFlowCalculation::PressureDifference<CalcMethod::DarcyWeisbach>(this, _velocity, _frictionFactor);
	};
	preferredCalcMethods.PressureDifferenceTurbulent = [this](decimal _velocity, decimal _frictionFactor) {
		return PipeFlowCalculation::PressureDifference<CalcMethod::DarcyWeisbach>(this, _velocity, _frictionFactor);
	};
	preferredCalcMethods.MeanVelocityLaminar = [this](decimal _deltaP) {
		return PipeFlowCalculation::MeanVelocity<CalcMethod::HagenPoiseuille>(this, _deltaP);
	};
	preferredCalcMethods.MeanVelocityTurbulent = [this](decimal _deltaP) {
		return PipeFlowCalculation::MeanVelocity<CalcMethod::Colebrook>(this, _deltaP);
	};
	preferredCalcMethods.FrictionFactorLaminar = [this](decimal _Re) {
		return PipeFlowCalculation::FrictionFactor<CalcMethod::Fanning>(this, _Re);
	};
	preferredCalcMethods.FrictionFactorTurbulent = [this](decimal _Re) {
		return PipeFlowCalculation::FrictionFactor<CalcMethod::Churchill>(this, _Re);
	};
}

/**<
* \return Reynolds Number
*
* \param _velocity : Mean Velocity in m/s
*/
decimal PipeFlow::ReynoldsNumber(decimal _velocity) {
	if (absoluteViscosity == 0) {
		return -1.0;
	}
	return density * _velocity * diameter / absoluteViscosity;
}

// END Class PipeFlow Definitions

static size_t count = 0;

// the setup function runs once when you press reset or power the board
void setup() {
    Serial.begin(115200);
	{
		Serial.println("Example with dry air.");
		Serial.println("Density and dynamic viscosity are from");
		Serial.println("VDI Heat Atlas, 2nd Edition.");

		//Geometric Data
		decimal diameter = 0.01;
		decimal length = 0.05;
		decimal roughness = 1e-6;

		//Fluid Data (These are for dry air (from VDI Heat Atlas))
		decimal density = 1.1124;//kg / m3
		decimal absoluteViscosity = 19.165e-6;// Pa * s

		PipeFlow pipeFlow;
		pipeFlow.Create(diameter, length, roughness, density, absoluteViscosity);

		decimal area = Pi * (diameter / 2.0) * (diameter / 2.0);
		decimal w0 = 0.0;
		//Case 1: deltaP is given
		//Volume flow calculated
		{
			Serial.println(" ");
			Serial.println("(1) Mean velocity and volume flow are calculated from pressure difference.");
			Serial.println(" ");

			decimal deltaP = 19.9993;
			decimal w = PipeFlowCalculation::MeanVelocity<CalcMethod::Any>(&pipeFlow, deltaP);
			decimal Re = pipeFlow.ReynoldsNumber(w);
			decimal vFlow = area * w;
			decimal f = deltaP / (length * density * w * w / 2 / diameter);

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

			decimal w = w0;
			decimal Re = pipeFlow.ReynoldsNumber(w);
			Serial.print("Mean Velocity (m/s) as calculated (1): ");
			Serial.println(w, 4);
			Serial.print("Calculated Reynolds Number (.): ");
			Serial.println(Re, 4);
			if (Re <= ReKrit) {

				Serial.println(" ");
				Serial.println("Laminar flow regime:");
				Serial.println(" ");

				//Calculations
				decimal vFlow = area * w;//Volume flow from velocity
				decimal f1 = PipeFlowCalculation::FrictionFactor<CalcMethod::Fanning>(&pipeFlow, Re);
				decimal f2 = PipeFlowCalculation::FrictionFactor < CalcMethod::Churchill>(&pipeFlow,Re);
				decimal deltaP1 = PipeFlowCalculation::PressureDifference<CalcMethod::DarcyWeisbach>(&pipeFlow,w,f1);
				decimal deltaP2 = PipeFlowCalculation::PressureDifference<CalcMethod::DarcyWeisbach>(&pipeFlow,w,f2);
				decimal deltaP3 = PipeFlowCalculation::PressureDifference<CalcMethod::HagenPoiseuille>(&pipeFlow,w);

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
				decimal f1 = PipeFlowCalculation::FrictionFactor<CalcMethod::Colebrook>(&pipeFlow, Re);
				decimal f2 = PipeFlowCalculation::FrictionFactor<CalcMethod::Churchill>(&pipeFlow, Re);
				decimal deltaP1 = PipeFlowCalculation::PressureDifference<CalcMethod::DarcyWeisbach>(&pipeFlow, w, f1);
				decimal deltaP2 = PipeFlowCalculation::PressureDifference<CalcMethod::DarcyWeisbach>(&pipeFlow, w, f2);

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
//Calculated Delta p with Churchill friction factor(Pa) : 20.1202