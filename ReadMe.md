## Alaska Herring Assessment Model
Steve Martell.
May 12, 2016
Jun 29, 2016
----

### Introduction
This Repository docments changes made to the Age-structured Model for Alaska herring stocks, VERSION 0.1, Jan 2015.  The authors of the assessment are Sherri Dressel, Sara Miller, and Kray Van Kirk.  The code in VERSION 0.1 was developed by Peter Hulson <pete.hulson@noaa.gov> .

To Do List

- [ ] Improve Numerical Stability:
	- [x] Parameter transformations to log-space
	- [x] Design matrix for estimated parameters (ival, lb, ub, phz, bayes)
	- [ ] Develop option to condition the model on F and fit to catch.

- [ ] Code reorganization:
	- [x] HEAD::initializeModelParameters
	- [x] BIOL::initializeAgeSchedule
		- [x] BIOL::getNaturalMortality
		- [x] OBSM::getSelex
	- [x] BIOL::initializeStateVariables
	- [x] BIOL::updateStateVariables
	- [x] BIOL::calcSpawningStockRecruitment
	- [x] OBSM::calcAgeCompResiduals
	- [x] OBSM::calcEggMiledaySurveyResiduals
	- [ ] OBSM::catchObservationModel
	- [ ] STAT::calcObjectiveFunction
		- [ ] STAT::penaltyFunctions
		- [ ] STAT::negativeLogLikelihoods
		- [ ] STAT::constraintFunctions
		- [ ] STAT::calculateDIC
	- [ ] FORE::

- [ ] Control file reorganization
	- [x] add design matrix to control parameter bounds and phases. 
	- [x] add contrls for time varying maturity
	- [x] add contrls for time varying natural mortality rate deviations.
	- [x] add design matrix for selectivity parameter controls.
	- [ ] add Miscellaneous controls for appending.


- [ ] Simulation model:
	- [x] Add command line argument to turn on simulation model.
	- [x] Add FUNCTION runSimulationModel
		- [x] Get and use True Parameter values
		- [x] Generate random variables | random number seed.
		- [ ] Run population dynamics model conditioned on catch & process errors.
		- [x] Cache observation model results into data variables.
		- [x] Add observation errors.
		- [x] Allow Alaska Herring Assessment to continue with parameter estimation.