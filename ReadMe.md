## Alaska Herring Assessment Model
Steve Martell.
May 12, 2015

----

### Introduction
This Repository docments changes made to the Age-structured Model for Alaska herring stocks, VERSION 0.1, Jan 2015.  The authors of the assessment are Sherri Dressel, Sara Miller, and Kray Van Kirk.  The code in VERSION 0.1 was developed by Peter Hulson <pete.hulson@noaa.gov> .

To Do List

- [ ] Improve Numerical Stability:
	- [ ] Parameter transformations to log-space
	- [ ] Design matrix for estimated parameters (ival, lb, ub, phz, bayes)

- [ ] Code reorganization:
	- [ ] HEAD::initializeModelParameters
	- [ ] BIOL::initializeAgeSchedule
		- [ ] BIOL::getNaturalMortality
		- [ ] OBSM::getSelex
	- [ ] BIOL::initializeStateVariables
	- [ ] BIOL::updateStateVariables
	- [ ] OBSM::
	- [ ] OBSM::catchObservationModel

- [ ] Simulation model:
	- [ ] Add command line argument to turn on simulation model.
	- [ ] Add FUNCTION runSimulationModel
		- [ ] Get and use True Parameter values
		- [ ] Generate random variables | random number seed.
		- [ ] Run population dynamics model conditioned on catch & process errors.
		- [ ] Cache observation model results into data variables.
		- [ ] Add observation errors.
		- [ ] Allow Alaska Herring Assessment to continue with parameter estimation.