## Alaska Herring Assessment Model
Steve Martell.
May 12, 2016
Jun 29, 2016
Jul 20, 2016
----

## Introduction
This Repository docments changes made to the Age-structured Model for Alaska herring stocks, VERSION 0.1, Jan 2015.  The authors of the assessment are Sherri Dressel, Sara Miller, and Kray Van Kirk.  The code in VERSION 0.1 was developed by Peter Hulson <pete.hulson@noaa.gov> .

## Installation
You can obtain the source code via cloning this project using Git, or you may also download a zip file with the latest code. 


Go to the [Github code repository for Ham](https://github.com/seastateinc/AlaskaHerring), and click on the Clone or download button.

![Select the Download Zip option.](https://github.com/seastateinc/AlaskaHerring/blob/develop/docs/CloneZip.png)

I would recommend using Git and Github, as they serve as valuable tools for version control, and tracking code changes over time.  Its perfectly fine to just download the zip file; however, using Git is far more efficient for incorporating code changes and working in groups.



----

## To Do List

- [x] Improve Numerical Stability:
	- [x] Parameter transformations to log-space
	- [x] Design matrix for estimated parameters (ival, lb, ub, phz, bayes)
	- [x] Develop option to condition the model on F and fit to catch.

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
	- [x] OBSM::catchObservationModel
	- [x] STAT::calcObjectiveFunction
		- [x] STAT::penaltyFunctions
		- [x] STAT::negativeLogLikelihoods
		- [x] STAT::constraintFunctions
		- [x] STAT::calculateDIC
	- [ ] FORE::runForecast
		- [ ] FORE::ghlCalc
		- [ ] FORE::calcTAC

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