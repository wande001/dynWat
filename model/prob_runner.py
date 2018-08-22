#!/usr/bin/env python
# -*- coding: utf-8 -*-

########PCRGLOB-WB with data assimilation######
#### Author N.Wanders (n.wanders@uu.nl)#####
#### Version 1.0 (25-08-2014) ###########

import sys
import logging
logger = logging.getLogger(__name__)

from pcraster.framework import DynamicModel
from pcraster.framework import MonteCarloModel
from pcraster.framework import DynamicFramework

from configuration import Configuration
from currTimeStep import ModelTime
from reporting import Reporting
from spinUp import SpinUp

from pcrglobwb import PCRGlobWB

from EnKF import *

import pcraster as pcr
import virtualOS as vos
import datetime
import numpy

logger = logging.getLogger(__name__)

class DeterministicRunner(DynamicModel, MonteCarloModel):

    def __init__(self, configuration, modelTime, initialState = None):
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)

        self._configuration = configuration
        self.modelTime = modelTime        
        self.model = PCRGlobWB(configuration, modelTime, initialState)
        self.reporting = Reporting(configuration, self.model, modelTime)
        self.cmdOptions = sys.argv
        print self.cmdOptions

        self.ReportTime = self._configuration.dataAssimilationOptions['filterTimeSteps']
        self.ReportTime = map(int, self.ReportTime.split(','))
        self.ReportTime.append(0)
        
        self.shortNames = ['f','g','p','n']
                
    def premcloop(self):

        pass
      
    def initial(self):
        pass

        
    def dynamic(self):
        #re-calculate current model time using current pcraster timestep value
        self.modelTime.update(self.currentTimeStep())

        #update model (will pick up current model time from model time object)
        
        self.model.read_forcings_forecasting()
        
        msg = "Date = " + str(self.modelTime.currTime.day)+ "-"+ str(self.modelTime.currTime.month)+ " " + str(self.currentSampleNumber())
        logger.info(msg)
        print self.model.landSurface.coverTypes
        self.model.update(report_water_balance=False)
        
        msg = 'Model update succesfull' + str(self.currentSampleNumber())
        logger.info(msg)

        #do any needed reporting for this time step        
        self.reporting.report()
        msg = 'Reporting succesfull' + str(self.currentSampleNumber())
        logger.info(msg)

        # writing the model states to disk 
        # - that will be re-used in the "resume" method:
          

def main():
    initial_state = None
    configuration = Configuration()
    
    spin_up = SpinUp(configuration)                   # object for spin_up
    
    currTimeStep = ModelTime() # timeStep info: year, month, day, doy, hour, etc
    
    # spinningUp
    noSpinUps = int(configuration.globalOptions['maxSpinUpsInYears'])
    if noSpinUps > 0:
        
        logger.info('Spin-Up #Total Years: '+str(noSpinUps))

        spinUpRun = 0 ; has_converged = False
        while spinUpRun < noSpinUps and has_converged == False:
            spinUpRun += 1
            currTimeStep.getStartEndTimeStepsForSpinUp(
                    configuration.globalOptions['startTime'],
                    spinUpRun, noSpinUps)
            logger.info('Spin-Up Run No. '+str(spinUpRun))
            deterministic_runner = DeterministicRunner(configuration, currTimeStep, initial_state)
            
            all_state_begin = deterministic_runner.model.getAllState() 
            
            dynamic_framework = DynamicFramework(deterministic_runner,currTimeStep.nrOfTimeSteps)
            dynamic_framework.setQuiet(True)
            dynamic_framework.run()
            
            all_state_end = deterministic_runner.model.getAllState() 
            
            has_converged = spin_up.checkConvergence(all_state_begin, all_state_end, spinUpRun, deterministic_runner.model.routing.cellArea)
            
            initial_state = deterministic_runner.model.getState()
    #
    # Running the deterministic_runner (excluding DA scheme)
    currTimeStep.getStartEndTimeSteps(configuration.globalOptions['startTime'],
                                      configuration.globalOptions['endTime'])
    
    logger.info('Transient simulation run started.')
    deterministic_runner = DeterministicRunner(configuration, currTimeStep, initial_state)
    
    dynamic_framework = DynamicFramework(deterministic_runner,currTimeStep.nrOfTimeSteps)
    dynamic_framework.setQuiet(True)
    if configuration.dataAssimilationOptions['method'] == "None":
        dynamic_framework.run()
    else:
        nrSamples = int(configuration.dataAssimilationOptions['nrSamples'])
        mcModel = MonteCarloFramework(dynamic_framework,nrSamples)
        mcModel.setForkSamples(True, nrCPUs=int(configuration.dataAssimilationOptions['nrCores']))
    if configuration.dataAssimilationOptions['method'] == "MonteCarlo":
        mcModel.run()
    if configuration.dataAssimilationOptions['method'] == "EnKF":
        ekfModel = EnsKalmanFilterFramework(mcModel)
        filterTime = configuration.dataAssimilationOptions['filterTimeSteps']
        filterTime = map(int, filterTime.split(','))
        ekfModel.setFilterTimesteps(filterTime)    #range(365,6900,30)
        ekfModel.run()

        
if __name__ == '__main__':
    sys.exit(main())

