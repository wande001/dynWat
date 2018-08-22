#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from pcraster.framework import *
import pcraster as pcr

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos
from ncConverter import *
import ETPFunctions as refPotET
import scipy.stats

class Meteo(object):

    def __init__(self,iniItems,landmask,spinUp):
        object.__init__(self)

        self.cloneMap = iniItems.cloneMap
        self.tmpDir = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
        
        # landmask/area of interest
        self.landmask = landmask
        if iniItems.globalOptions['landmask'] != "None":
           self.landmask = vos.readPCRmapClone(\
           iniItems.globalOptions['landmask'],
           self.cloneMap,self.tmpDir,self.inputDir)	

        self.preFileNC = iniItems.meteoOptions['precipitationNC']        # starting from 19 Feb 2014, we only support netcdf input files
        self.tmpFileNC = iniItems.meteoOptions['temperatureNC']

        self.refETPotMethod = iniItems.meteoOptions['referenceETPotMethod']
        if self.refETPotMethod == 'Hamon': self.latitudes = \
                                           pcr.ycoordinate(self.cloneMap) # needed to calculate 'referenceETPot'
        if self.refETPotMethod == 'Hargreaves': self.latitudes = \
                                           pcr.ycoordinate(self.cloneMap) # needed to calculate 'referenceETPot'
        if self.refETPotMethod == 'Input': self.etpFileNC = \
                             iniItems.meteoOptions['refETPotFileNC']              

        try:
            self.precipitationCorrectionFactor = float(iniItems.meteoOptions['precipitationCorrectionFactor'])
        except:
            self.precipitationCorrectionFactor = float(1.0)

        try:
            self.precipitationVarName = iniItems.meteoOptions['precipitationVarName']
        except:
            self.precipitationVarName = "precipitation"

        try:
            self.temperatureVarName = iniItems.meteoOptions['temperatureVarName']
        except:
            self.temperatureVarName = "temperature"
        try:
            self.evaporationVarName = iniItems.meteoOptions['evaporationVarName']
        except:
            self.evaporationVarName = "evaporation"            
        ### Read CDF for CDF matching seasonal forecasting
        try:
            self.matchingCDF = iniItems.meteoOptions['matchingCDF']
        except:
            self.matchingCDF = False
        if self.matchingCDF:
            self.matchingCDF = True
            try:
                self.filePrecipitationReferenceCDF = iniItems.meteoOptions['precipitationReferenceCDF']
            except:
                self.fileReferenceCDF = None
                self.matchingCDF = False
                logger.info("No precipitationReferenceCDF filepath provided.")
            try:
                self.fileTemperatureReferenceCDF = iniItems.meteoOptions['temperatureReferenceCDF']
            except:
                self.fileTemperatureReferenceCDF = None
                self.matchingCDF = False
                logger.info("No temperatureReferenceCDF filepath provided.")
            try:
                self.filePrecipitationInputCDF = iniItems.meteoOptions['precipitationInputCDF']
            except:
                self.filePrecipitationInputCDF = None
                self.matchingCDF = False
                logger.info("No precipitationInputCDF filepath provided.")
            try:
                self.fileTemperatureInputCDF = iniItems.meteoOptions['temperatureInputCDF']
            except:
                self.fileTemperatureInputCDF = None
                self.matchingCDF = False
                logger.info("No temperatureInputCDF filepath provided.")
        else:
            self.matchingCDF = False
        print self.matchingCDF
        try:
            self.nameModel = iniItems.globalOptions['model']
        except:
            self.nameModel = None
        self.randomTemp = 0.5
        self.randomPrec = 0.5
            
        # daily time step
        self.usingDailyTimeStepForcingData = False
        if iniItems.timeStep == 1.0 and iniItems.timeStepUnit == "day":
            self.usingDailyTimeStepForcingData = True
        
        # forcing downscaling options:
        self.forcingDownscalingOptions(iniItems)

        self.report = True
        try:
            self.outDailyTotNC = iniItems.meteoOptions['outDailyTotNC'].split(",")
            self.outMonthTotNC = iniItems.meteoOptions['outMonthTotNC'].split(",")
            self.outMonthAvgNC = iniItems.meteoOptions['outMonthAvgNC'].split(",")
            self.outMonthEndNC = iniItems.meteoOptions['outMonthEndNC'].split(",")
            self.outAnnuaTotNC = iniItems.meteoOptions['outAnnuaTotNC'].split(",")
            self.outAnnuaAvgNC = iniItems.meteoOptions['outAnnuaAvgNC'].split(",")
            self.outAnnuaEndNC = iniItems.meteoOptions['outAnnuaEndNC'].split(",")
        except:
            self.report = False
        if self.report == True:
            # daily output in netCDF files:
            self.outNCDir  = iniItems.outNCDir
            self.netcdfObj = PCR2netCDF(iniItems)
            #
            if self.outDailyTotNC[0] != "None":
                for var in self.outDailyTotNC:
                    # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_dailyTot.nc",\
                                                    var,"undefined")
            # MONTHly output in netCDF files:
            # - cummulative
            if self.outMonthTotNC[0] != "None":
                for var in self.outMonthTotNC:
                    # initiating monthlyVarTot (accumulator variable):
                    vars(self)[var+'MonthTot'] = None
                    # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_monthTot.nc",\
                                                    var,"undefined")
            # - average
            if self.outMonthAvgNC[0] != "None":
                for var in self.outMonthAvgNC:
                    # initiating monthlyTotAvg (accumulator variable)
                    vars(self)[var+'MonthTot'] = None
                    # initiating monthlyVarAvg:
                    vars(self)[var+'MonthAvg'] = None
                     # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_monthAvg.nc",\
                                                    var,"undefined")
            # - last day of the month
            if self.outMonthEndNC[0] != "None":
                for var in self.outMonthEndNC:
                     # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_monthEnd.nc",\
                                                    var,"undefined")
            # YEARly output in netCDF files:
            # - cummulative
            if self.outAnnuaTotNC[0] != "None":
                for var in self.outAnnuaTotNC:
                    # initiating yearly accumulator variable:
                    vars(self)[var+'AnnuaTot'] = None
                    # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_annuaTot.nc",\
                                                    var,"undefined")
            # - average
            if self.outAnnuaAvgNC[0] != "None":
                for var in self.outAnnuaAvgNC:
                    # initiating annualyVarAvg:
                    vars(self)[var+'AnnuaAvg'] = None
                    # initiating annualyTotAvg (accumulator variable)
                    vars(self)[var+'AnnuaTot'] = None
                     # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_annuaAvg.nc",\
                                                    var,"undefined")
            # - last day of the year
            if self.outAnnuaEndNC[0] != "None":
                for var in self.outAnnuaEndNC:
                     # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_annuaEnd.nc",\
                                                    var,"undefined")


    def forcingDownscalingOptions(self, iniItems):

        self.downscalePrecipitationOption  = False
        self.downscaleTemperatureOption    = False
        self.downscaleReferenceETPotOption = False

        if 'meteoDownscalingOptions' in iniItems.allSections:

            # downscaling options
            if iniItems.meteoDownscalingOptions['downscalePrecipitation']  == "True":
                self.downscalePrecipitationOption  = True  
                logger.info("Precipitation forcing will be downscaled to the cloneMap resolution.")

            if iniItems.meteoDownscalingOptions['downscaleTemperature']    == "True":
                self.downscaleTemperatureOption    = True  
                logger.info("Temperature forcing will be downscaled to the cloneMap resolution.")

            if iniItems.meteoDownscalingOptions['downscaleReferenceETPot'] == "True" and self.refETPotMethod != 'Hamon':
                self.downscaleReferenceETPotOption = True 
                logger.info("Reference potential evaporation will be downscaled to the cloneMap resolution.")

            # Note that for the Hamon method: referencePotET will be calculated based on temperature,  
            # therefore, we do not have to downscale it (particularly if temperature is already provided at high resolution). 

        if self.downscalePrecipitationOption or\
           self.downscaleTemperatureOption   or\
           self.downscaleReferenceETPotOption:

            # creating anomaly DEM
            highResolutionDEM = vos.readPCRmapClone(\
               iniItems.meteoDownscalingOptions['highResolutionDEM'],
               self.cloneMap,self.tmpDir,self.inputDir)
            highResolutionDEM = pcr.cover(highResolutionDEM, 0.0)
            highResolutionDEM = pcr.max(highResolutionDEM, 0.0)
            self.meteoDownscaleIds = vos.readPCRmapClone(\
               iniItems.meteoDownscalingOptions['meteoDownscaleIds'],
               self.cloneMap,self.tmpDir,self.inputDir,isLddMap=False,cover=None,isNomMap=True)
            self.cellArea = vos.readPCRmapClone(\
               iniItems.routingOptions['cellAreaMap'],
               self.cloneMap,self.tmpDir,self.inputDir)
            loweResolutionDEM = pcr.areatotal(pcr.cover(highResolutionDEM*self.cellArea, 0.0),\
                                              self.meteoDownscaleIds)/\
                                pcr.areatotal(pcr.cover(self.cellArea, 0.0),\
                                              self.meteoDownscaleIds)                  
            self.anomalyDEM = highResolutionDEM - loweResolutionDEM    # unit: meter  

            # temperature lapse rate (netCDF) file 
            self.temperLapseRateNC = vos.getFullPath(iniItems.meteoDownscalingOptions[\
                                        'temperLapseRateNC'],self.inputDir)                         
            self.temperatCorrelNC  = vos.getFullPath(iniItems.meteoDownscalingOptions[\
                                        'temperatCorrelNC'],self.inputDir)                    # TODO: Remove this criteria.                         

            # precipitation lapse rate (netCDF) file 
            self.precipLapseRateNC = vos.getFullPath(iniItems.meteoDownscalingOptions[\
                                        'precipLapseRateNC'],self.inputDir)
            self.precipitCorrelNC  = vos.getFullPath(iniItems.meteoDownscalingOptions[\
                                        'precipitCorrelNC'],self.inputDir)                    # TODO: Remove this criteria.                           

        else:
            logger.info("No forcing downscaling is implemented.")

        # forcing smoothing options: - THIS is still experimental. PS: MUST BE TESTED.
        self.forcingSmoothing = False
        if 'meteoDownscalingOptions' in iniItems.allSections and \
           'smoothingWindowsLength' in iniItems.meteoDownscalingOptions.keys():

            if float(iniItems.meteoDownscalingOptions['smoothingWindowsLength']) > 0.0:
                self.forcingSmoothing = True
                self.smoothingWindowsLength = vos.readPCRmapClone(\
                   iniItems.meteoDownscalingOptions['smoothingWindowsLength'],
                   self.cloneMap,self.tmpDir,self.inputDir)
                msg = "Forcing data are smoothed with 'windowaverage' using the window length:"+str(iniItems.meteoDownscalingOptions['smoothingWindowsLength'])
                logger.info(msg)   
 
    def perturb(self, name, **parameters):

        if name == "precipitation":

            # perturb the precipitation
            self.precipitation = self.precipitation * \
            pcr.min(pcr.max((1 + mapnormal() * parameters['standard_deviation']),0.01),2.0)
            #TODO: Please also make sure that precipitation >= 0
            #TODO: Add minimum and maximum 

        else:
            print("Error: only precipitation may be updated at this time")
            return -1


    def update(self,currTimeStep):
        #TODO: calculate  referencePotET
        pass

    def downscalePrecipitation(self, currTimeStep, useFactor = True, minCorrelationCriteria = 0.85):
        
        preSlope = 0.001 * vos.netcdf2PCRobjClone(\
                           self.precipLapseRateNC,"precipitation",\
                           currTimeStep.month, useDoy = "Yes",\
                           cloneMapFileName=self.cloneMap,\
                           LatitudeLongitude = True)
        preSlope = pcr.cover(preSlope, 0.0)
        preSlope = pcr.max(0.,preSlope)
        
        preCriteria = vos.netcdf2PCRobjClone(\
                     self.precipitCorrelNC,"precipitation",\
                     currTimeStep.month, useDoy = "Yes",\
                     cloneMapFileName=self.cloneMap,\
                     LatitudeLongitude = True)
        preSlope = pcr.ifthenelse(preCriteria > minCorrelationCriteria,\
                   preSlope, 0.0)             
        preSlope = pcr.cover(preSlope, 0.0)
    
        if useFactor == True:
            factor = pcr.max(0.,self.precipitation + preSlope*self.anomalyDEM)
            factor = factor / \
                     pcr.areaaverage(factor, self.meteoDownscaleIds)
            factor = pcr.cover(factor, 1.0)
            self.precipitation = factor * self.precipitation
        else:
            self.precipitation = self.precipitation + preSlope*self.anomalyDEM

        self.precipitation = pcr.max(0.0, self.precipitation)

    def downscaleTemperature(self, currTimeStep, useFactor = False, maxCorrelationCriteria = -0.75, zeroCelciusInKelvin = 273.15):
        
        tmpSlope = 1.000 * vos.netcdf2PCRobjClone(\
                           self.temperLapseRateNC,"temperature",\
                           currTimeStep.month, useDoy = "Yes",\
                           cloneMapFileName=self.cloneMap,\
                           LatitudeLongitude = True)
        tmpSlope = pcr.min(0.,tmpSlope)  # must be negative
        tmpCriteria = vos.netcdf2PCRobjClone(\
                      self.temperatCorrelNC,"temperature",\
                      currTimeStep.month, useDoy = "Yes",\
                      cloneMapFileName=self.cloneMap,\
                      LatitudeLongitude = True)
        tmpSlope = pcr.ifthenelse(tmpCriteria < maxCorrelationCriteria,\
                   tmpSlope, 0.0)             
        tmpSlope = pcr.cover(tmpSlope, 0.0)
    
        if useFactor == True:
            temperatureInKelvin = self.temperature + zeroCelciusInKelvin
            factor = pcr.max(0.0, temperatureInKelvin + tmpSlope * self.anomalyDEM)
            factor = factor / \
                     pcr.areaaverage(factor, self.meteoDownscaleIds)
            factor = pcr.cover(factor, 1.0)
            self.temperature = factor * temperatureInKelvin - zeroCelciusInKelvin
        else:
            self.temperature = self.temperature + tmpSlope*self.anomalyDEM

    def downscaleReferenceETPot(self, zeroCelciusInKelvin = 273.15):
        
        temperatureInKelvin = self.temperature + zeroCelciusInKelvin
        factor = pcr.max(0.0, temperatureInKelvin)
        factor = factor / \
                 pcr.areaaverage(factor, self.meteoDownscaleIds)
        factor = pcr.cover(factor, 1.0)
        self.referencePotET = pcr.max(0.0, factor * self.referencePotET)

    def read_forcings(self,currTimeStep):
        # reading precipitation:
        self.precipitation = vos.netcdf2PCRobjClone(\
                                  self.preFileNC,self.precipitationVarName,\
                                  str(currTimeStep.fulldate), 
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)

        self.precipitation = pcr.max(0.,self.precipitation*\
                self.precipitationCorrectionFactor)
        self.precipitation = pcr.cover( self.precipitation, 0.0)
        
        # ignore very small values of precipitation (less than 0.00001 m/day or less than 0.01 kg.m-2.day-1 )
        if self.usingDailyTimeStepForcingData:
            self.precipitation = pcr.rounddown(self.precipitation*100000.)/100000.

        # reading temperature
        self.temperature = vos.netcdf2PCRobjClone(\
                                 self.tmpFileNC,self.temperatureVarName,\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        if cellvalue(mapmaximum(self.temperature),1,1)[0] > 150.:
            self.temperature = self.temperature - 273.15


        if self.refETPotMethod == 'Hargreaves':
            self.temperatureMax = vos.netcdf2PCRobjClone(\
                                 self.tmpMaxFileNC,'tmax',\
                                 str(currTimeStep.fulldate),
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
            self.temperatureMin = vos.netcdf2PCRobjClone(\
                                 self.tmpMinFileNC,'tmin',\
                                 str(currTimeStep.fulldate),
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)

        # Downscaling precipitation and temperature
        if self.downscalePrecipitationOption: self.downscalePrecipitation(currTimeStep)
        if self.downscaleTemperatureOption: self.downscaleTemperature(currTimeStep)

        # calculate or obtain referencePotET
        if self.refETPotMethod == 'Hamon': self.referencePotET = \
                                  refPotET.HamonPotET(self.temperature,\
                                                      currTimeStep.doy,\
                                                      self.latitudes)

        if self.refETPotMethod == 'Hargreaves': self.referencePotET = \
                                  refPotET.HargreavesPotET(self.temperatureMax,\
                                                      self.temperatureMin,\
                                                      currTimeStep.doy,\
                                                      self.latitudes)

        if self.refETPotMethod == 'Input':
            self.referencePotET = vos.netcdf2PCRobjClone(\
                                     self.etpFileNC,self.evaporationVarName,\
                                     str(currTimeStep.fulldate), 
                                     useDoy = None,
                                      cloneMapFileName=self.cloneMap,\
                                      LatitudeLongitude = True)

        # Downscaling referenceETPot (based on temperature)
        if self.downscaleReferenceETPotOption: self.downscaleReferenceETPot()
        
        # smoothing:
        if self.forcingSmoothing == True:
            logger.info("Forcing data are smoothed.")   
            self.precipitation  = pcr.windowaverage(self.precipitation , self.smoothingWindowsLength)
            self.temperature    = pcr.windowaverage(self.temperature   , self.smoothingWindowsLength)
            self.referencePotET = pcr.windowaverage(self.referencePotET, self.smoothingWindowsLength)
        
        # define precipitation, temperature and referencePotET ONLY at landmask area (for reporting):
        self.precipitation  = pcr.ifthen(self.landmask, self.precipitation)
        self.temperature    = pcr.ifthen(self.landmask, self.temperature)
        self.referencePotET = pcr.ifthen(self.landmask, self.referencePotET)
 
        if self.report == True:
            timeStamp = datetime.datetime(currTimeStep.year,\
                                          currTimeStep.month,\
                                          currTimeStep.day,\
                                          0)
            # writing daily output to netcdf files
            timestepPCR = currTimeStep.timeStepPCR
            if self.outDailyTotNC[0] != "None":
                for var in self.outDailyTotNC:
                    self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_dailyTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,timestepPCR-1)

            # writing monthly output to netcdf files
            # -cummulative
            if self.outMonthTotNC[0] != "None":
                for var in self.outMonthTotNC:

                    # introduce variables at the beginning of simulation or
                    #     reset variables at the beginning of the month
                    if currTimeStep.timeStepPCR == 1 or \
                       currTimeStep.day == 1:\
                       vars(self)[var+'MonthTot'] = pcr.scalar(0.0)

                    # accumulating
                    vars(self)[var+'MonthTot'] += vars(self)[var]

                    # reporting at the end of the month:
                    if currTimeStep.endMonth == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'MonthTot'),\
                           vos.MV),timeStamp,currTimeStep.monthIdx-1)
            # -average
            if self.outMonthAvgNC[0] != "None":
                for var in self.outMonthAvgNC:
                    # only if a accumulator variable has not been defined: 
                    if var not in self.outMonthTotNC: 

                        # introduce accumulator at the beginning of simulation or
                        #     reset accumulator at the beginning of the month
                        if currTimeStep.timeStepPCR == 1 or \
                           currTimeStep.day == 1:\
                           vars(self)[var+'MonthTot'] = pcr.scalar(0.0)
                        # accumulating
                        vars(self)[var+'MonthTot'] += vars(self)[var]

                    # calculating average & reporting at the end of the month:
                    if currTimeStep.endMonth == True:
                        vars(self)[var+'MonthAvg'] = vars(self)[var+'MonthTot']/\
                                                     currTimeStep.day  
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthAvg.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'MonthAvg'),\
                           vos.MV),timeStamp,currTimeStep.monthIdx-1)
            #
            # -last day of the month
            if self.outMonthEndNC[0] != "None":
                for var in self.outMonthEndNC:
                    # reporting at the end of the month:
                    if currTimeStep.endMonth == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthEnd.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,currTimeStep.monthIdx-1)

            # writing yearly output to netcdf files
            # -cummulative
            if self.outAnnuaTotNC[0] != "None":
                for var in self.outAnnuaTotNC:

                    # introduce variables at the beginning of simulation or
                    #     reset variables at the beginning of the month
                    if currTimeStep.timeStepPCR == 1 or \
                       currTimeStep.doy == 1:\
                       vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)

                    # accumulating
                    vars(self)[var+'AnnuaTot'] += vars(self)[var]

                    # reporting at the end of the year:
                    if currTimeStep.endYear == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'AnnuaTot'),\
                           vos.MV),timeStamp,currTimeStep.annuaIdx-1)
            # -average
            if self.outAnnuaAvgNC[0] != "None":
                for var in self.outAnnuaAvgNC:
                    # only if a accumulator variable has not been defined: 
                    if var not in self.outAnnuaTotNC: 
                        # introduce accumulator at the beginning of simulation or
                        #     reset accumulator at the beginning of the year
                        if currTimeStep.timeStepPCR == 1 or \
                           currTimeStep.doy == 1:\
                           vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)
                        # accumulating
                        vars(self)[var+'AnnuaTot'] += vars(self)[var]
                    #
                    # calculating average & reporting at the end of the year:
                    if currTimeStep.endYear == True:
                        vars(self)[var+'AnnuaAvg'] = vars(self)[var+'AnnuaTot']/\
                                                     currTimeStep.doy  
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaAvg.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'AnnuaAvg'),\
                           vos.MV),timeStamp,currTimeStep.annuaIdx-1)
            #
            # -last day of the year
            if self.outAnnuaEndNC[0] != "None":
                for var in self.outAnnuaEndNC:
                    # reporting at the end of the year:
                    if currTimeStep.endYear == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaEnd.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,currTimeStep.annuaIdx-1)

    def read_forcings_forecasting(self,currTimeStep,sampleNumber):
        # reading precipitation:
        #try:
            #precipitationCorrectionFactor = float(iniItems.meteoOptions['precipitationCorrectionFactor'])
            #print "Nice"
        #except:
            #precipitationCorrectionFactor = float(1.0)
        #print iniItems.meteoOptions['precipitationCorrectionFactor']
        if self.nameModel == "FLOR" or self.nameModel == "CCSM":
            self.precipitation = vos.netcdf2PCRobjClone(\
                                  self.generate_Ensemble_Forcing_Name(self.preFileNC,sampleNumber),self.precipitationVarName,\
                                  int(currTimeStep.timeStepPCR), 
                                  useDoy = "Yes",
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        elif self.nameModel == "PGF":
            deltaYear = 0
            if currTimeStep.startTime.year != currTimeStep.currTime.year:
              deltaYear = 1
            try:
              EPSdate = datetime.datetime((self.sampleYears[sampleNumber-1]+deltaYear),currTimeStep.currTime.month,currTimeStep.currTime.day)
              EPSdate = '%d-%.2d-%.2d' %((self.sampleYears[sampleNumber-1]+deltaYear),currTimeStep.currTime.month,currTimeStep.currTime.day)
            except:
              EPSdate = datetime.datetime((self.sampleYears[sampleNumber-1]+deltaYear),currTimeStep.currTime.month,currTimeStep.currTime.day-1)
              EPSdate = '%d-%.2d-%.2d' %((self.sampleYears[sampleNumber-1]+deltaYear),currTimeStep.currTime.month,currTimeStep.currTime.day-1)
            logger.info('EPS using forcing from '+EPSdate)
            self.precipitation = vos.netcdf2PCRobjClone(\
                                  self.preFileNC,self.precipitationVarName,\
                                  EPSdate,
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        elif self.nameModel[0:6] == "revPGF":
            logger.info('revEPS using forcing from '+str(currTimeStep.fulldate))
            self.precipitation = vos.netcdf2PCRobjClone(\
                                  self.preFileNC,self.precipitationVarName,\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        elif self.nameModel == "Weighted" or self.nameModel == "WeightedEqual":
            self.precipitation = vos.netcdf2PCRobjClone(\
                                  self.preFileNC,self.precipitationVarName,\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
            self.precipitationVar = vos.netcdf2PCRobjClone(\
                                  self.preFileNC,self.precipitationVarName+str("_var"),\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
            self.randomPrec = np.minimum(np.maximum(self.randomPrec + np.random.uniform(-0.05,0.05,1),0.1),0.9)
            self.precipitation = self.precipitation + scipy.stats.norm.ppf(self.randomPrec)[0]*self.precipitationVar**0.5
        else:    
            self.precipitation = vos.netcdf2PCRobjClone(\
                                  self.generate_Ensemble_Forcing_Name(self.preFileNC,sampleNumber),self.precipitationVarName,\
                                  str(currTimeStep.fulldate), 
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)

        self.precipitation = pcr.max(0.,self.precipitation*\
                self.precipitationCorrectionFactor)
        self.precipitation = pcr.cover( self.precipitation, 0.0)
        
        # correct CDF if option is selected
        if self.matchingCDF:
            self.precipitation = self.matchCDF(self.precipitation, self.precipitationForcingCDF, self.precipitationReferenceCDF, var="prec")
            msg = "Precipiation forcing data is succesfully match to reference CDF for ensemble member "+str(sampleNumber)
            logger.info(msg)

        # ignore very small values of precipitation (less than 0.00001 m/day or less than 0.01 kg.m-2.day-1 )
        if self.usingDailyTimeStepForcingData:
            self.precipitation = pcr.rounddown(self.precipitation*100000.)/100000.

        # reading temperature
        if self.nameModel == "FLOR" or self.nameModel == "CCSM":
            self.temperature = vos.netcdf2PCRobjClone(\
                                 self.generate_Ensemble_Forcing_Name(self.tmpFileNC,sampleNumber),self.temperatureVarName,\
                                 int(currTimeStep.timeStepPCR), 
                                 useDoy = "Yes",
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        elif self.nameModel == "PGF":
            deltaYear = 0
            if currTimeStep.startTime.year != currTimeStep.currTime.year:
              deltaYear = 1
            try:
              EPSdate = datetime.datetime((self.sampleYears[sampleNumber-1]+deltaYear),currTimeStep.currTime.month,currTimeStep.currTime.day)
              EPSdate = '%d-%.2d-%.2d' %((self.sampleYears[sampleNumber-1]+deltaYear),currTimeStep.currTime.month,currTimeStep.currTime.day)
            except:
              EPSdate = datetime.datetime((self.sampleYears[sampleNumber-1]+deltaYear),currTimeStep.currTime.month,currTimeStep.currTime.day-1)
              EPSdate = '%d-%.2d-%.2d' %((self.sampleYears[sampleNumber-1]+deltaYear),currTimeStep.currTime.month,currTimeStep.currTime.day-1)
            logger.info('EPS using forcing from '+EPSdate)
            self.temperature = vos.netcdf2PCRobjClone(\
                                  self.tmpFileNC,self.temperatureVarName,\
                                  EPSdate,
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        elif self.nameModel[0:6] == "revPGF":
            logger.info('revEPS using forcing from '+str(currTimeStep.fulldate))
            self.temperature = vos.netcdf2PCRobjClone(\
                                  self.tmpFileNC,self.temperatureVarName,\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        elif self.nameModel == "Weighted" or self.nameModel == "WeightedEqual":
            self.temperature = vos.netcdf2PCRobjClone(\
                                  self.tmpFileNC,self.temperatureVarName,\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
            self.temperatureVar = vos.netcdf2PCRobjClone(\
                                  self.tmpFileNC,self.temperatureVarName+str("_var"),\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
            self.randomTemp = np.minimum(np.maximum(self.randomTemp + np.random.uniform(-0.02,0.02,1),0.01),0.99)
            self.temperature = self.temperature + scipy.stats.norm.ppf(self.randomTemp)[0]*self.temperatureVar**0.5
        else:
            self.temperature = vos.netcdf2PCRobjClone(\
                                 self.generate_Ensemble_Forcing_Name(self.tmpFileNC,sampleNumber),self.temperatureVarName,\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)

        # correct CDF if option is selected
        if self.matchingCDF:
            self.temperature = self.matchCDF(self.temperature, self.temperatureForcingCDF, self.temperatureReferenceCDF, var="tas")
            msg = "Temperature forcing data is succesfully match to reference CDF for ensemble member "+str(sampleNumber)
            logger.info(msg)
        self.orgTemperature = self.temperature
            
        # Downscaling precipitation and temperature
        self.orgPrecipitation = self.precipitation
        if self.downscalePrecipitationOption: self.downscalePrecipitation(currTimeStep)
        if self.downscaleTemperatureOption: self.downscaleTemperature(currTimeStep)

        # calculate or obtain referencePotET
        if self.refETPotMethod == 'Hamon': self.referencePotET = \
                                  refPotET.HamonPotET(self.temperature,\
                                                      currTimeStep.doy,\
                                                      self.latitudes)
        if self.refETPotMethod == 'Input': 
            self.referencePotET = vos.netcdf2PCRobjClone(\
                                     self.etpFileNC,'evapotranspiration',\
                                     str(currTimeStep.fulldate), 
                                     useDoy = None,
                                      cloneMapFileName=self.cloneMap,\
                                      LatitudeLongitude = True)

        # Downscaling referenceETPot (based on temperature)
        if self.downscaleReferenceETPotOption: self.downscaleReferenceETPot()
        
        # smoothing:
        if self.forcingSmoothing == True:
            logger.info("Forcing data are smoothed.")   
            self.precipitation  = pcr.windowaverage(self.precipitation , self.smoothingWindowsLength)
            self.temperature    = pcr.windowaverage(self.temperature   , self.smoothingWindowsLength)
            self.referencePotET = pcr.windowaverage(self.referencePotET, self.smoothingWindowsLength)
        
        # define precipitation, temperature and referencePotET ONLY at landmask area (for reporting):
        self.precipitation  = pcr.ifthen(self.landmask, self.precipitation)
        self.temperature    = pcr.ifthen(self.landmask, self.temperature)
        self.referencePotET = pcr.ifthen(self.landmask, self.referencePotET)
 
        if self.report == True:
            timeStamp = datetime.datetime(currTimeStep.year,\
                                          currTimeStep.month,\
                                          currTimeStep.day,\
                                          0)
            # writing daily output to netcdf files
            timestepPCR = currTimeStep.timeStepPCR
            if self.outDailyTotNC[0] != "None":
                for var in self.outDailyTotNC:
                    self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_dailyTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,timestepPCR-1)

            # writing monthly output to netcdf files
            # -cummulative
            if self.outMonthTotNC[0] != "None":
                for var in self.outMonthTotNC:

                    # introduce variables at the beginning of simulation or
                    #     reset variables at the beginning of the month
                    if currTimeStep.timeStepPCR == 1 or \
                       currTimeStep.day == 1:\
                       vars(self)[var+'MonthTot'] = pcr.scalar(0.0)

                    # accumulating
                    vars(self)[var+'MonthTot'] += vars(self)[var]

                    # reporting at the end of the month:
                    if currTimeStep.endMonth == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'MonthTot'),\
                           vos.MV),timeStamp,currTimeStep.monthIdx-1)
            # -average
            if self.outMonthAvgNC[0] != "None":
                for var in self.outMonthAvgNC:
                    # only if a accumulator variable has not been defined: 
                    if var not in self.outMonthTotNC: 

                        # introduce accumulator at the beginning of simulation or
                        #     reset accumulator at the beginning of the month
                        if currTimeStep.timeStepPCR == 1 or \
                           currTimeStep.day == 1:\
                           vars(self)[var+'MonthTot'] = pcr.scalar(0.0)
                        # accumulating
                        vars(self)[var+'MonthTot'] += vars(self)[var]

                    # calculating average & reporting at the end of the month:
                    if currTimeStep.endMonth == True:
                        vars(self)[var+'MonthAvg'] = vars(self)[var+'MonthTot']/\
                                                     currTimeStep.day  
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthAvg.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'MonthAvg'),\
                           vos.MV),timeStamp,currTimeStep.monthIdx-1)
            #
            # -last day of the month
            if self.outMonthEndNC[0] != "None":
                for var in self.outMonthEndNC:
                    # reporting at the end of the month:
                    if currTimeStep.endMonth == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthEnd.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,currTimeStep.monthIdx-1)

            # writing yearly output to netcdf files
            # -cummulative
            if self.outAnnuaTotNC[0] != "None":
                for var in self.outAnnuaTotNC:

                    # introduce variables at the beginning of simulation or
                    #     reset variables at the beginning of the month
                    if currTimeStep.timeStepPCR == 1 or \
                       currTimeStep.doy == 1:\
                       vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)

                    # accumulating
                    vars(self)[var+'AnnuaTot'] += vars(self)[var]

                    # reporting at the end of the year:
                    if currTimeStep.endYear == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'AnnuaTot'),\
                           vos.MV),timeStamp,currTimeStep.annuaIdx-1)
            # -average
            if self.outAnnuaAvgNC[0] != "None":
                for var in self.outAnnuaAvgNC:
                    # only if a accumulator variable has not been defined: 
                    if var not in self.outAnnuaTotNC: 
                        # introduce accumulator at the beginning of simulation or
                        #     reset accumulator at the beginning of the year
                        if currTimeStep.timeStepPCR == 1 or \
                           currTimeStep.doy == 1:\
                           vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)
                        # accumulating
                        vars(self)[var+'AnnuaTot'] += vars(self)[var]
                    #
                    # calculating average & reporting at the end of the year:
                    if currTimeStep.endYear == True:
                        vars(self)[var+'AnnuaAvg'] = vars(self)[var+'AnnuaTot']/\
                                                     currTimeStep.doy  
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaAvg.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'AnnuaAvg'),\
                           vos.MV),timeStamp,currTimeStep.annuaIdx-1)
            #
            # -last day of the year
            if self.outAnnuaEndNC[0] != "None":
                for var in self.outAnnuaEndNC:
                    # reporting at the end of the year:
                    if currTimeStep.endYear == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaEnd.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,currTimeStep.annuaIdx-1)

    def generate_Ensemble_Forcing_Name(self,orgFile,sampleNumber):
        splitName = str.split(orgFile, "_")
        splitName[4] = "r"+str(sampleNumber)+"i1p1"
        newFile = "_".join(splitName)
        return newFile


    def matchCDF(self,dataPCR, orgDataCDF, refDataCDF, var="prec"):
        nx, ny, orgMax = orgDataCDF.shape
        refMax = refDataCDF.shape[2]
        data = pcr2numpy(dataPCR,-999.)
    
        if var == "prec":
            nonZero = orgDataCDF > 0.0
            nonZeroOrg = nonZero.argmax(2)-1
            
            nonZero = refDataCDF > 0.0
            nonZeroRef = nonZero.argmax(2)-1
            
        optDif = np.zeros((nx,ny))+9e9
        lowDif = np.zeros((nx,ny))
        highDif = np.zeros((nx,ny))
        optVal = np.zeros((nx,ny))+9e10
        out = np.zeros((nx,ny))
        
        for p in range(orgMax):
            absDif = np.abs(orgDataCDF[:,:,p] - data)
            improveVal = optDif > absDif
            optVal[improveVal] = p
            optDif[improveVal] = absDif[improveVal]
            lowDif[improveVal] = np.abs(orgDataCDF[:,:,np.max([p-1, 0])] - data)[improveVal]
            highDif[improveVal] = np.abs(orgDataCDF[:,:,np.min([p+1, orgMax-1])] - data)[improveVal]
        
        upInt = highDif+1e-05 < lowDif
        lowInt = highDif > lowDif+1e-05
        out = optVal
        out[upInt] = optVal[upInt] + optDif[upInt]/(highDif[upInt] + optDif[upInt])
        out[lowInt] = optVal[lowInt] - optDif[lowInt]/np.maximum(lowDif[lowInt] + optDif[lowInt], 1e-10)
        out = out/(orgMax-1.)
        out[out < 0.0] = 0.0
        out[out > 100.0] = 100.0

        transData= np.zeros((nx,ny))
        out = out * (refMax-1.)
        for p in range(refMax):
            selVal = np.floor(out) == p
            maxVal = refDataCDF[:,:,np.min([p+1,refMax-1])][selVal]
            minVal = refDataCDF[:,:,np.max([p,0])][selVal]
            lowData = minVal + (out[selVal] - p) * (maxVal - minVal)
            highData = minVal + ((p+1) - out[selVal]) * (maxVal - minVal)
            transData[selVal] = (lowData + highData)/2
            if p == 0 and var == "prec":
                randomRainChance = np.maximum(nonZeroOrg - nonZeroRef,0.0)/np.maximum(nonZeroOrg, 1e-10)
                randomRainChance[randomRainChance >= 0.99999] = 2
                randomRain = (np.random.random((nx,ny)) - (1-randomRainChance)) / np.maximum(randomRainChance, 1e-10)
                randomRain[randomRain < 0] = 0.0
                randomRain[randomRain > 1.0] = 0.0
                rainPercentile = (nonZeroRef-0.001) + np.ceil((nonZeroOrg - nonZeroRef) * randomRain)
                out[selVal] = rainPercentile[selVal]
        
        outPCRdata = pcr.numpy2pcr(pcr.Scalar, transData, -999)
        return(outPCRdata)
      
    def update_forcings_CDF(self,currTimeStep):
        if self.matchingCDF:
            self.precipitationReferenceCDF = vos.netcdf2PCRobjCloneMultiDim(\
                              self.filePrecipitationReferenceCDF,'prec',\
                              str(currTimeStep.fulldate), 
                              useDoy = 'month',
                              cloneMapFileName=self.cloneMap,\
                              LatitudeLongitude = True)

            self.precipitationForcingCDF = vos.netcdf2PCRobjCloneMultiDim(\
                              self.filePrecipitationInputCDF,self.precipitationVarName,\
                              str(currTimeStep.fulldate), 
                              useDoy = 'month',
                              cloneMapFileName=self.cloneMap,\
                              LatitudeLongitude = True)

            self.temperatureReferenceCDF = vos.netcdf2PCRobjCloneMultiDim(\
                              self.fileTemperatureReferenceCDF,'tas',\
                              str(currTimeStep.fulldate), 
                              useDoy = 'month',
                              cloneMapFileName=self.cloneMap,\
                              LatitudeLongitude = True)

            self.temperatureForcingCDF = vos.netcdf2PCRobjCloneMultiDim(\
                              self.fileTemperatureInputCDF,self.temperatureVarName,\
                              str(currTimeStep.fulldate), 
                              useDoy = 'month',
                              cloneMapFileName=self.cloneMap,\
                              LatitudeLongitude = True)
