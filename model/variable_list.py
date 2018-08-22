'''
List of variables.

Created on July 7, 2014

@author: Edwin H. Sutanudjaja
'''

netcdf_short_name = {}
netcdf_unit       = {}
netcdf_monthly_total_unit = {} 
netcdf_yearly_total_unit  = {}
netcdf_long_name  = {}
description       = {}
comment           = {}
latex_symbol      = {}

# precipitation
pcrglobwb_variable_name = 'precipitation'
netcdf_short_name[pcrglobwb_variable_name] = 'precipitation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# temperature
pcrglobwb_variable_name = 'temperature'
netcdf_short_name[pcrglobwb_variable_name] = 'temperature'
netcdf_unit[pcrglobwb_variable_name]       = 'degrees Celcius'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'mean_air_temperature'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# referencePotET
pcrglobwb_variable_name = 'referencePotET'
netcdf_short_name[pcrglobwb_variable_name] = 'reference_potential_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# directRunoff                    
pcrglobwb_variable_name = 'directRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'direct_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# interflowTotal                  
pcrglobwb_variable_name = 'interflowTotal'
netcdf_short_name[pcrglobwb_variable_name] = 'interflow'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# baseflow                  
pcrglobwb_variable_name = 'baseflow'
netcdf_short_name[pcrglobwb_variable_name] = 'baseflow'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# waterBodyActEvaporation
pcrglobwb_variable_name = 'waterBodyActEvaporation'
netcdf_short_name[pcrglobwb_variable_name] = 'water_body_actual_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Flux values given are over the entire cell area (not only over surface water body fraction).'
latex_symbol[pcrglobwb_variable_name]      = None

# waterBodyPotEvaporation
pcrglobwb_variable_name = 'waterBodyPotEvaporation'
netcdf_short_name[pcrglobwb_variable_name] = 'water_body_potential_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Flux values given are over the entire cell area (not only over surface water body fraction).'
latex_symbol[pcrglobwb_variable_name]      = None

# totalEvaporation
pcrglobwb_variable_name = 'totalEvaporation'
netcdf_short_name[pcrglobwb_variable_name] = 'total_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Including from water bodies.'
latex_symbol[pcrglobwb_variable_name]      = None

# runoff
pcrglobwb_variable_name = 'runoff'
netcdf_short_name[pcrglobwb_variable_name] = 'land_surface_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "direct_runoff + interflow + baseflow, but not including local runoff from water bodies."
latex_symbol[pcrglobwb_variable_name]      = None

# accuRunoff
pcrglobwb_variable_name = 'accuRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'accumulated_land_surface_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm3.s-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "direct_runoff + interflow + baseflow, but not including local runoff from water bodies."
latex_symbol[pcrglobwb_variable_name]      = None

# accuBaseflow
pcrglobwb_variable_name = 'accuBaseflow'
netcdf_short_name[pcrglobwb_variable_name] = 'accumulated_land_surface_baseflow'
netcdf_unit[pcrglobwb_variable_name]       = 'm3.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = "m3.month-1" 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = "m3.year-1"
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# discharge
pcrglobwb_variable_name = 'discharge'
netcdf_short_name[pcrglobwb_variable_name] = 'discharge'
netcdf_unit[pcrglobwb_variable_name]       = 'm3.s-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# totalRunoff
pcrglobwb_variable_name = 'totalRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'total_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Including local changes at water bodies."
latex_symbol[pcrglobwb_variable_name]      = None

# local_water_body_flux
pcrglobwb_variable_name = 'local_water_body_flux'
netcdf_short_name[pcrglobwb_variable_name] = 'local_water_body_flux'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# accuTotalRunoff
pcrglobwb_variable_name = 'accuTotalRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'accumulated_total_surface_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm3.s-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Including runoff from water bodies."
latex_symbol[pcrglobwb_variable_name]      = None

# surfaceWaterStorage
pcrglobwb_variable_name = 'surfaceWaterStorage'
netcdf_short_name[pcrglobwb_variable_name] = 'surface_water_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Negative values may be reported, due to excessive demands.'
latex_symbol[pcrglobwb_variable_name]      = None

# waterBodyStorage 
pcrglobwb_variable_name = 'waterBodyStorage'
netcdf_short_name[pcrglobwb_variable_name] = 'lake_and_reservoir_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm3'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'The values given are for every lake and reservoir ids (not per cells) and after lake/reservoir releases/outflows.'
latex_symbol[pcrglobwb_variable_name]      = None

# dynamic water fraction
pcrglobwb_variable_name = 'dynamicFracWat'
netcdf_short_name[pcrglobwb_variable_name] = 'dynamic_water_fraction'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Fraction_of_cell_flooded'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Flooded fraction is a combination of channel, waterbody and floodplain'
latex_symbol[pcrglobwb_variable_name]      = None

# Water temperature
pcrglobwb_variable_name = 'waterTemp'
netcdf_short_name[pcrglobwb_variable_name] = 'waterTemperature'
netcdf_unit[pcrglobwb_variable_name]       = 'K'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Temperature_of_surface_water'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Surface water temperature assuming fully mixed conditions'
latex_symbol[pcrglobwb_variable_name]      = None

# ice Thickness
pcrglobwb_variable_name = 'iceThickness'
netcdf_short_name[pcrglobwb_variable_name] = 'iceThickness'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Thickness_of_ice_layer'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Thickness of ice layer on channel network'
latex_symbol[pcrglobwb_variable_name]      = None

# channel storage
pcrglobwb_variable_name = 'channelStorage'
netcdf_short_name[pcrglobwb_variable_name] = 'channelStorage'
netcdf_unit[pcrglobwb_variable_name]       = 'm3'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Total_in_channel_storage'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Riverine channel storage, values include reservoirs and lakes'
latex_symbol[pcrglobwb_variable_name]      = None

# water height
pcrglobwb_variable_name = 'waterHeight'
netcdf_short_name[pcrglobwb_variable_name] = 'waterHeight'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Total_in_channel_water_height'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Water height in channel, values include reservoirs and lakes'
latex_symbol[pcrglobwb_variable_name]      = None

#~ # remove/clear pcrglobwb_variable_name 
#~ pcrglobwb_variable_name = None
#~ del pcrglobwb_variable_name
