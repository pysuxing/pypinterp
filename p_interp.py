import os.path
import glob
import numpy as np
import netCDF4 as nc4

namelist = dict(
  path_to_input  = '.',
  input_name     = 'wrfout_d01_2000-01*00',
  path_to_output = '.',
  process        = 'all',
  fields         = ['PRES','TT','GHT','RH'],
  met_em_output  = False, 
  split_output   = False,
  debug          = False,
  bit64          = False, 
  interp_levels  = [1000.,987.5,975.,962.5,950.,937.5,925., 
                    900.,875.,850.,825.,800.,750.,700.,650.,  
                    600.,550.,500.,450.,400.,350.,300.,250., 
                    225.,200.,175.,150.,137.5,125.,112.5,100., 
                    87.5,75.,62.5,50.,37.5,25.,12.5],
  unstagger_grid = True,
  # extrapolate    = 0,          # set values below ground and above model top to missing (default)
    extrapolate    = 1,          # extrapolate below ground, and set above model top to model top values
  interp_method  = 1,          # linear in p interpolation (default)
  # interp_method  = 2,          # linear in log p interpolation
)

class Interpolator:
  def __init__(self, namelist):
    self.namelist = namelist
  def canonicalize_namelist(self):
    for k, d in self.namelist.items():
      if isinstance(d):
        self.namelist[k] = d.strip()
    self.namelist.path_to_input = self.namelist.path_to_input.strip().rstrip(os.sep)
    self.namelist.path_to_output = self.namelist.path_to_output.strip().rstrip(os.sep)


if __name__ == '__main__':
  pass
