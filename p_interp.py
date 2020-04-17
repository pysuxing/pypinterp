import os.path
import glob
import datetime
import numpy as np
import netCDF4 as nc4

class PInterpError(Exception):
  def __init__(self, msg):
    self.msg = msg
  def __str__(self):
    return msg

class NameList:
  fields = ['PRES','TT','GHT','RH']
  met_em_fields = [
    dict(name = 'LANDUSEF', unit='category', desc = '24-category USGS landuse', order = 'XYZ'),
    dict(name = 'SOILCTOP', unit='category', desc = '16-category top-layer soil type', order = 'XYZ'),
    dict(name = 'SOIL_LAYERS', unit='cm', desc = '', order = 'XYZ'),
    dict(name = 'SOILHGT', unit='m', desc = 'Terrain field of source analysis', order = 'XY'),
    dict(name = 'PMSL', unit='Pa', desc = 'Sea-level Pressure', order = 'XY')]
  
  def __init__(self, nldict):
    self.__dict__ = nldict
    self.canonicalize()

  def canonicalize(self):
    for k, d in self.__dict__.items():
      if isinstance(d, str):
        self.__dict__[k] = d.strip()
    if self.process not in ['all', 'list']:
      raise PInterpError('Invalide setting for process')
    for f in self.fields:
      if not f in NameList.fields:
        raise PInterpError('Invalide setting for fields ')
    if self.met_em_output:  # line 324
      if self.process != 'all':
        raise PInterpError('Process must be "all" if met_em_output==True')
      if self.interp_levels[0] < 950:
        raise PInterpError('''
        ERROR:  Lowest pressure level set in interp_levels in the
        nldict is above 950 mb.  met_em needs surface data.
        Include a pressure level close to the surface: 1000 mb
        and other mandatory levels: 925, 850, 700, 500, 400, 300, 250, 200, 150, 100''') # line 334
      self.extrapolate = True
      self.split_output = True
      self.unstagger_grid = False
      surface = self.interp_levels[0]+5
      self.interp_levels = [surface] ++ self.interp_levels
    self.interp_levels = np.array(self.interp_levels)
    self.path_to_input = self.path_to_input.rstrip(os.sep)
    self.path_to_output = self.path_to_output.rstrip(os.sep)
    if self.process == 'all':
      self.fields = NameList.fields

class Interpolator:
  def __init__(self, namelist):
    self.namelist = NameList(namelist)
    self.filenames = glob.glob(os.path.join(self.namelist.path_to_input, self.namelist.input_name))
    if not self.filenames:      # line 345
      raise PInterpError('No input file specified')
    self.collect_mefields()     # add extra fields for met_em_output
    # output summary if required
    # self.pressures = self.namelist.interp_levels * 100
    for ifilename in self.filenames:
      print(ifilename)
      self.interp_file(ifilename)
    
  def collect_mefields(self):
    if self.namelist.process == 'all' and self.namelist.met_em_output:
      ncfile = nc4.Dataset(self.filenames[0])
      mes = [me for me in NameList.met_em_fields if me['name'] not in ncfile.variables]
      self.namelist.fields.extend(me['name'] for me in NameList.met_em_fields
                                  if me['name'] not in ncfile.variables)
      ncfile.close()
      
  def interp_file(self, filename):
    ifilename = os.path.basename(filename)
    # num_metgrid_levels = self.namelist.interp_levels.size
    incfile = nc4.Dataset(filename)
    dimensions, variables, gattributes = incfile.dimensions, incfile.variables, incfile.__dict__ # .copy()
    self.debug('input file has {} dimensions, {} variables, and {} global attributes'
               .format(len(dimensions), len(variables), len(gattributes)))

    diag_processed = False
    met_em_processed = False
    timevar = incfile.variables['Times']
    timestrs = np.empty(timevar.shape, dtype='U20')
    timestrs = timevar[:]
    grid_id = gattributes.get('GRID_ID')
    if not grid_id:
      grid_id = gattributes.get('grid_id')
    if (not grid_id or grid_id > 99):
      raise PInterpError('invalid grid id')

    noutput = 1
    ntimes = timevar.size
    if self.namelist.met_em_output or self.namelist.split_output:
      noutput, ntimes = ntimes, noutput
    ofilenames = []
    for o in range(noutput):
      if self.namelist.met_em_output:
        ofilename = 'met_em.d' + str(grid_id).zfill(2) + '.' + timestrs[o].tobytes().decode() + '.nc'
      else:
        if not self.namelist.output_name:
          if ifilename.startswith('wrfout'):
            ofilename = 'wrfout_d' + str(grid_id).zfill(2) + '_' + timestrs[o].tobytes().decode()
          elif ifilename.startswith('wrfinp'):
            ofilename = 'wrfinput_d' + str(grid_id).zfill(2) + '_' + timestrs[o].tobytes().decode()
          else:
            ofilename = 'p_interp_d' + str(grid_id).zfill(2) + '_' + timestrs[o].tobytes().decode()
        else:
          ofilename = self.namelist.output_name + str(grid_id).zfill(2) + '_' + timestrs[o].tobytes().decode()
      # raise if already exists
      ofilenames.append(os.path.join(self.namelist.path_to_output, ofilename))
    if self.namelist.bit64:
      fmt = 'NETCDF3_64BIT_OFFSET'
    else:
      fmt = 'NETCDF4'
    oncfiles = [nc4.Dataset(ofilename, mode='w', clobber=True, diskless=True, persist=True, format=fmt) # clober = False
                for ofilename in ofilenames]

    welen = gattributes['WEST-EAST_GRID_DIMENSION']
    snlen = gattributes['SOUTH-NORTH_GRID_DIMENSION']
    btlen = gattributes['BOTTOM-TOP_GRID_DIMENSION']
    
    map_proj = gattributes.get('MAP_PROJ')
    truelat1 = gattributes.get('TRUELAT1')
    truelat2 = gattributes.get('TRUELAT2')
    stand_lon = gattributes.get('STAND_LON')
    ogridmapping = map_proj and truelat1 and truelat2 and stand_lon

    times = self.parse_times(timestrs)
    for onc in oncfiles:
      # create dimensions in output
      self.def_output_dimensions(onc, dimensions)
      # create global attributes in output
      self.def_output_gattributes(onc, gattributes)
      # create variables
      self.def_output_variables(onc, times)
    
    ptop = incfile.variables['P_TOP'][0]
    usertop = self.namelist.interp_levels[-1]*100
    if self.namelist.extrapolate and ptop > usertop:
      print('WARNING: Highest requested pressure level is above PTOP.'
            'Use all pressure level data above {} mb with caution.'.format(ptop))

  def def_output_variables(self, onc, times):
    self.def_output_var(onc, 'Time', np.float64, 'Time', times['timestamp'],
                        long_name='Time', units='hours since 1997-01-01 00:00:00', calendar='standard')
    self.def_output_var(onc, 'DateTime', np.int32, 'Time', times['ymdh'], long_name='Date and Time')
    self.def_output_var(onc, 'year', np.int32, 'Time', times['year'], long_name='year')
    self.def_output_var(onc, 'month', np.int32, 'Time', times['month'], long_name='month')
    self.def_output_var(onc, 'day', np.int32, 'Time', times['day'], long_name='day')
    self.def_output_var(onc, 'hour', np.int32, 'Time', times['hour'], long_name='hour')
    self.def_output_var(onc, 'minute', np.int32, 'Time', times['minute'], long_name='minute')
    self.def_output_var(onc, 'second', np.int32, 'Time', times['second'], long_name='second')
    # create the pressure variable
    self.def_output_var(onc, 'pressures', np.float32, 'num_metgrid_levels', self.namelist.interp_levels,
                          long_name='Pressure Levels', standard_name='air_pressure', units='hPa', positive='down')
  def def_output_var(self, onc, name, xtype, dims, values, **attrs):
    var = onc.createVariable(name, xtype, dims)
    var.setncatts(attrs)
    
  def parse_times(self, timestrs):
    times = np.empty(timestrs.size, dtype=[('year', np.int32),
                                           ('month', np.int32),
                                           ('day', np.int32),
                                           ('hour', np.int32),
                                           ('minute', np.int32),
                                           ('second', np.int32),
                                           ('ymdh', np.int32),
                                           ('timestamp', np.float64)])
    for i, ts in enumerate(timestrs):
      year = int(ts[0:4].tobytes())
      month = int(ts[5:7].tobytes())
      day = int(ts[8:10].tobytes())
      hour = int(ts[11:13].tobytes())
      minute = int(ts[14:16].tobytes())
      second = int(ts[17:19].tobytes())
      ymdh = year * 1000000 + month * 100000 + day * 100 + hour
      timestamp = datetime.datetime(year, month, day, hour, minute, second).timestamp()/3600
      times[i] = (year, month, day, hour, minute, second, ymdh, timestamp)
    return times
  def def_output_dimensions(self, onc, dimensions):
    if 'Time' in dimensions:
      onc.createDimension('Time', 0) # line 723
    onc.createDimension('num_metgrid_levels', self.namelist.interp_levels.size)
    if self.namelist.met_em_output:
      onc.createDimension('z-dimension0024', 24)
      onc.createDimension('z-dimension0016', 16)
      onc.createDimension('z-dimension0012', 12)
      onc.createDimension('num_st_layers', dimensions['soil_layers_stag'].size)
  def def_output_gattributes(self, onc, gattributes):
    for attrname, attr in gattributes.items():
      if attrname == 'BOTTOM-TOP_PATCH':
        continue
      if attrname == 'TITLE':
        attr += ' - ON PRES LEVELS'
      elif attrname == 'BOTTOM-TOP_GRID_DIMENSION':
        attr = self.namelist.interp_levels.size
      onc.setncattr(attrname, attr)
      if self.namelist.met_em_output:
        attrdict = {'FLAG_SLP': 1,
                    'FLAG_PSFC': 1,
                    'FLAG_SOILHGT': 1,
                    'FLAG_MF_XY': 1,
                    'FLAG_METGRID': 1,
                    'FLAG_P_INTERP': 1}
        surf_physics = gattributes['SF_SURFACE_PHYSICS']
        if surf_physics in [1,2,7]:
          attrdict['FLAG_SOIL_LAYERS'] = 1
        elif surf_physics == 3:
          attrdict['FLAG_SOIL_LEVELS'] = 1
        elif surf_physics == 0:
          attrdict['FLAG_SOIL_LEVELS'] = 0
          attrdict['FLAG_SOIL_LAYERS'] = 0
        onc.setncatts(attrdict)
    
  def debug(self, msg):
    if self.namelist.debug:
      print(msg)
    
          
namelist = dict(
  path_to_input  = '../wrfout',
  input_name     = 'wrfout_d01_2017-01-02_*00',
  path_to_output = '.',
  output_name     = '',
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

if __name__ == '__main__':
  Interpolator(namelist)
