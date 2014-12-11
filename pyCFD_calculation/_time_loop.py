"""
module for time loops
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import numpy
import ConfigParser
import os
import pyCFD_config.config as Config_
import pyCFD_output.output as output

__OUTDIR__          = Config_.__FIELDDIR__
__OUTVTKDIR__       = Config_.__OUTITERDIR__
__runSettingsFile__ = Config_.__runSettingsFile__
#__OUTDIR__        = '_FIELD_FILES/'
#__OUTVTKDIR__     = '_OUTPUT/'

class TimeLoop:
    """
    basic class for time loops
    """
    
    def __init__(self, field_list, start_time=0., stop_time=0., dt=0.):
        """
        default constructor for the TimeLoop class
    
        :param field_list:      list of volume fields
        :type field_list:       pyCFD_fields.fields.VolumeField
        :param start_time:      physical time of first time step
        :type start_time:       float        
        :param stop_time:       physical time of first time step
        :type stop_time:        float
        :param dt:              time step (fixed or list)
        :type dt:               float or list
        
        .. note::
            
            syntax for changing dt:
            
            dt  = [ ( up_to_0, dt_0 ), ( up_to_1, dt_1 ) ]
        
            e.g.:
            
            dt  = [ (  3.5, 0.005), ( 20.0, 0.01 ) ]
        
            last dt value is used if calculating further.        
        """
        # attributes
        self.fields = field_list
        """list of fields"""
        self.saveTimes = []
        """list of timesteps where fields are saved"""
        self.saveStep = []
        """uniform saving timestep"""
        self.startTime = start_time
        """first timestep"""
        self.lastTime = stop_time
        """last timestep"""
        self.time = self.startTime
        """current time"""
        self.times = []
        """all timesteps"""
        self.vtkStart = 0
        """Write vtk file vtk_start+1. Useful when restarting"""
#        self.variable_dt = []
#        """"""
        self.saveNow = False
        """if output should be save in current time step"""
        self.timeStep = 0.
        """value or list of timestep(s)"""

        self.load_run_settings()
        self.update_timestep_list()
        
    def update_timestep_list(self):
        if isinstance(self.timeStep, float) == True:
            self.dt = self.timeStep
            """current timestep"""
            self.times = numpy.arange(self.startTime+self.dt, self.stopTime+self.dt, self.dt)
        else:
            self.dt = self.find_variable_dt()
            for i_ in xrange(len(self.variable_dt)):
                if i_ == 0:
                    self.times = numpy.append(self.times[0:-2],
                                              numpy.arange(
                                                            self.startTime+self.timeStep[0][1],
                                                            self.timeStep[0][0]+self.timeStep[0][1],
                                                            self.timeStep[0][1]
                                                          )
                                             )
                                             
                else:
                    next_start = self.times[-1]
                    self.times = numpy.append(self.times,
                                              numpy.arange(
                                                           next_start+self.variable_dt[i_][1],
                                                           self.timeStep[i_][0]+self.timeStep[i_][1],
                                                           self.timeStep[i_][1]
                                                          )
                                             )
                self.times = [self.times[i] for i in range(len(self.times)) if self.times[i] <= self.stopTime] 
            
    def find_variable_dt(self):
        """
        function find list of time steps from a given lookup table
        """
        if isinstance(self.timeStep, float) == True:
            return self.timeStep
        else:
            if self.time < self.timeStep[0][0]:
                print self.timeStep
                print self.timeStep[0]
                input()
                return self.timeStep[0][1]
            elif self.time > self.timeStep[-1][0] or self.time == self.timeStep[-1][0]:
                return self.timeStep[-1][1]
            else:                
                for i_ in xrange(len(self.timeStep)):
                    if self.time < self.timeStep[i_][0] and self.time > self.timeStep[i_-1][0]:
                        return self.timeStep[i_][1]
                    elif self.time == self.timeStep[i_][0]:
                        return self.timeStep[i_+1][1]
        
    def uniform_save_times(self, save_dt, save_start = 0.):
        """
        set up save times with uniform save steps
        
        :param save_dt:     save step
        :type save_dt:      float
        :param save_start:  default: 0.0, first timestep to save
        :type save_start:   float
        """
        self.saveTimes = numpy.arange(save_start, self.lastTime+self.dt, save_dt)[1:]
        
    def save_time(self):
        """
        save timestep if given in save times
        """
        if numpy.any(abs(self.saveTimes+self.startTime-self.time) < 1e-8):
            index_ = (abs(self.saveTimes+self.startTime-self.time) < 1e-8).tolist().index(True) + 1 + self.vtk_start
            output.write_mesh_file_with_fields(self.fields, str(index_))
            for field_ in self.fields:
                numpy.save(__OUTDIR__+field_.name, field_.V)
                for patch_ in field_.patches:
                    numpy.save(__OUTDIR__+patch_.name, patch_.values)
        
        if abs(self.time - self.lastTime) < self.dt/1000.: # self.time == self.lastTime:
            if self.startTime == 0.:
                output.write_pvd_collection(self.saveTimes)        
            else:
                output.write_standalone_pvd_collection()
            
    def save_current(self, iter_):
        """
        save current timestep
        """
        output.write_mesh_file_with_fields(self.fields, str(iter_))
            
    def print_time(self):
        """
        print iteration counter, physical time and current time step to screen
        
        >>>
        ================================================
        Timestep 81: t = 3.405s | dt = 0.005
        ================================================
        """
        time_step = 0
        for step_i,time_ in enumerate(self.times):
            if numpy.allclose(time_, self.time):
                time_step = step_i + 1
        print "================================================"
        print "Timestep "+str(time_step)+": t = "+str(self.time)+"s | dt = "+str(self.dt)
        print "================================================"
        
    def load_run_settings(self):
        run_settings = ConfigParser.ConfigParser()
        run_settings.read(__runSettingsFile__)
        self.startTime = float(run_settings.get("run_settings", "start_time"))
        self.stopTime  = float(run_settings.get("run_settings", "stop_time" ))
        self.timeStep  = float(run_settings.get("run_settings", "time_step" ))
        self.saveStep  = float(run_settings.get("run_settings", "save_step" ))
        self.vtkStart  = float(run_settings.get("run_settings", "vtk_start" ))
        self.saveNow   = bool( run_settings.get("run_settings", "save_now"  ))
        self.print_run_settings()
        (mode, ino, dev, nlink, uid, gid, size, atime, mtime_, ctime_) = os.stat(__runSettingsFile__)
        Config_.__runSettingsModTime__ = mtime_
        
    def print_run_settings(self):
        print "\nRUN SETTINGS:\n"
        print self.startTime
        print self.stopTime
        print self.timeStep
        print self.saveStep
        print self.vtkStart
        print self.saveNow
        print "\n"
        
    def load_run_settings_if_changed(self):
        (mode, ino, dev, nlink, uid, gid, size, atime, mtime_new, ctime_) = os.stat(__runSettingsFile__)
        if mtime_new > Config_.__runSettingsModTime__:
            run_settings = ConfigParser.ConfigParser()
            run_settings.read(__runSettingsFile__)
            self.startTime = float(run_settings.get("run_settings", "start_time"))
            self.stopTime  = float(run_settings.get("run_settings", "stop_time" ))
            self.timeStep  = float(run_settings.get("run_settings", "time_step" ))
            self.saveStep  = float(run_settings.get("run_settings", "save_step" ))
            self.vtkStart  = float(run_settings.get("run_settings", "vtk_start" ))
            self.saveNow   = bool( run_settings.get("run_settings", "save_now"  ))
            self.print_run_settings()
            Config_.__runSettingsModTime__ = mtime_new
            