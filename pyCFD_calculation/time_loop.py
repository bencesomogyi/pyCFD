"""
module for time loops
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import pyCFD_config.config as Config_
import numpy
import pyCFD_output.output as output


__OUTDIR__    = Config_.__FIELDDIR__
__OUTVTKDIR__ = Config_.__OUTITERDIR__
#__OUTDIR__        = '_FIELD_FILES/'
#__OUTVTKDIR__     = '_OUTPUT/'

class TimeLoop:
    """
    basic class for time loops
    """
    
    def __init__(self, field_list, start_time, stop_time, dt):
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
        self.startTime = start_time
        """first timestep"""
        self.lastTime = stop_time
        """last timestep"""
        self.time = self.startTime
        """current time"""
        self.times = []
        """all timesteps"""
        self.vtk_start = 0
        """Write vtk file vtk_start+1. Useful when restarting"""
        self.variable_dt = []
        
        if isinstance(dt, float) == True:
            self.dt = dt
            """current timestep"""
            self.times = numpy.arange(start_time+dt, stop_time+dt, dt)
        else:
            self.variable_dt = dt
            self.dt = self.find_variable_dt()
            for i_ in xrange(len(self.variable_dt)):
                if i_ == 0:
                    self.times = numpy.append(self.times[0:-2],
                                              numpy.arange(
                                                            self.startTime+self.variable_dt[0][1],
                                                            self.variable_dt[0][0]+self.variable_dt[0][1],
                                                            self.variable_dt[0][1]
                                                          )
                                             )
                                             
                else:
                    next_start = self.times[-1]
                    self.times = numpy.append(self.times,
                                              numpy.arange(
                                                           next_start+self.variable_dt[i_][1],
                                                           self.variable_dt[i_][0]+self.variable_dt[i_][1],
                                                           self.variable_dt[i_][1]
                                                          )
                                             )
                self.times = [self.times[i] for i in range(len(self.times)) if self.times[i] <= stop_time]
            
    def find_variable_dt(self):
        """
        function find list of time steps from a given lookup table
        """
        if self.variable_dt == []:
            return self.dt
        else:
            if   self.time < self.variable_dt[0][0]:
                return self.variable_dt[0][1]
            elif self.time > self.variable_dt[-1][0] or self.time == self.variable_dt[-1][0]:
                return self.variable_dt[-1][1]
            else:                
                for i_ in xrange(len(self.variable_dt)):
                    if self.time < self.variable_dt[i_][0] and self.time > self.variable_dt[i_-1][0]:
                        return self.variable_dt[i_][1]
                    elif self.time == self.variable_dt[i_][0]:
                        return self.variable_dt[i_+1][1]
        
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
        for field_ in self.fields:
                numpy.save(__OUTDIR__+field_.name, field_.V)
                for patch_ in field_.patches:
                    numpy.save(__OUTDIR__+patch_.name, patch_.values)
            
    def save_current_fields_as(self, iter_, dir_name):
        """
        save current fields into directory
        """
        for field_ in self.fields:
                numpy.save(__OUTDIR__+dir_name+"/"+field_.name, field_.V)
                for patch_ in field_.patches:
                    numpy.save(__OUTDIR__+dir_name+"/"+patch_.name, patch_.values)
            
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