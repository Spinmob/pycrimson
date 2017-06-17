from gnuradio import gr  as _gr
from gnuradio import uhd as _uhd

import thread      as _thread
import time        as _time
import numpy       as _n

import spinmob     as _s
import spinmob.egg as _egg
import pycrimson   as _pc


class data_streamer():
    """
    Simple GUI for streaming / analyzing / saving data from the Crimson.    
    """

    def __init__(self, rx_channels=[0,1,2,3], buffer_size=500):
        
        # Keep the list of available channels.
        self._rx_channels = rx_channels

        # Build the GUI first, connecting all the signals and loading the 
        # previous settings.
        self._build_gui()
        
        # Create the pycrimson scripted interface (and connect!)
        self.crimson = _pc.crimson(rx_channels, buffer_size)
        
        # Update settings from the Crimson itself
        self._get_settings_from_crimson()      
        
        # Set the buffer size to the previous value
        self.crimson.buffer.set_size(self.settings['Software/Collect/buffer_size'])        
        
        # Show it!
        self.show()
        
        # Start the crimson.
        self.crimson.start()
        
        # Start the collection timer
        self.timer_collect.start()
    
    def _build_gui(self):
        """
        Places all the controls and plots etc, loads previous config.
        """        
        # Embeddable docker to hold GUI stuff
        self.docker = _egg.gui.Docker('Crimson Data Streamer', autosettings_path='data_streamer_window')
        
        # Top controls
        self.g_top = self.docker.place_object(_egg.gui.GridLayout(False))
        self.b_reset         = self.g_top.place_object(_egg.gui.Button("Reset"))
        self.n_trace_counter = self.g_top.place_object(_egg.gui.NumberBox(0))
        self.n_psd_counter   = self.g_top.place_object(_egg.gui.NumberBox(0))
        
        # button functions
        self.b_reset.signal_clicked.connect(self._b_reset_clicked)
        
        # Data collection timer
        self.timer_collect = _egg.pyqtgraph.QtCore.QTimer()
        self.timer_collect.setInterval(1)
        self.timer_collect.timeout.connect(self._timer_collect_tick)
        
        # Create Settings tree and parameters
        self.docker.new_autorow()
        self.g_stuff = self.docker.place_object(_egg.gui.GridLayout(False), alignment=0)
        self.settings = self.g_stuff.place_object(_egg.gui.TreeDictionary('data_streamer_settings')).set_width(300)       

        self.settings.add_parameter('Crimson/sample_rate',      1e5,  type='float', limits=(1e4,325e6), decimals=12, suffix='Hz', siPrefix=True, dec=True)
        self.settings.add_parameter('Crimson/settle_time',  1e-3, type='float', limits=(0,None), dec=True, suffix='s', siPrefix=True)      

        self.settings.add_parameter('Crimson/RXA', False, type='bool')
        self.settings.add_parameter('Crimson/RXA/center_frequency', 7e6,  type='float', limits=(0,7e9),     decimals=12, suffix='Hz', siPrefix=True, step=0.1e6)
        self.settings.add_parameter('Crimson/RXA/gain',             0,    type='float', limits=(0,31.5), step=0.5,    suffix=' dB')

        self.settings.add_parameter('Crimson/RXB', False, type='bool')
        self.settings.add_parameter('Crimson/RXB/center_frequency', 7e6,  type='float', limits=(0,7e9),     decimals=12, suffix='Hz', siPrefix=True, step=0.1e6)
        self.settings.add_parameter('Crimson/RXB/gain',             0,    type='float', limits=(0,31.5), step=0.5,    suffix=' dB')

        self.settings.add_parameter('Crimson/RXC', False, type='bool')
        self.settings.add_parameter('Crimson/RXC/center_frequency', 7e6,  type='float', limits=(0,7e9),     decimals=12, suffix='Hz', siPrefix=True, step=0.1e6)
        self.settings.add_parameter('Crimson/RXC/gain',             0,    type='float', limits=(0,31.5), step=0.5,    suffix=' dB')

        self.settings.add_parameter('Crimson/RXD', False, type='bool')
        self.settings.add_parameter('Crimson/RXD/center_frequency', 7e6,  type='float', limits=(0,7e9),     decimals=12, suffix='Hz', siPrefix=True, step=0.1e6)
        self.settings.add_parameter('Crimson/RXD/gain',             0,    type='float', limits=(0,31.5), step=0.5,    suffix=' dB')

        self.settings.add_parameter('Software/Collect',                False, type='bool')        
        self.settings.add_parameter('Software/Collect/buffer_size',      1e3, type='int',   limits=(1,None), suffix=' packets', dec=True)
        self.settings.add_parameter('Software/Collect/calibration',2.053e-05, type='float')
        self.settings.add_parameter('Software/Collect/time',            0.01, type='float', limits=(1e-6,None), suffix='s', siPrefix=True, dec=True)
        self.settings.add_parameter('Software/Collect/keep_extra',     False, type='bool')        
        self.settings.add_parameter('Software/Collect/trigger_channel', None, type='list', values=[None]+self.get_channel_names())        
        
        self.settings.add_parameter('Software/PSD',            False, type='bool')        
        self.settings.add_parameter('Software/PSD/averages',       0, type='float')
        self.settings.add_parameter('Software/PSD/window',     "none",type='str')
        self.settings.add_parameter('Software/PSD/pow2',       True,  type='bool')
        

        # Load and set previous settings (if any)        
        self.settings.load()
        
        # Set it up so that any changes are sent to the Crimson
        self.settings.connect_signal_changed('Crimson/sample_rate',      self._settings_sample_rate_changed)
        self.settings.connect_signal_changed('Crimson/RXA/center_frequency', self._settings_center_frequency_changed)
        self.settings.connect_signal_changed('Crimson/RXA/gain',             self._settings_gain_changed)
        
        self.settings.connect_signal_changed('Software/Collect',                 self._settings_collect_changed)        
        self.settings.connect_signal_changed('Software/Collect/buffer_size',     self._settings_buffer_size_changed)
        
        # for saving
        self.settings.connect_any_signal_changed(self._any_setting_changed)

        # Plotter tabs
        self.tabs_plots = self.g_stuff.place_object(_egg.gui.TabArea(False, 'data_streamer_tabs_plots'), alignment=0)
        
        # Raw time trace    
        self.tab_raw = self.tabs_plots.add_tab('Raw')                
        self.p_raw = self.tab_raw.place_object(_egg.gui.DataboxPlot("*.txt", "data_streamer_p_raw"), alignment=0)
        self.p_raw_trigger_level = _egg.pyqtgraph.InfiniteLine(0, movable=True, angle=0)
        self.p_raw.ROIs.append([self.p_raw_trigger_level])
        self.p_raw.load_gui_settings()
        
        # "Calibrated" time trace
        self.tab_volts = self.tabs_plots.add_tab('Volts')
        self.p_volts   = self.tab_volts.place_object(_egg.gui.DataboxPlot("*.txt", "data_streamer_p_volts"), alignment=0)
        
        # PSD
        self.tab_psd = self.tabs_plots.add_tab('PSD')        
        self.p_psd   = self.tab_psd.place_object(_egg.gui.DataboxPlot("*.txt", "data_streamer_p_psd"), alignment=0)
        
        # Info grid
        self.pb_buffer       = self.g_top.place_object(_egg.pyqtgraph.QtGui.QProgressBar())        
        self.b_overrun       = self.g_top.place_object(_egg.gui.Button('Buffer Overrun', True))
        
        # load last selected tab
        self.tabs_plots.load_gui_settings()
        
        # Disable and hide the unused channels.
        for n in range(4):
            if not n in self._rx_channels:
                self.settings.block_events()
                name = self.get_channel_names()[n]
                self.settings.hide_parameter(name)
                self.settings.unblock_events()

    def _any_setting_changed(self, *a): 
        
        # Reset the PSD averages and save the settings
        self.reset_acquisition()        
        self.settings.save()

    def _settings_collect_changed(self, *a):
        
        # If we're taking new data, reset the PSD, flush the buffer, reset the overrun        
        if self.settings['Software/Collect']: self.reset_acquisition()
            
    def _b_reset_clicked(self, *a): 
        self.reset_acquisition()

    def _get_settings_from_crimson(self):
        """
        Asks the crimson what all the settings currently are, and updates
        the GUI. This also blocks all user-defined signals for changing the GUI.
        """
        self.settings.set_value('Crimson/sample_rate',      self.crimson.get_sample_rate(),      True)
        self.settings.set_value('Crimson/RXA/center_frequency', self.crimson.get_center_frequency(), True)
        self.settings.set_value('Crimson/RXA/gain',             self.crimson.get_gain(),             True)
        

    def _settings_sample_rate_changed(self, *a): 
        self.timer_collect.stop()
        
        # Tell the Crimson then get the actual value
        self.crimson.set_sample_rate(self.settings['Crimson/sample_rate'])
        self.settings.set_value('Crimson/sample_rate', self.crimson.get_sample_rate(), True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/settle_time'])
        self.timer_collect.start()
        
    def _settings_center_frequency_changed(self, *a):
        self.timer_collect.stop()
        
        # Tell the Crimson then get the actual value
        self.crimson.set_center_freq(self.settings['Crimson/RXA/center_frequency'])
        self.settings.set_value('Crimson/RXA/center_frequency', self.crimson.get_center_freq(), True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/settle_time'])
        self.timer_collect.start()
            
    def _settings_gain_changed(self, *a):
        self.timer_collect.stop()
        
        # Tell the Crimson then get the actual value
        self.crimson.set_gain(self.settings['Crimson/RXA/gain'])
        self.settings.set_value('Crimson/RXA/gain', self.crimson.get_gain(), True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/settle_time'])        
        self.timer_collect.start()

    def _settings_buffer_size_changed(self, *a):
        self.crimson.buffer.set_size(self.settings['Software/Collect/buffer_size'])
    
    def start(self): return self._top_block.start()
    
    def get_channel_names(self):
        """
        Returns the channel names.
        """
        return ['Crimson/RXA','Crimson/RXB','Crimson/RXC','Crimson/RXD']
    
    def get_enabled_channels(self):
        """
        Returns the indices of the enabled channels.
        """
        names = self.get_channel_names()
        a = []
        for n in range(4):
            if self.settings[names[n]]: a.append(n)
        return a    
    
    def _timer_collect_tick(self, *a):
        """
        Called every time the data collection timer ticks.
        """
        # Update the number of packets
        #self.pb_buffer.setValue(int(_n.round(100.0*len(self.crimson.buffer)/self.crimson.buffer.get_size())))

        # Update the user if the buffer overran
        if self.crimson.buffer.get_overruns(): self.b_overrun.set_checked(True)
                            
        if self.settings['Software/Collect']: 

            # Get the time step
            dt = 1.0/self.settings['Crimson/sample_rate']

            # Total time and number of samples
            T = self.settings['Software/Collect/time']                        
            N = int(_n.round(T/dt))
            
            # Get all the packets we need            
            samples = self.crimson.get_data(N, keep_all=self.settings['Software/Collect/keep_extra'])
            
            # If we timed out, the last element is False.
            if samples[-1]: return

            # Success! Increment the counter
            self.n_trace_counter.increment()
            
            # Store the data
            
            # Remove old data and put a placeholder in for time
            self.p_raw.clear()
            self.p_raw['t'] = []
            
            # Store the enabled channels
            for n in self.get_enabled_channels():
                
                # By name                
                self.p_raw[self.get_channel_names()[n]+"/X"] = samples[n][0]
                self.p_raw[self.get_channel_names()[n]+"/Y"] = samples[n][1]
                
                # Check for overruns
                if samples[n][2]: self.b_overrun.set_checked()   
            
            # If we're not supposed to trigger, use time from 0 to max            
            name = self.settings['Software/Collect/trigger_channel']
            if name == None:
                self.p_raw['t'] = _n.arange(0, (len(samples[0][0])-0.5)*dt, dt)
            
            # Otherwise, find the first index at which
            # it crosses the threshold
            else:
                # Choose the channel, and use the x quadrature for now
                x = self.p_raw[name+"/X"]                
                
                # Trigger level
                l = self.p_raw_trigger_level.value()

                # All values where we're below the trigger level
                ns0 = _n.where(x<l)[0]            
                
                # If we have any values below the trigger and we're supposed to 
                # trigger, find the crossing
                if len(ns0) and self.b_trigger.is_checked():
                    
                    # first value below trigger
                    n0  = ns0[0]
    
                    # indices after n0 where we're above the trigger
                    ns1 = _n.where(x[n0:]>l)[0]+n0
    
                    # If we have a crossing                
                    if(len(ns1)):  self.p_raw['t'] = _n.arange(-ns1[0]*dt,(-ns1[0]+len(x)-0.5)*dt, dt) 
                    else:          self.p_raw['t'] = _n.arange(0, (len(samples[0][0])-0.5)*dt, dt)   
                

            # Plot Raw
            self.settings.send_to_databox_header(self.p_raw)            
            self.p_raw.plot()

#            # Plot Volts
#            self.p_volts['t'] = t
#            Vx = x*self.settings['Software/Collect/calibration']*10**(-0.05*self.settings['Crimson/RXA/gain']) 
#            Vy = y*self.settings['Software/Collect/calibration']*10**(-0.05*self.settings['Crimson/RXA/gain'])
#            self.p_volts['Vx'] = Vx
#            self.p_volts['Vy'] = Vy
#            self.settings.send_to_databox_header(self.p_volts)            
#            self.p_volts.plot()
#            
#            # PSD
#            if self.settings['Software/PSD']:
#                
#                # If we want the FFT to be extra efficient
#                if self.settings['Software/PSD/pow2']: N = 2**int(_n.log2(N))              
#                
#                # Loop over the data to get as many PSD's as possible from it                
#                for k in range(int(len(t)/int(N))):
#                    
#                    # Do the PSD
#                    result = _s.fun.fft(t[k*N:(k+1)*N], Vx[k*N:(k+1)*N]+1j*Vy[k*N:(k+1)*N], 
#                                          pow2=False, window=self.settings['Software/PSD/window'])
#                    
#                    # Happens if there is an invalid window
#                    if result == None: return                
#
#                    # Convert to PSD
#                    f, fft = result
#                    df = f[1]-f[0]
#                    f = f+self.settings['Crimson/RXA/center_frequency']                
#                    P = 0.5*abs(fft)**2 / df
#                    
#                    # Reset the counter if the length or the length has changed
#                    if len(self.p_psd) and not len(self.p_psd[0]) == len(f): self.reset_acquisition()
#    
#                    # Make sure we have the data arrays
#                    if len(self.p_psd) == 0: 
#                        self.p_psd['f'] = f
#                        self.p_psd['P'] = P
#    
#                    # Do the average
#                    n = self.n_psd_counter.get_value()
#                    self.p_psd['f'] = f
#                    self.p_psd['P'] = (self.p_psd['P']*n + P)/(n+1)
#                    self.settings.send_to_databox_header(self.p_psd)
#                    self.p_psd.plot()
#                
#                    # increment the psd counter
#                    if self.settings['Software/PSD/averages'] == 0 or n+1 <= self.settings['Software/PSD/averages']:
#                        self.n_psd_counter.increment()
#                    
#                    # Otherwise we're done!                    
#                    else:
#                        self.b_collect.set_checked(False)
                    
            self.docker.process_events()
            
            
    def hide(self):  return self.docker.hide()
    def show(self):  return self.docker.show()
    
    def reset_acquisition(self):
        
        # Resets counters, clears plots, clears buffers        
        self.n_psd_counter.set_value(0)
        self.n_trace_counter.set_value(0)
        self.b_overrun.set_checked(False)
        self.crimson.buffer.flush_buffer()
        self.p_psd.clear()
        self.p_raw.clear()
        self.p_volts.clear()



if __name__ == '__main__':   
    self = data_streamer([0,3])
    