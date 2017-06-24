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
    "Simple" GUI for streaming / analyzing / saving data from the Crimson.    
    """

    def __init__(self, rx_channels=[0,1,2,3], buffer_size=500, simulation_mode=False):
        
        # Note if we're in simulation mode.
        self._simulation_mode = simulation_mode
        
        # Keep the list of available channels.
        self._rx_channels = rx_channels

        # Build the GUI first, connecting all the signals and loading the 
        # previous settings.
        self._build_gui()
        
        # Create the pycrimson scripted interface (and connect!)
        self.crimson = _pc.crimson(rx_channels, buffer_size, simulation_mode)
        
        # Update settings from the Crimson itself
        self._get_settings_from_crimson()      
        
        # Set the buffer size to the previous value
        if not self._simulation_mode:
            self.crimson.buffer.set_size(self.settings['Software/Acquisitions/buffer_size'])        
        
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
        self.b_collect       = self.g_top.place_object(_egg.gui.Button("Go!", True, False))
        self.b_reset         = self.g_top.place_object(_egg.gui.Button("Reset"))
        
        # button functions
        self.b_collect.signal_clicked.connect(self._b_collect_clicked)
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

        #self.b_collect = self.settings.add_button('Software/Acquisitions', True, False)
        #self.settings.add_parameter('Software/Acquisitions',                False, type='bool')        
        self.settings.add_parameter('Software/Acquisitions',                    0, type='int')
        self.settings.add_parameter('Software/Acquisitions/count',              0, type='int')
        self.settings.add_parameter('Software/Acquisitions/buffer_size',      1e3, type='int',   limits=(1,None), suffix=' packets', dec=True)
        self.settings.add_parameter('Software/Acquisitions/calibration',2.053e-05, type='float')
        self.settings.add_parameter('Software/Acquisitions/time',            0.01, type='float', limits=(1e-6,None), suffix='s', siPrefix=True, dec=True)
        self.settings.add_parameter('Software/Acquisitions/keep_extra',     False, type='bool')        
        self.settings.add_parameter('Software/Acquisitions/trigger_channel', None, type='list', values=['None']+self.get_channel_names())        
        
        self.settings.add_parameter('Software/PSD',            False, type='bool')        
        self.settings.add_parameter('Software/PSD/count',   0,     type='int')
        self.settings.add_parameter('Software/PSD/window',     "None",type='str')
        self.settings.add_parameter('Software/PSD/pow2',       True,  type='bool')
        
        # Load and set previous settings (if any)        
        self.settings.load()
        self.settings['Software/PSD/count'] = 0
        
        # Set it up so that any changes are sent to the Crimson
        self.settings.connect_signal_changed('Crimson/sample_rate',      self._settings_sample_rate_changed)
        self.settings.connect_signal_changed('Crimson/RXA/center_frequency', self._settings_center_frequency_changed)
        self.settings.connect_signal_changed('Crimson/RXA/gain',             self._settings_gain_changed)
        
        self.settings.connect_signal_changed('Software/Acquisitions/buffer_size',     self._settings_buffer_size_changed)
        
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
                
                # Get the name, e.g. "RXA"                
                name = self.get_channel_names()[n]
                
                # Uncheck it.                
                self.settings[name] = False
                
                # Hide the GUI for this channel.                
                self.settings.hide_parameter(name)
                
                self.settings.unblock_events()
            

    def _any_setting_changed(self, *a): 
        
        # Reset the PSD averages and save the settings
        self.settings.save()

    def _b_collect_clicked(self, *a):
        
        # If we're taking new data, reset the PSD, flush the buffer, reset the overrun        
        if self.b_collect.get_value(): self.reset_acquisition()
            
    def _b_reset_clicked(self, *a): 
        self.reset_acquisition()

    def _get_settings_from_crimson(self):
        """
        Asks the crimson what all the settings currently are, and updates
        the GUI. This also blocks all user-defined signals for changing the GUI.
        """
        self.settings.set_value('Crimson/sample_rate',          self.crimson.get_sample_rate(),      True)
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
        self.crimson.set_center_frequency(self.settings['Crimson/RXA/center_frequency'])
        self.settings.set_value('Crimson/RXA/center_frequency', self.crimson.get_center_frequency(), True)
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
        self.crimson.buffer.set_size(self.settings['Software/Acquisitions/buffer_size'])
    
    def start(self): return self._top_block.start()
    
    def get_channel_names(self):
        """
        Returns the channel names.
        """
        return ['Crimson/RXA','Crimson/RXB','Crimson/RXC','Crimson/RXD']
    
    def get_gui_enabled_channels(self):
        """
        Returns the indices of the enabled channels (as per the check boxes).

        See self.get_crimson_enabled_channels() for the list of channels that are
        always streaming.
        """
        names = self.get_channel_names()
        a = []
        for n in range(4):
            if self.settings[names[n]]: a.append(n)
        return a    

    def get_crimson_enabled_channels(self):
        """
        Returns the list of channels that can be enabled, as per the 
        initialization of this object.
        """
        return list(self._rx_channels)


    def get_and_plot_data(self):
        """
        Empties the buffer, plots and calibrates the data according to the 
        settings panel. Stores and plots the raw and calibrated data in
           
           self.p_raw
           self.p_volts

        respectively. If PSD is enabled, performs the PSD and plots it in
        
           self.p_psd.
           
        If you don't want this to plot, you can disable the plots in the gui.
        """
        
        # Get the time step
        dt = 1.0/self.settings['Crimson/sample_rate']

        # Total time and number of samples
        T = self.settings['Software/Acquisitions/time']                        
        N = int(_n.round(T/dt))
        
        # Get all the packets we need          
        samples = self.crimson.get_data(N, keep_all=self.settings['Software/Acquisitions/keep_extra'])
        
        # If we timed out, the last element is False.
        if samples[-1]: return

        # Remove old data
        self.p_raw.clear()
        self.p_volts.clear()
        
        # Default for the time column
        t = _n.arange(0, (len(samples[0][0])-0.5)*dt, dt)
        self.p_raw['t']   = t
        self.p_volts['t'] = t
        
        # Store the enabled channels
        enabled   = self.get_gui_enabled_channels()
        available = self.get_crimson_enabled_channels()         
        for n in enabled: 
            
            # Get the column labels
            nameX = self.get_channel_names()[n]+"/X"            
            nameY = self.get_channel_names()[n]+"/Y"            
            
            # Store the raw data by name             
            self.p_raw[nameX] = samples[available.index(n)][0]
            self.p_raw[nameY] = samples[available.index(n)][1]
            
            # Calibrate to volts
            Vx = self.p_volts[nameX] = self.p_raw[nameX]*self.settings['Software/Acquisitions/calibration']*10**(-0.05*self.settings['Crimson/RXA/gain']) 
            Vy = self.p_volts[nameY] = self.p_raw[nameY]*self.settings['Software/Acquisitions/calibration']*10**(-0.05*self.settings['Crimson/RXA/gain'])  
            
            # Check for overruns
            if samples[available.index(n)][2]: self.b_overrun.set_checked()   
        
            # PSD
            if self.settings['Software/PSD']:
                
                # If we want the FFT to be extra efficient, set the number of
                # points to a power of 2.
                if self.settings['Software/PSD/pow2']: N = 2**int(_n.log2(N))              
                
                # Loop over the data to get as many PSD averages as possible from it                
                for k in range(int(len(t)/int(N))):
                    
                    # Do the PSD
                    result = _s.fun.fft(t[0:N], Vx[k*N:(k+1)*N]+1j*Vy[k*N:(k+1)*N], 
                                        pow2=False, window=self.settings['Software/PSD/window'])
                    
                    # None comes back if there is an invalid window, for example
                    if not result == None:                 
        
                        # Convert to PSD
                        f, fft = result
                        df = f[1]-f[0]
                        f = f+self.settings['Crimson/RXA/center_frequency']                
                        P = 0.5*abs(fft)**2 / df
                        
                        # Reset the counter if the there exists data that has
                        # the wrong format. This should never happen as the
                        # Acquisition is reset every time a setting changes.
                        if len(self.p_psd) and (not len(self.p_psd[0]) == len(f) or not (self.p_psd[0][0]==f[0] and self.p_psd[0][-1]==f[-1])): 
                            print("WARNING: I noticed a change in the frequency domain. This shouldn't happen at this point in the code. Acquisition should be reset elsewhere.")                            
                            self.p_psd.clear()
                            self.settings['Software/PSD/count'] = 0                            
                            
                        # Get the name of the column
                        nameP = self.get_channel_names()[n]+"/P"          
        
                        # Make sure we have the data arrays (e.g. after a reset)
                        if not nameP in self.p_psd.ckeys: 
                            self.p_psd['f']   = f
                            self.p_psd[nameP] = P
        
                        # Do the average
                        n_psd = self.settings['Software/PSD/count']
                        self.p_psd[nameP] = (self.p_psd[nameP]*n_psd + P)/(n_psd+1)

                        # increment the psd counter
                        self.settings['Software/PSD/count'] += 1
        
        # If we're not supposed to trigger, use time from 0 to max            
        name = self.settings['Software/Acquisitions/trigger_channel']
            
        # Otherwise, find the first index at which
        # it crosses the threshold
        if not name == 'None':
            
            # Choose the raw channel, and use the x quadrature for now
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

                # If we have a crossing, subract the time offset
                if(len(ns1)):  
                    self.p_raw['t']   = self.p_raw['t']   - ns1[0]*dt 
                    self.p_volts['t'] = self.p_volts['t'] - ns1[0]*dt

        # Plot Raw
        self.settings.send_to_databox_header(self.p_raw)            
        self.p_raw.plot()

        # Plot Volts
        self.settings.send_to_databox_header(self.p_volts)            
        self.p_volts.plot()

        # Plot PSD
        self.settings.send_to_databox_header(self.p_psd)
        self.p_psd.plot()

        # Find out if we're done
        self.settings['Software/Acquisitions/count'] += 1
        if self.settings['Software/Acquisitions'] > 0 and \
           self.settings['Software/Acquisitions/count'] >= self.settings['Software/Acquisitions']:
               self.b_collect.set_checked(False)
        
        # Update the gui if it's not already
        self.docker.process_events()

    
    def _timer_collect_tick(self, *a):
        """
        Called every time the data collection timer ticks.
        """
        # Update the number of packets
        if not self._simulation_mode:
            self.pb_buffer.setValue(int(_n.round(100.0*len(self.crimson.buffer)/self.crimson.buffer.get_size()[0])))
        
            # Update the user if the buffer overran
            if sum(self.crimson.buffer.get_overruns()): self.b_overrun.set_checked(True)
               
        # We need to enable at least one channel and have "Collect" checked
        # to collect / plot any data
        enabled = False
        for n in self.get_gui_enabled_channels(): 
            if self.settings[self.get_channel_names()[n]]: enabled = True
                
        if self.b_collect.is_checked() and enabled: self.get_and_plot_data()



                    
            
            
    def hide(self):  return self.docker.hide()
    def show(self):  return self.docker.show()
    
    def reset_acquisition(self):
        
        # Resets counters, clears plots, clears buffers        
        self.settings['Software/PSD/count'] = 0
        self.settings['Software/Acquisitions/count'] = 0
        self.b_overrun.set_checked(False)
        if not self._simulation_mode: self.crimson.buffer.flush_buffer()
        self.p_psd.clear()
        self.p_raw.clear()
        self.p_volts.clear()



if __name__ == '__main__':   
    self = data_streamer([0], simulation_mode=True)
    