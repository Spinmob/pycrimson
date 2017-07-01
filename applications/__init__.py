#!/usr/bin/env python

from gnuradio import gr  as _gr
from gnuradio import uhd as _uhd

import thread      as _thread
import time        as _time
import numpy       as _n

import spinmob     as _s
import spinmob.egg as _egg
import pycrimson   as _pc



class rx_streamer():
    """
    GUI for streaming / analyzing / saving data from the Crimson. Requires
    spinmob and pyqtgraph python libraries.

    Parameters
    ----------
    rx_channels
        list of physically enabled channels
    buffer_size
        how many packets the buffer should hold (default 346 samples/packet)
    simulation_mode
        whether to run this without connecting to the crimson (janky right now)
    block_command_line
        whether to block the terminal; set to False if running within a python
        console if you'd like to interact with it.
    """

    def __init__(self, rx_channels=[0,2,3], buffer_size=500, simulation_mode=False, block_command_line=False):
        
        # Note if we're in simulation mode.
        self._simulation_mode = simulation_mode
        
        # Keep the list of available channels.
        self._rx_channels   = rx_channels
        self._channel_names = ['RXA','RXB','RXC','RXD']

        # Build the GUI, connecting all the signals and loading the previous settings.
        self._build_gui()

        # Create the pycrimson scripted interface (and connect!)
        self.crimson = _pc.crimson(rx_channels, buffer_size, simulation_mode, self.docker.process_events)
        
        # Update settings from the Crimson itself
        self._get_settings_from_crimson()      
        
        # Set the buffer size to the previous value
        if not self._simulation_mode:
            self.crimson.buffer.set_size(self.settings['Software/Acquisitions/buffer_size'])        
        
        # Start the crimson.
        self.crimson.start()
        
        # Start the collection timer
        self.timer_collect.start()

        # Show it!
        self.show(block_command_line)        

    
    def _build_gui(self):
        """
        Places all the controls and plots etc, loads previous config.
        """        
        # Embeddable docker to hold GUI stuff
        self.docker = _egg.gui.Docker('Crimson RX Streamer', autosettings_path='data_streamer_window')
        self.docker.set_size([1000,700])
        
        # Top controls
        self.g_top = self.docker.place_object(_egg.gui.GridLayout(False))
        self.b_collect       = self.g_top.place_object(_egg.gui.Button("Go!", True, False))
        self.b_reset         = self.g_top.place_object(_egg.gui.Button("Reset"))
        
        # Info grid
        self.pb_buffer       = self.g_top.place_object(_egg.pyqtgraph.QtGui.QProgressBar())        
        self.b_overrun       = self.g_top.place_object(_egg.gui.Button('Buffer Overrun', True))
        self.b_timeout       = self.g_top.place_object(_egg.gui.Button('Timeout', True))        
        
        # button functions
        self.b_collect.signal_clicked.connect(self._b_collect_clicked)
        self.b_reset.signal_clicked.connect(self._b_reset_clicked)
        
        # Data collection timer
        self.timer_collect = _egg.gui.Timer(1)
        self.timer_collect.signal_tick.connect(self._timer_collect_tick)
        
        # Create Settings tree and parameters
        self.docker.new_autorow()
        self.g_stuff = self.docker.place_object(_egg.gui.GridLayout(False), alignment=0)
        self.settings = self.g_stuff.place_object(_egg.gui.TreeDictionary('data_streamer_settings')).set_width(300)       

        self.settings.add_parameter('Crimson/sample_rate',  1e5,  type='float', limits=(1e4,325e6), decimals=12, suffix='Hz', siPrefix=True, dec=True)
        self.settings.add_parameter('Crimson/post_command_settle',  1e-3, type='float', limits=(0,None), dec=True, suffix='s', siPrefix=True)      

        # Create the channel controls
        for n in range(4):

            # Assemble the channel name
            channel_name = 'Crimson/'+self._channel_names[n]
            
            # Create the channel and its parameters
            self.settings.add_parameter(channel_name, True, type='bool')
            self.settings.add_parameter(channel_name+'/center_frequency', 7e6, type='float', limits=(0,7e9),  decimals=12, suffix='Hz', siPrefix=True, step=0.1e6)
            self.settings.add_parameter(channel_name+'/gain',               0, type='float', limits=(0,31.5), step=0.5,    suffix=' dB')
        
            # Disable and hide the unused channels. Hopefully one day we can 
            # just have all 4 always visible and properly enable/disable them,
            # but currently it is not possible to stop and restart top_block
            # to update the connections without a freeze.
            if not n in self._rx_channels:
                self.settings.block_events()
                
                # Hide it.                
                self.settings.hide_parameter(channel_name)
                
                self.settings.unblock_events()

            # Otherwise, connect the signal changes
            else:
                self.settings.connect_signal_changed(channel_name+'/center_frequency', self._settings_center_frequency_changed)
                self.settings.connect_signal_changed(channel_name+'/gain',             self._settings_gain_changed)

        self.settings.add_parameter('Software/Acquisitions',                    0, type='int')
        self.settings.add_parameter('Software/Acquisitions/count',              0, type='int')
        self.settings.add_parameter('Software/Acquisitions/buffer_size',      1e3, type='int',   limits=(1,None), suffix=' packets', dec=True)
        self.settings.add_parameter('Software/Acquisitions/calibration',2.053e-05, type='float')
        self.settings.add_parameter('Software/Acquisitions/time',            0.01, type='float', limits=(1e-6,None), suffix='s', siPrefix=True, dec=True)
        self.settings.add_parameter('Software/Acquisitions/timeout',           10, type='float', limits=(0.1,None),  suffix='s', siPrefix=True, dec=True)        
        self.settings.add_parameter('Software/Acquisitions/keep_extra',     False, type='bool')        
        self.settings.add_parameter('Software/Acquisitions/trigger_channel', None, type='list', values=['None']+self._channel_names)        
        
        self.settings.add_parameter('Software/PSD',            False, type='bool')        
        self.settings.add_parameter('Software/PSD/window',     "None",type='str')
        self.settings.add_parameter('Software/PSD/pow2',       True,  type='bool')
        
        # Load and set previous settings (if any)        
        self.settings.load(ignore_errors=True, block_user_signals=True)
        
        # Add a PSD counter for each channel
        for n in range(4):
            self.settings.add_parameter('Software/PSD/averages_'+self._channel_names[n], 0, type='int', limits=(0,None))
        
        # Set it up so that any changes are sent to the Crimson
        self.settings.connect_signal_changed('Crimson/sample_rate',               self._settings_sample_rate_changed)
        self.settings.connect_signal_changed('Software/Acquisitions/buffer_size', self._settings_buffer_size_changed)
        self.settings.connect_signal_changed('Software/Acquisitions/time',        self._settings_acquisition_time_changed)
        
        # for saving
        self.settings.connect_any_signal_changed(self._any_setting_changed)

        # Plotter tabs
        self.tabs_plots = self.g_stuff.place_object(_egg.gui.TabArea(False, 'data_streamer_tabs_plots'), alignment=0)
        
        # Raw time trace    
        self.tab_raw = self.tabs_plots.add_tab('Raw')                
        self.p_raw = self.tab_raw.place_object(_egg.gui.DataboxPlot("*.raw", "data_streamer_p_raw"), alignment=0)
        self.p_raw_trigger_level = _egg.pyqtgraph.InfiniteLine(0, movable=True, angle=0)
        self.p_raw.ROIs.append([self.p_raw_trigger_level])
        self.p_raw.load_gui_settings()
        
        # "Calibrated" time trace
        self.tab_volts = self.tabs_plots.add_tab('Volts')
        self.p_volts   = self.tab_volts.place_object(_egg.gui.DataboxPlot("*.volts", "data_streamer_p_volts"), alignment=0)
        
        # FFT
        self.tab_fft = self.tabs_plots.add_tab('FFT')
        self.p_fft   = self.tab_fft.place_object(_egg.gui.DataboxPlot("*.fft", "data_streamer_p_fft", autoscript=3), alignment=0)
        
        # PSD
        self.tab_psd = self.tabs_plots.add_tab('PSD')        
        self.p_psd   = self.tab_psd.place_object(_egg.gui.DataboxPlot("*.psd", "data_streamer_p_psd", autoscript=2), alignment=0)
        
        # load last selected tab
        self.tabs_plots.load_gui_settings()
        
        
            

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

        # Sample rate is global
        self.settings.set_value('Crimson/sample_rate', self.crimson.get_sample_rate(), True)
        for n in range(len(self._rx_channels)):
            channel_name = self.crimson._channel_names[self._rx_channels[n]]
            self.settings.set_value('Crimson/'+channel_name+'/center_frequency', self.crimson.get_center_frequency(n), True)
            self.settings.set_value('Crimson/'+channel_name+'/gain',             self.crimson.get_gain(n),             True)
        

    def _settings_sample_rate_changed(self, *a): 
        self.timer_collect.stop()
        
        # Tell the Crimson then get the actual value
        self.crimson.set_sample_rate(self.settings['Crimson/sample_rate'])
        self.settings.set_value('Crimson/sample_rate', self.crimson.get_sample_rate(), True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/post_command_settle'])
        self.timer_collect.start()
        
    def _settings_center_frequency_changed(self, *a):
        self.timer_collect.stop()

        # get the channel info
        name = a[0].parent().name()
        key  = 'Crimson/'+name+'/center_frequency' 
        n    = self._rx_channels.index(self._channel_names.index(name))
        
        # Tell the Crimson then get the actual value
        self.crimson.set_center_frequency(self.settings[key], n)
        self.settings.set_value(name, self.crimson.get_center_frequency(n), True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/post_command_settle'])
        self.timer_collect.start()
            
    def _settings_gain_changed(self, *a):
        self.timer_collect.stop()
        
        # get the channel info
        name = a[0].parent().name()
        key  = 'Crimson/'+name+'/gain' 
        n    = self._rx_channels.index(self._channel_names.index(name))
    
        # Tell the Crimson then get the actual value
        self.crimson.set_gain(self.settings[key], n)
        self.settings.set_value(name, self.crimson.get_gain(n), True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/post_command_settle'])        
        self.timer_collect.start()

    def _settings_acquisition_time_changed(self, *a):
        self.reset_acquisition()

    def _settings_buffer_size_changed(self, *a):
        self.crimson.buffer.set_size(self.settings['Software/Acquisitions/buffer_size'])
    
    def start(self): return self._top_block.start()    

    def _get_gui_enabled_channels(self):
        """
        Returns a list of channel indices for each one that is both
        present and checked.
        """
        x = []
        for n in self._rx_channels: # [0,3], e.g.
            if self.settings['Crimson/'+self._channel_names[n]]: x.append(n)
        return x

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
        
        # Get the data        
        data, overruns, timeout_reached = self.crimson.get_data(N, 
                    keep_all = self.settings['Software/Acquisitions/keep_extra'],
                    timeout  = self.settings['Software/Acquisitions/timeout'])
        
        # If we timed out, quit
        if timeout_reached: 
            self.b_timeout.set_checked(True)
            return

        # Remove old data
        self.p_raw.clear()
        self.p_volts.clear()
        self.p_fft.clear()
        
        # Default for the time column
        t = _n.arange(0, (len(data[0][0])-0.5)*dt, dt)
        self.p_raw['t']   = t
        self.p_volts['t'] = t
        
        # Store the data for each enabled channel
        enabled = self._get_gui_enabled_channels()         # e.g. [0,2]
        for n in enabled: 
            
            # Get the column labels
            nameX = self._channel_names[n]+"/X"            
            nameY = self._channel_names[n]+"/Y"            
            
            # Store the raw data by name             
            self.p_raw[nameX] = data[self._rx_channels.index(n)][0]
            self.p_raw[nameY] = data[self._rx_channels.index(n)][1]
            
            # Calibrate to volts
            Vx = self.p_volts[nameX] = self.p_raw[nameX]*self.settings['Software/Acquisitions/calibration']*10**(-0.05*self.settings['Crimson/RXA/gain']) 
            Vy = self.p_volts[nameY] = self.p_raw[nameY]*self.settings['Software/Acquisitions/calibration']*10**(-0.05*self.settings['Crimson/RXA/gain'])  
            
            # Check for overruns
            if (_n.array(overruns) > 0).any(): self.b_overrun.set_checked()   
        
            # PSD
            if self.settings['Software/PSD']:

                psd_counter_name = 'Software/PSD/averages_'+self._channel_names[n]                
                
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
                        f = f+self.settings['Crimson/'+self._channel_names[n]+'/center_frequency']                
                        P = 0.5*abs(fft)**2 / df
                                                    
                        # Get the name of the column
                        namef = self._channel_names[n]+"_f(Hz)"
                        nameF = self._channel_names[n]+"_FFT(V)"
                        nameP = self._channel_names[n]+"_PSD(V^2/Hz)"
                        
                        # Store the single-shot fft
                        self.p_fft[namef] = f
                        self.p_fft["Re["+nameF+"]"] = _n.real(fft)
                        self.p_fft["Im["+nameF+"]"] = _n.imag(fft)
                        
                        # Make sure we have the data arrays (e.g. after a reset)
                        if not nameP in self.p_psd.ckeys: 
                            self.p_psd[namef] = f
                            self.p_psd[nameP] = P
        
                        # Do the average
                        n_psd = self.settings[psd_counter_name]
                        self.p_psd[nameP] = (self.p_psd[nameP]*n_psd + P)/(n_psd+1)

                        # increment the psd counter
                        self.settings[psd_counter_name] += 1
        
        # If we're not supposed to trigger, use time from 0 to max            
        trigger_channel = self.settings['Software/Acquisitions/trigger_channel']
        if not trigger_channel == 'None':
            
            # Choose the raw channel, and use the x quadrature for now
            x = self.p_raw[trigger_channel+"/X"]                
            
            # Trigger level
            l = self.p_raw_trigger_level.value()

            # All values where we're below the trigger level
            ns0 = _n.where(x<l)[0]            
            
            # If we have any values below the trigger and we're supposed to 
            # trigger, find the crossing
            if len(ns0):
                
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
        self.settings.send_to_databox_header(self.p_fft)
        self.p_fft.plot()
        
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
                
        # If the button's enabled and we've enabled at least one channel
        if self.b_collect.is_checked() and len(self._get_gui_enabled_channels()): 
            self.get_and_plot_data()


    def hide(self):                            return self.docker.hide()
    def show(self, block_command_line=False):  return self.docker.show(block_command_line)
    
    def reset_acquisition(self, reset_psd=True):
        
        # Resets counters, clears plots, clears buffers     
        self.settings['Software/Acquisitions/count'] = 0
        self.p_raw.clear()
        self.p_volts.clear()
        self.p_fft.clear()

        self.b_overrun.set_checked(False)
        if not self._simulation_mode: self.crimson.buffer.flush()
                
        if reset_psd:
            for n in range(4): self.settings['Software/PSD/averages_'+self._channel_names[n]] = 0
            self.p_psd.clear()
        




# Code to run if imported / executed as a standalone script
if __name__ == '__main__': self = rx_streamer([0,1,2,3], block_command_line=True)
    
