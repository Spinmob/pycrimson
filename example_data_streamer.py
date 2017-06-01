from gnuradio import gr  as _gr
from gnuradio import uhd as _uhd

import thread      as _thread
import time        as _time
import numpy       as _n

import spinmob     as _s
import spinmob.egg as _egg


class jax_data_streamer():

    def __init__(self):
        
        # Create the gnuradio top block
        self._top_block = _gr.top_block("Data Streamer")        
        
        # Create the faucet from the Crimson
        self._crimson = _uhd.usrp_source("crimson",
        	_uhd.stream_args(cpu_format="sc16", args='sc16', channels=([0])))
        
        # Buffer defined above the __init__        
        # Connect the faucet to the buffer
        self._top_block.connect((self._crimson, 0), (self.buffer, 0))

        # Build the gui, load previous settings, connect signals
        self._build_gui()
        
        # Update settings from the Crimson itself
        self._get_settings_from_crimson()      
    
        # Set the buffer size
        self.buffer.set_size(self.settings['Software/buffer_size'])        
        
        # Show it!
        self.show()
        
        # Start the crimson.
        self._top_block.start()
        
        # Start the collection timer
        self.timer_collect.start()
    
    def _build_gui(self):
        """
        Places all the controls and plots etc, loads previous config.
        """        
        # Embeddable docker to hold GUI stuff
        self.docker = _egg.gui.Docker('RXA via Buffer', autosettings_path='rxa_via_buffer_window')
        
        # Top controls
        self.g_top = self.docker.place_object(_egg.gui.GridLayout(False))
        self.b_collect       = self.g_top.place_object(_egg.gui.Button('Collect Data!' , True, True))
        self.b_trigger       = self.g_top.place_object(_egg.gui.Button('Trigger',        True, False))        
        self.b_reset         = self.g_top.place_object(_egg.gui.Button('Reset'))
        self.n_trace_counter = self.g_top.place_object(_egg.gui.NumberBox(0))
        self.n_psd_counter   = self.g_top.place_object(_egg.gui.NumberBox(0))
        
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
        self.settings = self.g_stuff.place_object(_egg.gui.TreeDictionary('rxa_via_buffer_settings')).set_width(300)       
        self.settings.add_parameter('Crimson/sample_rate',      1e6,  type='float', limits=(1e4,325e6), decimals=12, suffix='Hz', siPrefix=True, dec=True)
        self.settings.add_parameter('Crimson/center_frequency', 7e6,  type='float', limits=(0,6e9),     decimals=12, suffix='Hz', siPrefix=True, step=0.1e6)
        #self.settings.add_parameter('bandwidth',        1e5, type='float', limits=(0,325e6/2), decimals=12, suffix='Hz', siPrefix=True)        
        self.settings.add_parameter('Crimson/gain',             0,    type='float', limits=(0,31.5), step=0.5,    suffix=' dB')
        self.settings.add_parameter('Crimson/settle_time',      1e-3, type='float', limits=(0,None), dec=True, suffix='s', siPrefix=True)      
        self.settings.add_parameter('Software/buffer_size',      1e3, type='int',   limits=(1,None), suffix=' packets', dec=True)
        self.settings.add_parameter('Software/voltage_cal',2.053e-05, type='float')
        self.settings.add_parameter('Software/sample_time',     0.01, type='float', limits=(1e-6,None), suffix='s', siPrefix=True, dec=True)
        self.settings.add_parameter('Software/PSD_enabled',    False, type='bool')        
        self.settings.add_parameter('Software/PSD_averages',       0, type='float')
        self.settings.add_parameter('Software/PSD_window',     "none",type='str')
        self.settings.add_parameter('Software/PSD_pow2',       True,  type='bool')

        # Load and set previous settings (if any)        
        self.settings.load()
        
        # Set it up so that any changes are sent to the Crimson
        self.settings.connect_signal_changed('Crimson/sample_rate',      self._settings_sample_rate_changed)
        self.settings.connect_signal_changed('Crimson/center_frequency', self._settings_center_frequency_changed)
        #self.settings.connect_signal_changed('bandwidth',        self._settings_bandwidth_changed)
        self.settings.connect_signal_changed('Crimson/gain',             self._settings_gain_changed)
        self.settings.connect_signal_changed('Software/buffer_size',     self._settings_buffer_size_changed)
        
        # for saving
        self.settings.connect_any_signal_changed(self._any_setting_changed)

        # Plotter tabs
        self.tabs_plots = self.g_stuff.place_object(_egg.gui.TabArea(False, 'rxa_via_buffer_tabs_plots'), alignment=0)
        
        # Raw time trace    
        self.tab_raw = self.tabs_plots.add_tab('Raw')                
        self.p_raw = self.tab_raw.place_object(_egg.gui.DataboxPlot("*.txt", "rxa_via_buffer_p_raw"), alignment=0)
        self.p_raw_trigger_level = _egg.pyqtgraph.InfiniteLine(0, movable=True, angle=0)
        self.p_raw.ROIs.append([self.p_raw_trigger_level])
        self.p_raw.load_gui_settings()
        
        # "Calibrated" time trace
        self.tab_volts = self.tabs_plots.add_tab('Volts')
        self.p_volts   = self.tab_volts.place_object(_egg.gui.DataboxPlot("*.txt", "rxa_via_buffer_p_volts"), alignment=0)
        
        # PSD
        self.tab_psd = self.tabs_plots.add_tab('PSD')        
        self.p_psd   = self.tab_psd.place_object(_egg.gui.DataboxPlot("*.txt", "rxa_via_buffer_p_psd"), alignment=0)
        
        # Info grid
        self.docker.new_autorow()
        self.g_bottom = self.docker.place_object(_egg.gui.GridLayout(False))        
        self.pb_buffer       = self.g_bottom.place_object(_egg.pyqtgraph.QtGui.QProgressBar())        
        self.b_overrun       = self.g_bottom.place_object(_egg.gui.Button('Buffer Overrun', True))
        
        # load last selected tab
        self.tabs_plots.load_gui_settings()

    def _any_setting_changed(self, *a): 
        
        # Reset the PSD averages and save the settings
        self.reset_acquisition()        
        self.settings.save()

    def _b_collect_clicked(self, *a):
        
        # If we're taking new data, reset the PSD, flush the buffer, reset the overrun        
        if self.b_collect.is_checked(): self.reset_acquisition()
            
    def _b_reset_clicked(self, *a): 
        self.reset_acquisition()

    def _get_settings_from_crimson(self):
        """
        Asks the crimson what all the settings currently are, and updates
        the GUI. This also blocks all user-defined signals for changing the GUI.
        """
        self.settings.set_value('Crimson/sample_rate',      self._crimson.get_samp_rate(),   True)
        self.settings.set_value('Crimson/center_frequency', self._crimson.get_center_freq(), True)
        self.settings.set_value('Crimson/gain',   (126-self._crimson.get_gain())/4,          True)
        

    def _settings_sample_rate_changed(self, *a): 
        self.timer_collect.stop()
        
        # Tell the Crimson then get the actual value
        self._crimson.set_samp_rate(self.settings['Crimson/sample_rate'])
        self.settings.set_value('Crimson/sample_rate', self._crimson.get_samp_rate(), True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/settle_time'])
        self.timer_collect.start()
        
    def _settings_center_frequency_changed(self, *a):
        self.timer_collect.stop()
        
        # Tell the Crimson then get the actual value
        self._crimson.set_center_freq(self.settings['Crimson/center_frequency'])
        self.settings.set_value('Crimson/center_frequency', self._crimson.get_center_freq(), True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/settle_time'])
        self.timer_collect.start()
            
    def _settings_gain_changed(self, *a):
        self.timer_collect.stop()
        
        # Tell the Crimson then get the actual value
        self._crimson.set_gain(self.settings['Crimson/gain'])
        self.settings.set_value('Crimson/gain', (126-self._crimson.get_gain())/4, True)
        self.reset_acquisition()        
        
        self.docker.sleep(self.settings['Crimson/settle_time'])        
        self.timer_collect.start()

    def _settings_buffer_size_changed(self, *a):
        self.buffer.set_size(self.settings['Software/buffer_size'])
    
    def start(self): return self._top_block.start()
        
    
    def _timer_collect_tick(self, *a):
        """
        Called every time the data collection timer ticks.
        """
        # Update the number of packets
        self.pb_buffer.setValue(_n.round(100.0*len(self.buffer)/self.buffer.get_size()))

        # Update the user if the buffer overran
        if self.buffer.get_overruns(): self.b_overrun.set_checked(True)
                            
        if self.b_collect.is_checked(): 

            # Get the time spacing
            dt = 1.0/self.settings['Crimson/sample_rate']

            # Total time
            T = self.settings['Software/sample_time']                        
            N = int(_n.round(T/dt))
            
            # Get all the packets we need            
            packets, overruns = self.get_packets(N, 1)
            
            # If we timed out
            if len(packets) == 0: 
                return

            # Success! Increment the counter
            self.n_trace_counter.increment()
            
            # Assemble the packets            
            [x,y] = _n.concatenate(packets).transpose()

            # If we're supposed to trigger, find the first index at which
            # it crosses the threshold
            l   = self.p_raw_trigger_level.value()
            
            ns0 = _n.where(x<l)[0]
            if len(ns0) and self.b_trigger.is_checked():
                
                # first value below trigger
                n0  = ns0[0]

                # indices after n0 where we're above the trigger
                ns1 = _n.where(x[n0:]>l)[0]+n0

                # If we have a crossing                
                if(len(ns1)):  t  = _n.arange(-ns1[0]*dt,(-ns1[0]+len(x)-0.5)*dt, dt) 
                else:          t = _n.arange(0, (len(x)-0.5)*dt, dt)    
            else:
                t = _n.arange(0, (len(x)-0.5)*dt, dt)

            # Plot Raw
            self.p_raw['t'] = t
            self.p_raw['Nx'] = x
            self.p_raw['Ny'] = y
            self.settings.send_to_databox_header(self.p_raw)            
            self.p_raw.plot()

            # Plot Volts
            self.p_volts['t'] = t
            Vx = x*self.settings['Software/voltage_cal']*10**(-0.05*self.settings['Crimson/gain']) 
            Vy = y*self.settings['Software/voltage_cal']*10**(-0.05*self.settings['Crimson/gain'])
            self.p_volts['Vx'] = Vx
            self.p_volts['Vy'] = Vy
            self.settings.send_to_databox_header(self.p_volts)            
            self.p_volts.plot()
            
            # PSD
            if self.settings['Software/PSD_enabled']:
                
                # If we want the FFT to be extra efficient
                if self.settings['Software/PSD_pow2']: N = 2**int(_n.log2(N))              
                
                # Loop over the data to get as many PSD's as possible from it                
                for k in range(int(len(t)/int(N))):
                    
                    # Do the PSD
                    result = _s.fun.fft(t[k*N:(k+1)*N], Vx[k*N:(k+1)*N]+1j*Vy[k*N:(k+1)*N], 
                                          pow2=False, window=self.settings['Software/PSD_window'])
                    
                    # Happens if there is an invalid window
                    if result == None: return                

                    # Convert to PSD
                    f, fft = result
                    df = f[1]-f[0]
                    f = f+self.settings['Crimson/center_frequency']                
                    P = 0.5*abs(fft)**2 / df
                    
                    # Reset the counter if the length or the length has changed
                    if len(self.p_psd) and not len(self.p_psd[0]) == len(f): self.reset_acquisition()
    
                    # Make sure we have the data arrays
                    if len(self.p_psd) == 0: 
                        self.p_psd['f'] = f
                        self.p_psd['P'] = P
    
                    # Do the average
                    n = self.n_psd_counter.get_value()
                    self.p_psd['f'] = f
                    self.p_psd['P'] = (self.p_psd['P']*n + P)/(n+1)
                    self.settings.send_to_databox_header(self.p_psd)
                    self.p_psd.plot()
                
                    # increment the psd counter
                    if self.settings['Software/PSD_averages'] == 0 or n+1 <= self.settings['Software/PSD_averages']:
                        self.n_psd_counter.increment()
                    
                    # Otherwise we're done!                    
                    else:
                        self.b_collect.set_checked(False)
                    
            self.docker.process_events()
            
            
    def hide(self):  return self.docker.hide()
    def show(self):  return self.docker.show()
    
    def reset_acquisition(self):
        
        # Resets counters, clears plots, clears buffers        
        self.n_psd_counter.set_value(0)
        self.n_trace_counter.set_value(0)
        self.b_overrun.set_checked(False)
        self.buffer.flush_buffer()
        self.p_psd.clear()
        self.p_raw.clear()
        self.p_volts.clear()



if __name__ == '__main__':   
    self = crimson()
    self.enable_rx_channels([0,1])
    self.start()
    self.get_packets()
