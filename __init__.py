from gnuradio import gr  as _gr
from gnuradio import uhd as _uhd

import thread      as _thread
import time        as _time
import numpy       as _n



class _data_buffer(_gr.sync_block):
    """
    A simple GNU Radio data buffer.
    """
    
    def __init__(self, size=1000, channels=2):
        """
        Thread-safe data buffer that will hold the number of data packets
        specified by "size" in the specified number of channels (size can be
        an integer or list of sizes, one for each channel). "payload" is the
        number of bytes per packet, which for the crimson means
        
           (payload-16)/4/channels samples per packet
        """
        
        # Initialize the base object
        _gr.sync_block.__init__(self, name="_data_buffer", 
                in_sig=[(_n.int16,2)]*channels, out_sig=None)
            
        # Some internal variables to keep track of
        self._channel_count = channels
         
        # Thread locks for the elements read / written on the fly
        self._buffer_lock          = _thread.allocate_lock()
        self._buffer_overruns_lock = _thread.allocate_lock()
        self._size_lock            = _thread.allocate_lock()                
        
        # Creates self._buffer and self._buffer_overruns to track
        # the actual data and count overruns, respectively
        self._buffer          = []
        self._buffer_overruns = []
        self.flush()
        
        # Save this for size checking
        self.set_size(size)        
        
        # Counter for diagnostics
        self._n = 0

    def __len__(self): return max(self.get_packet_counts())         
        
    def flush(self, reset_overruns=True):
        """
        Clears the buffer, returning it and the overruns count.
        """
        # We can't do [[]]*N because it doesn't create separate instances!
        self._buffer_lock.acquire()
        self._buffer_overruns_lock.acquire()

        # Keep these safe, making a copy of overruns because we don't
        # necessarily unlink it from the internal object
        b = self._buffer
        o = list(self._buffer_overruns)

        # Reset the buffer
        self._buffer = list()
        for n in range(self.get_channel_count()): self._buffer.append(list())
        
        # Reset the overruns if we're supposed to
        if reset_overruns: self._buffer_overruns = [0]*self.get_channel_count()

        self._buffer_lock.release()
        self._buffer_overruns_lock.release()        
            
        # Send the data back
        return b, o

    def get_channel_count(self): 
        """
        Returns the number of channels.
        """
        return self._channel_count

    def get_overruns(self, reset=False):
        """
        Returns the number of overruns. If specified, resets the overrun
        counts to zero.
        """
        self._buffer_overruns_lock.acquire()
        
        # Make a copy to avoid referencing issues
        x = self._buffer_overruns
        if reset: self._buffer_overruns = [0]*self.get_channel_count()
        
        self._buffer_overruns_lock.release()     
        
        return x
    
    def get_data(self, samples=1024, keep_all=False, timeout=1.0):
        """
        Waits for enough packets to have at least the specified number of 
        samples (or the timeout), then returns all packets, the number of 
        overruns for each channel, and the timeout status (True for timeout), 
        clearing the buffer and resetting the overrun count.
        
        samples       Specifies the minimum number of samples to wait for.
                      Can also be a list of numbers matching the number of
                      channels.
        
        keep_all      If False, this will keep only the newest packets
                      such that there are at least the specified number of 
                      samples. If True, returns whatever was in the buffer
                      (this is ideal for taking long continuous data).
        
        timeout       Specifies the "give up" time
        
        Returns [[x1,y1],[x2,y2],...], overrun_counts, timeout_reached
        
        where x1, y1, x2, etc are numpy arrays for the two quadratures of each
        channel.
        """
        
        # Make sure samples is a numpy array (for easy comparison)
        if not type(samples) == list: samples = [samples]*self.get_channel_count()
        samples = _n.array(samples)

        # create the thing to hold the data
        d = []
        for n in range(self.get_channel_count()): d.append([[],[]])
        sample_counts = _n.zeros(self.get_channel_count())
        overruns      = _n.zeros(self.get_channel_count())        

        # Get the starting time for measuring timeout
        t0 = _time.time()
                
        # Wait until we get (at least) the right number of samples in each channel
        while _time.time()-t0 < timeout and (sample_counts < samples).any(): 
            
            # Get and clear the buffer
            ps, o = self.flush()
            
            # Convert packets (p) into data [[x1,y1],[x2,y2],...]        
            for n in range(self.get_channel_count()):
                
                # If we have packets on this channel
                if len(ps[n]):
                    
                    # Loop over the packets for this channel, store the data
                    # and increment the sample count
                    for p in ps[n]:
                        d[n][0].append(p[:,0])
                        d[n][1].append(p[:,1])
                        
                        # Update the sample counts
                        sample_counts[n] += len(d[n][0][-1])
        
                # Update the overruns count
                overruns[n] += o[n]
        
            # Save some cpu cycles
            _time.sleep(0.001)
            
            # Update the GUI (if the user has supplied a function for this)
            self.process_events()
        
        # At this point we've got at least enough or have timed out.

        # If we reached the timeout
        timeout_reached = False
        if _time.time()-t0 >= timeout: 
            print("WARNING: get_data() timeout")
            timeout_reached = True

        # Tidy up the data
        for n in range(self.get_channel_count()):
            
            # Concatenate the data lists (we only want to do this once)
            if len(d[n][0]):           
                d[n][0] = _n.concatenate(d[n][0])
                d[n][1] = _n.concatenate(d[n][1])
                
            # Throw away extras if we're supposed to
            if not keep_all: 
                # Keep the most recent                
                d[n][0] = d[n][0][len(d[n][0])-samples[n]:]
                d[n][1] = d[n][1][len(d[n][1])-samples[n]:]

        return d, overruns, timeout_reached

    def process_events(self):
        """
        Overwrite this function with a valid process-events function in order
        to make your gui responive while running a long get_data().
        """
        
    def reset_overruns(self):
        """
        Sets the number of overruns to zero.
        """
        self.get_overruns(True)

    def set_size(self, size):
        """
        Changes the size of the buffer; size can be an integer or list of sizes,
        one for each channel.
        """
        # Make sure we have a size for each channel
        if not type(size)==list: size = [size]*self.get_channel_count()

        # Idiot proofing
        if not len(size) == self.get_channel_count():
            print "ERROR set_size(): size list length doesn't match number of channels."
            return
        
        self._size_lock.acquire()
        self._size = size
        self._size_lock.release()
    
    def get_size(self):
        """
        Returns the maximum size of the buffer.
        """
        self._size_lock.acquire()
        x = list(self._size)
        self._size_lock.release()
        return x

    def get_packet_counts(self):
        """
        Returns the current number of packets in each channel of the buffer. 
        len(self) will return the maximum of this list.
        """
        Ns = []
        self._buffer_lock.acquire()
        for ps in self._buffer: Ns.append(len(ps))
        self._buffer_lock.release()
        
        return Ns

    def work(self, input_items, output_items):
        
        # Acquire the buffer lock for the whole function.
        self._buffer_lock.acquire()

        # Loop over channels        
        for n in range(len(input_items)):
            
#            if self._n%500==1: 
#                print "work", self._n, len(input_items), len(input_items[0]), len(self._buffer[0])
#            self._n +=1            
            
            # Append the new data
            self._buffer[n].append(_n.array(input_items[n]))

            self._size_lock.acquire()
        
            # If we've overrun
            while len(self._buffer[n]) > self._size[n]:
    
                # Remove the oldest data point            
                self._buffer[n].pop()
                
                # Increment the overrun
                self._buffer_overruns_lock.acquire()
                self._buffer_overruns[n] += 1
                self._buffer_overruns_lock.release()
            
            self._size_lock.release()
        
        # release the buffer lock
        self._buffer_lock.release()        
        
        return len(input_items[0])

class crimson():
    """
    Scripted interface to the Crimson (no GUI). Typical workflow:
    
      c = crimson()

      c.start()
       
      c.get_data()
      
    Parameters
    ----------    
      rx_channels : list    
                      list of enabled RX channels      
      buffer_size : int    
                      how many packets (default 346 samples) for the buffer to keep
      simulation_mode : bool   
                      don't connect to the crimson, but send fake data (not well-developed yet)
      process_events : function
                      optional process_events function that will be run during
                      long "get_data()" operations to allow, e.g., GUI updates.

      
    """
    
    def __init__(self, rx_channels=[0,1], buffer_size=500, simulation_mode=False, process_events=None):
        if simulation_mode: print("WARNING: Simulation mode!")

        self._enabled_rx_channels = rx_channels
        self._top_block           = None
        self._crimson             = None
        self._simulation_mode     = simulation_mode
        self._connected_points    = []
        self._channel_names       = ["RXA","RXB", "RXC", "RXD"]
        self.buffer               = _data_buffer(buffer_size)
        
        # Enable the specified channels and buffer
        if not self._simulation_mode:
            self._initialize_crimson(rx_channels, buffer_size)        
        
        # Provide the function to process events during get_data()
        if not process_events == None: self.buffer.process_events = process_events
        
        return        
        
    def _initialize_crimson(self, channels=[0,2], buffer_size=500):
        """
        Enables the specified channels.
        
        Behind the scenes, this creates a top block and a GNU radio UHD 
        usrp_source object, then connects the usrp_source to a data buffer 
        (self.buffer).
        """
        
        # Create the buffer
        self.buffer = _data_buffer(buffer_size, len(channels))
        
        print("pycrimson - Initializing GNU Radio with these channels enabled:")
        for n in channels: print("  "+self._channel_names[n])
        
        # Create the gnuradio top block
        self._top_block = _gr.top_block("Crimson")

        # Create the crimson data faucet  
        self._crimson = _uhd.usrp_source("crimson",
        	_uhd.stream_args(cpu_format="sc16", args='sc16', channels=(channels)))

        # Connect the faucet to the buffer
        for n in range(len(channels)): self._connect(((self._crimson, n), (self.buffer, n)))

    
    def _lock(self):   return self._top_block.lock()
    def _unlock(self): return self._top_block.unlock()

    def _connect(self, *point_pairs):
        """
        Takes any number of point pairs (a,b) and connects them in sequence.
        """
        for p in point_pairs:

            # Store the connections for safe keeping
            self._connected_points.append((p[0],p[1]))

            # Connect them
            self._top_block.connect(p[0],p[1])

    def _disconnect(self, *point_pairs):
        """
        Takes any number of point pairs (a,b) and disconnects them. 
        """
        for p in point_pairs:
            if not p in self._connected_points: 
                print "ERROR _disconnect(): not in list of connected points."

            # Find the index of this connection
            i = self._connected_points.index(p)
            
            # Disconnect and forget
            self._top_block.disconnect(*p)
            self._connected_points.pop(i)
    
    def _disconnect_all(self):
        """
        Disconnects all crimson points from the buffer.
        """
        self._disconnect(*self._connected_points)
       

    def get_enabled_rx_channels(self):
        """
        Returns the list of enabled rx channels.
        """
        return self._enabled_rx_channels
    
    def get_sample_rate(self):
        """
        Gets the sampling rate (global).
        """
        if self._simulation_mode: return 17.777e6
        return self._crimson.get_samp_rate()
    
    def get_center_frequency(self, channel=0):
        """
        Returns the current center frequency of the specified channel.
        """
        if self._simulation_mode: return 17.777e6
        return self._crimson.get_center_freq(channel)
    
    def get_data(self, samples=1024, keep_all=False, timeout=1.0):
        """
        Shortcut to self.buffer.get_data()
        """
        return self.buffer.get_data(samples, keep_all, timeout)
    
    def get_gain(self, channel=0):
        """
        Returns the current gain of the specified channel.
        
        NOTE: IS CURRENTLY INCORRECT FOR HIGH-BAND OPERATION.
        """
        if self._simulation_mode: return 7
        return (126-self._crimson.get_gain(channel))/4

    def set_sample_rate(self, sample_rate=1e5):
        """
        Gets the sampling rate (global).
        """
        if self._simulation_mode: return 17.777e6
        return self._crimson.set_samp_rate(sample_rate)
    
    def set_center_frequency(self, center_frequency=7e7, channel=0):
        """
        Returns the current center frequency of the specified channel.
        """
        if self._simulation_mode: return 17.777e6
        return self._crimson.set_center_freq(center_frequency, channel)
    
    def set_gain(self, gain=0, channel=0):
        """
        Returns the current gain of the specified channel.
        
        NOTE: IS CURRENTLY INCORRECT FOR HIGH-BAND OPERATION.
        """
        if self._simulation_mode: return 7
        return self._crimson.set_gain(gain,channel)


    def start(self):
        """
        Start the data flow.
        """
        if self._simulation_mode: return self
        self._top_block.start()
        return self
    
    def stop(self):
        """
        Stop the data flow.
        """
        if self._simulation_mode: return self
        self._top_block.stop()
        return self

if __name__ == '__main__':   

# Super simple example of getting data
#    self = crimson([0,2])
#    self.start()
#    d, o, z = self.get_data()


# Basic egg GUI to stream data.
   
    import spinmob.egg as _egg
    import pycrimson   as _pc
    
    channels = [0,2]
    single_shot = False

    # Basic GUI
    w = _egg.gui.Window(autosettings_path="w")
    t = _egg.gui.Timer(1, single_shot=single_shot)
    p = w.place_object(_egg.gui.DataboxPlot())  
    w.load_gui_settings()
    
    # Crimson
    c = _pc.crimson(channels,1000)
    c.set_sample_rate(1e4)
    for n in range(len(channels)): c.set_center_frequency(7e6, n)
    c.start()
        
    # Timer tick function
    def tick(*a):
        
        # Get the data
        ds,o,z = c.get_data(1000, False)
        
        # Write the plotter columns of data
        p['t'] = range(len(ds[0][0]))
        for n in range(len(ds)):
            p['x'+str(n)] = ds[n][0]
            p['y'+str(n)] = ds[n][1]
            
        # Plot and update the GUI
        p.plot()
        w.process_events()
    
    t.signal_tick.connect(tick)

    ## Code that cleans up the front-panel RX LEDs
    def reset_rx(*a):
        t.stop()
        c.stop()
    w.event_close = reset_rx

    # Show the window and start the timer.
    w.show(); t.start()
    
