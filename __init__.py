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
        an integer or list of sizes, one for each channel).
        """
        
        # Initialize the base object
        _gr.sync_block.__init__(self, name="_data_buffer", 
                in_sig=[(_n.int16,2)]*channels, out_sig=None)
        
        # Create the buffer
        self._buffer  = [[]]*channels
        self._buffer_lock = _thread.allocate_lock()

        # Save this for size checking
        self._size_lock = _thread.allocate_lock()                
        self.set_size(size)        
        
        # Flag for whether the buffer has overrun since the last reset 
        self._buffer_overruns = [0]*self.get_channels()
        self._buffer_overruns_lock = _thread.allocate_lock()
        
    
    def __len__(self):
        """
        Returns the number of packets currently in the buffer.
        """
        Ns = []
        self._buffer_lock.acquire()
        for n in range(self.get_channels()): Ns.append(len(self._buffer[n]))
        self._buffer_lock.release()
        return Ns
    
    def _at_least_n_packets(self, packets, n=1):
        """
        Checks the supplied packet list to make sure there is at least n
        packets in each channel. n can be an integer or list of integers
        (one for each channel).
        """
        # Make sure we have a number for each channel
        if not type(n)==list: n = [n]*self.get_channels()        

        # Check each channel        
        for i in range(len(packets)): 
            if len(packets[i]) < n[i]: return False
        return True
    
    def flush_buffer(self):
        """
        Clears the buffer.
        """
        self._buffer_lock.acquire()
        self._buffer = [[]]*self.get_channels()
        self._buffer_lock.release()
        
        self._buffer_overruns_lock.acquire()
        self._buffer_overruns = [0]*self.get_channels()
        self._buffer_overruns_lock.release()        

    def get_channels(self): 
        """
        Returns the number of channels.
        """
        return len(self._buffer)

    def get_overruns(self, reset=False):
        """
        Returns the number of overruns. If specified, resets the overrun
        counts to zero.
        """
        self._buffer_overruns_lock.acquire()
        
        # Make a copy to avoid referencing issues
        x = list(self._buffer_overruns)
        if reset: self._buffer_overruns = [0]*self.get_channels()
        
        self._buffer_overruns_lock.release()     
        
        return x
    
    def get_packets(self, samples=1, keep_all=False, timeout=1.0):
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
        """

        # Make sure there is a number for each channel.
        if not type(samples) == list: samples = [samples]*self.get_channels()

        # Idiot proofing
        if not len(samples) == self.get_channels():
            print "ERROR get_packets(): length of samples is not equal to the number of channels."
        
        t0 = _time.time()
        packets = [[]]*self.get_channels()

        # Timeout flag
        timeout_reached = False        
        
        # Wait until we get (at least) our first packet
        while _time.time()-t0 < timeout and not self._at_least_n_packets(packets,1): 
            
            self._buffer_lock.acquire()
            
            # Append all data from each channel and clear the buffer
            for n in range(self.get_channels()): 
                packets[n] = packets[n] + self._buffer[n]
                self._buffer[n] = []

            self._buffer_lock.release()

            # Save some cpu cycles
            _time.sleep(0.001)
        
        # Find required number of packets
        packet_counts = []
        for n in range(len(packets)):
            packet_counts.append(int(_n.ceil(1.0*samples[n]/len(packets[n][0]))))        
        
        # Now we can find the size of the packets, so wait for the specified 
        # number of samples.
        while _time.time()-t0 < timeout and not self._at_least_n_packets(packets, packet_counts): 
            
            self._buffer_lock.acquire()
            
            # Return all data and clear the buffer
            packets[n] = packets[n] + self._buffer[n]
            self._buffer[n] = []

            self._buffer_lock.release()

            # Save some cpu cycles
            _time.sleep(0.001)

        # Get the number of overruns and reset
        overruns = self.get_overruns(True)

        # Throw away extras if we're supposed to
        if not keep_all: 
            for n in range(len(packets)):            
                
                # Keep the most recent                
                packets[n] = packets[n][len(packets[n])-packet_counts[n]:]

        # If we reached the timeout
        if _time.time()-t0 >= timeout: 
            print "WARNING: get_packets() timeout"
            timeout_reached = True

        return packets, overruns, timeout_reached

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
        if not type(size)==list: size = [size]*self.get_channels()

        # Idiot proofing
        if not len(size) == self.get_channels():
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

    def work(self, input_items, output_items):
        
        # Acquire the buffer lock for the whole function.
        self._buffer_lock.acquire()

        # Loop over channels        
        for n in range(len(input_items)):
            
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
        
        return 1

class crimson():
    """
    Interface to the Crimson. Typical workflow:
    
    c = crimson()
    
    c.enable_rx_channels([0,1])

    c.start()

    c.get_packets()
    """
    
    def __init__(self):
        """
        Interface to the Crimson. Typical workflow:
        
        c = crimson()
        c.enable_rx_channels([0,1])
        c.start()
        c.get_packets()
        """
        self._enabled_rx_channels = None
        self.buffer               = None
        self._top_block           = None
        self._crimson             = None
        
        return        
        
    def enable_rx_channels(self, channels=[0,1], buffer_size=500):
        """
        Enables the specified channels.
        
        Behind the scenes, this creates a top block and a GNU radio UHD 
        usrp_source object, then connects the usrp_source to a data buffer 
        (stored in self.buffer).
        """
        # Stop and clear any old processes
        if not self._top_block == None: 
            self._top_block.lock()            
            self._top_block.stop()
        
        # Store for safe keeping
        self._enabled_rx_channels = channels        
        
        # Create the buffer
        self.buffer = _data_buffer(500, len(channels))
        
        # Create the gnuradio top block
        self._top_block = _gr.top_block("Crimson")

        # Create the crimson data faucet        
        self._crimson = _uhd.usrp_source("crimson",
        	_uhd.stream_args(cpu_format="sc16", args='sc16', channels=(channels)))
        
        # Connect the faucet to the buffer
        for n in range(len(channels)):
            self._top_block.connect((self._crimson, n), (self.buffer, n))

    def get_enabled_rx_channels(self):
        """
        Returns the list of enabled rx channels.
        """
        return self._enabled_rx_channels
    
    def get_samples(self, N=1024, keep_all=False, timeout=1.0):
        """
        Gets the most recent N samples from each channel, plus a list of 
        overrun counts for each. Resets the overrun counts. Return format:
        
        [(X0,Y0,overruns0), (X1,Y1,overruns1), ...]
        
        N           Number of samples to get from each channel. Can be an 
                    integer or a list of integers (one for each channel).
        
        keep_all    If True, this will return all the data from the buffer,
                    not just the most recent N samples. Use this option if you
                    need to take long continuous blocks of data or want to 
                    maximize your duty cycle.
                    
        timeout     Number of seconds to wait for the N samples.
        """
        # Make sure we have a sample number for each channel
        if not type(N) == list: N = [N]*len(self.get_enabled_rx_channels())
        
        # Idiot proof
        if not len(N) == len(self.get_enabled_rx_channels()):
            print "ERROR get_samples(): size of N list not equal to number of channels."
        
        # Get the packets from the buffer
        packets, overruns, timeout_reached = self.buffer.get_packets(N, keep_all, timeout)
        
        # For each channel, get the x and y quadrature data
        data = []
        for n in range(len(packets)):

            if len(packets[n]):
                # make it the right shape            
                x,y = _n.concatenate(packets[n]).transpose()
                
                # If we're not keeping everything, throw away the old stuff
                if not keep_all: 
                    x = x[len(x)-N[n]:]
                    y = y[len(y)-N[n]:]
                
            else:
                x = _n.array([], dtype=_n.int16)
                y = _n.array([], dtype=_n.int16)

            # Store it
            data.append((x,y,overruns[n]))
        
        # Return it
        return data

    def start(self):
        """
        Start the data flow.
        """
        self._top_block.start()
    


if __name__ == '__main__':   
    self = crimson()
    self.enable_rx_channels([0,1])
    self.start()
    self.get_samples()
