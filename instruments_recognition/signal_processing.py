"""
signal_processing.py
"""

import scipy
from scipy import signal
import numpy as np
import glob
import wave

#for testing
import matplotlib.pyplot as plt 


class FSignal:
    """Represents an array of floats together with a samplerate"""
    #fsignal
    #samplerate

    def __init__(self, floats, samplerate):
        self.fsignal = np.array(floats)
        self.samplerate = samplerate

    @classmethod
    def from_wav_file(cls, wav_file):
        """
        input: a .wav file
        output: an fsignal
        """
        wf = wave.open(wav_file)
        samplerate = wf.getframerate()
        bitdepth = wf.getsampwidth() * 8
        # print("the depth is: " + str(bitdepth))
        channels = wf.getnchannels()
        # print("the number of channels is: " + str(channels))
        length_sample = wf.getnframes()
        # print("the number of samples is: " + str(length_sample))
        bitdepthnp = np.int8 if bitdepth == 8 else np.int16
        imported = np.fromstring(wf.readframes(length_sample*channels), bitdepthnp)
        wf.close()
        if channels == 1:
            imported1ch = imported
        else :
            # keep only left channel
            imported1ch = imported[::2]
        # normalize
        tofloats = imported1ch/(2**float(bitdepth))
        return cls(tofloats, samplerate)

    # create from a wave file and: remove possible silence, normalize, and apply window
    @classmethod
    def from_wav_file_and_clean(cls, wav_file):
        """
        input: a .wav file
        output: an fsignal with silences removed, normalized, and windowed (fade in, fade out)
        """
        a_signal = cls.from_wav_file(wav_file)
        a_signal.kill_silence()
        a_signal.normalize_signal()
        a_signal.window()
        return a_signal

    def kill_silence(self):
        """
        effect: kills the silence at the begining and end of signal
        """
        # TODO: justify this, pretty hardcoded
        threshold = 0.000458

        first = 0
        last = len(self.fsignal)-1
        while (abs(self.fsignal[first]) < threshold) :
            first += 1
        while (abs(self.fsignal[last]) < threshold) :
            last -= 1
        if first >= last :
            raise ValueError('This signal seems to be pure noise!')
        self.fsignal = self.fsignal[first:last]

    def window(self):
        """
        effect: applies a window (fade in, fade out)
        """
        floats = self.fsignal
        #parameter for tukey window
        # apply tukey window
        self.fsignal = signal.tukey(len(floats),alpha=0.05) * floats
 
    def normalize_signal(self):
        """
        effect: normalizes the signal so that the maximum is 1
        """
        self.fsignal = self.fsignal/max(abs(self.fsignal))
    
    def many_windows(self, howmany, duration):
        """
        returns a list of signals
        """
        nsamples_over2 = int(duration * self.samplerate/2)
        nsamples = nsamples_over2 * 2
        # pad the signal so it doesn't fail if it is too short
        totalsignal = np.pad(self.fsignal,[(0,howmany*nsamples)],'constant')
        return [ FSignal(totalsignal[n*nsamples_over2:nsamples+n*nsamples_over2], self.samplerate)
                 for n in range(0,howmany) ]

    def pad_a_little(self, sec):
        """
        effect: pads the signal sec seconds at the beggining and end
        """
        howmanysamples = int(sec*self.samplerate)
        self.fsignal = np.pad(self.fsignal,[(howmanysamples,howmanysamples)],'constant')

    def unpad_a_little(self, sec):
        """
        effect: reverts the effect of pad_a_little
        """
        howmanysamples = int(sec*self.samplerate)
        self.fsignal = self.fsignal[howmanysamples:len(self.fsignal)-howmanysamples]


    #def copy(self):
        #TODO


    def local_volume(self):
        """
        output: an fsignals containing the "local volume" in each moment of self
        side effect: pads the signal
        """
        # the "localness" is 30 ms
        window_len = 0.03
        # pad so we can inegrate easily (this is reverted later)
        self.pad_a_little(window_len)
        num_samples_window = self.samplerate * window_len
        n_samp_over2 = int(num_samples_window/2)
        squared_signal = self.fsignal**2


        len_paded = len(squared_signal)
        window_len_samples = window_len * self.samplerate
        squared_signal_integral = np.empty(len_paded) 
        squared_signal_integral[0] = squared_signal[0]
        for pos in range(1,len_paded-n_samp_over2) :
            squared_signal_integral[pos] = squared_signal_integral[pos-1] + squared_signal[pos]
        

#        plt.plot(squared_signal_integral)
#        plt.show()

        local_volume = squared_signal_integral[n_samp_over2*2:len_paded-n_samp_over2*2] - squared_signal_integral[:len_paded-n_samp_over2*4]

#        print("how many samples in the window: " + str(n_samp_over2))
#        print(squared_signal_integral[:len_paded-n_samp_over2*2])
#        print(squared_signal_integral[n_samp_over2*2:len_paded-1])
#        print(squared_signal_integral[n_samp_over2*2:n_samp_over2*4])

#        #plt.plot(local_volume)
#        plt.plot(squared_signal_integral[::len_paded-n_samp_over2*2])
#        plt.show()
#        plt.plot(squared_signal_integral[n_samp_over2*2::])
#        plt.show()
        # leave the input sample as we find it
        self.unpad_a_little(window_len)

        res = FSignal(local_volume,self.samplerate)
        #res.unpad_a_little(window_len)

        #local_volume = [ squared_signal_integral[pos+n_samp_over2] - squared_signal_integral[pos-n_samp_over2]
        #               for pos in range(0,len(squared_signal_integral)) ]

        ##        integral_from_zero = [ 
        return res

    #def only_from_to(sec_start, sec_end):
        #TODO

#    def auto_crop(self):
#        """
#        output: a list of fsignals
#        """

            


relevant_range_min = 100.
relevant_range_max = 15000.
howmuchtopad = 1000

class FSpectrum:
    #spectrum
    #freq_lbls
    #samplerate

    def __init__(self, fsignal):
        #window the signal
        fsignal.window()

        floats = fsignal.fsignal
        samplerate = fsignal.samplerate
        self.samplerate = samplerate

        # pad with some zeros
        padded = np.pad(floats,[(howmuchtopad,howmuchtopad)],'constant')
        # transform to frequency domain (since this is a real signal the transformed vector has half of the length)
        transformed_complex = np.fft.rfft(padded)
        # take absolute values
        self.spectrum = list(map(abs,transformed_complex))
        # compute frequencies for the xlabel
        self.freq_lbls = np.fft.rfftfreq(len(padded),1/samplerate)
        
        # TODO: Plancherel's theorem does not seem to be holding :(
        #print "norm of signal: " + str(rms(ex_windowed))
        #print "norm of transformed: " + str(rms(ex_transformed))

    def find_tonic_pos(self, relevant_range_min, relevant_range_max):
        """
        input: a transformed signal
        output: the position in the array of the tonic (not the frequency of the tonic!)
        """
        spec = self.spectrum
        # the total spectrum has len(spec) samples and the spectrum has 44100/2 frequencies
        max_frequency = self.samplerate/2
        specinitial = int(len(spec)/max_frequency * relevant_range_min)
        specfinal = int(len(spec)/max_frequency * relevant_range_max)
            # print("the minum relevant frequency is: " + str(frq_lbl[specinitial]))
            # print("this is the actual maximum: " + str(frq_lbl[len(spec)-1]))
            # print("this is the length of the specample: " + str(len(spec)))
        # find tonic
        tonic_pos = np.argmax(spec[specinitial:specfinal]) + specinitial
        #print("the tonic is in: " + str(tonic_pos))
        #print("the volume is: " + str(spec[tonic_pos]))
        return tonic_pos

    def find_tonic(self, relevant_range_min, relevant_range_max):
        pos = self.find_tonic_pos(relevant_range_min, relevant_range_max)

        spec = self.spectrum
        max_frequency = self.samplerate/2
        tonic_freq = pos/len(spec) * max_frequency
        # print("the tonic is: " + str(tonic_freq) + " Hz")
        return tonic_freq

    def frequency_energy(self, frequency):
        max_frequency = self.samplerate/2
        freq_pos = len(self.spectrum) * frequency/max_frequency
        return sum(interval_to_integrate(self.spectrum,freq_pos))



num_harm = 10
class Harmonics:
    """Represents a vector of volumes of harmonics"""
    #volumes
    #frequencies

    def __init__(self, a_spectrum, a_tonic, num_harm):
       self.frequencies = [ a_tonic * n for n in range(1,num_harm) ]
       volumes_abs = [ a_spectrum.frequency_energy(frq) for frq in self.frequencies ]
       volume_tonic = volumes_abs[0]
       self.volumes = volumes_abs/volume_tonic


howmanywindows = 4
windowsduration = 0.12
def signal_to_harmonics(a_signal):
    """
    input: a signal
    output: a list of harmonics
    """
    a_spectrum = FSpectrum(a_signal)
    tonic_freq = a_spectrum.find_tonic(relevant_range_min, relevant_range_max)

    windows = a_signal.many_windows(howmanywindows, windowsduration)

    the_volumes = [ rms(win.fsignal) for win in windows ]
    the_volumes_rel = the_volumes/the_volumes[0]

    the_harmonics = [ Harmonics(FSpectrum(win), tonic_freq, num_harm) for win in windows ]
    return the_harmonics#, the_volumes_rel


def signal_to_volumes(a_signal):
    """
    input: a signal
    output: a list of volumes (of harmonics), ommiting the firstone, which is 1
    """
    return [ har.volumes[1:] for har in signal_to_harmonics(a_signal) ]


### AUX functions

# interval of 200Hz centered at the note
width_of_average = 200.
def interval_to_integrate(ex, peak) :
    # the total spectrum has len(ex) samples and the spectrum has 44100/2 frequencies
    max_frequency = samplerate/2
    hz_in_samples = len(ex)/max_frequency * width_of_average
    windowinitial = int(peak-hz_in_samples/2) 
    windowfinal = int(peak+hz_in_samples/2)
    #print "the length in samples of the average window is: " + str(windowfinal-windowinitial)
    return ex[windowinitial:windowfinal]


# root mean square of a vector
def rms(vec) :
    return np.sqrt(np.average(np.power(vec,2)))

def sumofsquares(vec) :
    return sum(np.power(vec,2))

# l^2 norm of a vector
def l2n(vec) :
    return np.sqrt(sumofsquares(vec))


### IMPORTING synthetic sinusoids (for testing)

samplerate = 44100.
Ts = 1.0/samplerate

def gimme_sinusoid(freq) :
    t = np.arange(0,1,Ts)
    return np.sin(2*np.pi * freq * t)
    
def gimme_sinusoids(freq,harmonics_amp) :
    n = 2
    res = gimme_sinusoid(freq)
    for amp in harmonics_amp : 
        res += amp * gimme_sinusoid(n*freq)
        n += 1
    return res

def gimme_sinusoids_noise(freq,harmonics_amp) :
    signal = gimme_sinusoids(freq,harmonics_amp)
    return signal + np.random.normal(0,.5,len(signal))
