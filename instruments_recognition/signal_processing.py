"""
signal_processing.py
"""

import scipy
from scipy import signal
import numpy as np
import glob
import wave
import re
import os

#for testing
import matplotlib.pyplot as plt 


class FSignal:
    """Represents an array of floats together with a samplerate"""
    #fsignal
    #samplerate
    #name

    def __init__(self, floats, samplerate, name):
        self.fsignal = np.array(floats)
        self.samplerate = samplerate
        self.name = name

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
        tofloats = imported1ch/(2**(float(bitdepth)-1))

        name = re.sub(r'.wav', r'', os.path.basename(wav_file))
        
        return cls(tofloats, samplerate, name)

    # create from a wave file and: remove possible silence, normalize, and apply window
    @classmethod
    def from_wav_file_and_clean(cls, wav_file) :
        """
        input: a .wav file
        output: an fsignal with silences removed, normalized, and windowed (fade in, fade out)
        """
        a_signal = cls.from_wav_file(wav_file)
        a_signal.kill_silence()
        a_signal.normalize_signal()
        a_signal.simpl_window()
        return a_signal

    def to_wav_file(self, directory) :
        """
        input: a directory (must have a slash at the end!)
        effect: creates a mono 16bit .wav file with the signal
        """
        bitdepth = 16.

        out_file = wave.open(directory + self.name + ".wav", 'wb')

        # mono
        out_file.setnchannels(1)
        # samplerate
        out_file.setframerate(self.samplerate)
        # 16 bit
        out_file.setsampwidth(2)

        # this cast to integer is not the best thing in the world
        # should apply dithering...
        rawdata_int = (self.fsignal * (2**(float(bitdepth)-1))).astype(np.int16)
        rawdata_byte = rawdata_int.tostring()
        out_file.writeframes(rawdata_byte)

        out_file.close()


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

    def window(self,alpha):
        """
        effect: applies a window (fade in, fade out)
        """
        floats = self.fsignal
        window = signal.tukey(len(floats),alpha = alpha)
        # apply gaussian window
        self.fsignal = window * floats
 
    def simpl_window(self):
        floats = self.fsignal
        # parameter for tukey window
        alpha_t = 0.05
        # apply tukey window
        self.fsignal = signal.tukey(len(floats),alpha=alpha_t) * floats
 

    def normalize_signal(self):
        """
        effect: normalizes the signal so that the maximum is 1
        """
        self.fsignal = self.fsignal/max(abs(self.fsignal))
    
    def many_windows(self, howmany, duration):
        """
        returns a list of signals
        """
        orig_name = self.name

        nsamples_over2 = int(duration * self.samplerate/2)
        nsamples = nsamples_over2 * 2
        # pad the signal so it doesn't fail if it is too short
        totalsignal = np.pad(self.fsignal,[(0,howmany*nsamples)],'constant')
        return [ FSignal(totalsignal[n*nsamples_over2:nsamples+n*nsamples_over2], self.samplerate, orig_name + "_window_" + str(n))
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

        # leave the input sample as we find it
        self.unpad_a_little(window_len)

        res = FSignal(local_volume,self.samplerate, self.name + "_local_volume")
        return res

    def from_to_samples(self, sample_start, sample_end):
        """
        out: the cropped sample from sample sample_start to sample sample_end
        """
        a_signal = self.fsignal[sample_start:sample_end]
        return FSignal(a_signal, self.samplerate, self.name + "_crop_samp" + str(sample_start) + "to" + str(sample_end))

    def from_to_sec(self, sec_start, sec_end):
        """
        out: the cropped sample from second sec_start to second sec_end
        """
        sample_start = self.samplerate * sec_start
        sample_end = self.samplerate * sec_end
        return self.from_to_samples(sample_start, sample_end)

    def find_beats(self, threshold):
        """
        output: a list of pairs (bit_start, bit_end), where bit_start and bit_end
                are the positions in the sample (in samples, *not* seconds)
        """
        local_volume_signal = self.local_volume().fsignal
        N = len(local_volume_signal)

        # attack = 30ms
        attack = 0.03
        samples_attack = int(attack * self.samplerate)
        # hold = 15ms
        hold = 0.015
        samples_hold = int(hold * self.samplerate)
        # release = 30ms
        release = 0.03
        samples_release = int(release * self.samplerate)

        beats = []
        pos = 0
        while pos < N-1 :
            if local_volume_signal[pos] > threshold :
                beat_start = max(0, pos - samples_attack)

                #print("found beat at pos: " + str(pos))

                while pos < N-1 and np.max(local_volume_signal[pos:min(pos + samples_hold, N-1)]) > threshold :
                    # hold 15 ms
                    pos = min(pos + samples_hold, N-1)

                # release
                pos = min(pos + samples_release, N-1)

                beat_end = pos
                beats.append([beat_start, beat_end])
            else :
                pos += 1
        return beats


    def auto_crop(self, threshold):
        """
        output: a list of fsignals
        """
        beats_beg_end = self.find_beats(threshold)

        # attack = 30ms
        attack = 0.030
        samples_attack = int(attack * self.samplerate)

        # release = 30ms
        release = 0.030
        samples_release = int(release * self.samplerate)


        res = [ self.from_to_samples(max(0,b-samples_attack),min(e+samples_release, len(self.fsignal)-1)) for b,e in beats_beg_end ]
        # window
        [ r.window(0.3) for r in res ]
        # pad
        [ r.pad_a_little(0.1) for r in res ]
        return res



            


relevant_range_min = 100.
relevant_range_max = 15000.
howmuchtopad = 1000

class FSpectrum:
    #spectrum
    #freq_lbls
    #samplerate

    def __init__(self, fsignal):
        #window the signal
        fsignal.simpl_window()

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