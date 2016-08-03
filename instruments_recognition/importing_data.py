"""
importing_data.py
"""

import scipy
from scipy import signal
from scipy.io import wavfile
import numpy as np
import glob
import matplotlib.pyplot as plt

import wave


class FSignal:
    """Represents an array of floats with a samplerate"""
    #fsignal
    #samplerate

    def __init__(self, floats, samplerate):
        self.fsignal = np.array(floats)
        self.samplerate = samplerate

    @classmethod
    def from_wav_file(cls, wav_file):
        """
        input: a .wav file
        """
        wf = wave.open(wav_file)
        samplerate = wf.getframerate()
        bitdepth = wf.getsampwidth() * 8
        print("the depth is: " + str(bitdepth))
        channels = wf.getnchannels()
        print("the number of channels is: " + str(channels))
        length_sample = wf.getnframes()
        print("the number of samples is: " + str(length_sample))
        bitdepthnp = np.int8 if bitdepth == 8 else np.int16
        imported = np.fromstring(wf.readframes(length_sample*channels), bitdepthnp)
        wf.close()
        if channels == 1:
            imported1ch = imported
        else :
            # keep only left channel
            imported1ch = [n for n,m in imported]
        # normalize
        tofloats = [ sample/(2**float(bitdepth)) for sample in imported1ch]
        return cls(tofloats, samplerate)

    @classmethod
    def from_wav_file_and_clean(cls, wav_file):
        a_signal = cls.from_wav_file(wav_file)
        a_signal.kill_silence()
        a_signal.normalize_signal()
        a_signal.window()
        return a_signal

    def kill_silence(self):
        """
        kills the silence at the begining and end of signal
        """
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
        floats = self.fsignal
        #parameter for tukey window
        alpha_tukey = 0.05
        # apply tukey window
        self.fsignal = signal.tukey(len(floats),alpha=0.05) * floats
 
    def normalize_signal(self):
        """
        normalizes the signal so that the maximum is 1
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


relevant_range_min = 100.
relevant_range_max = 15000.
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

        howmuchtopad = 1000
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
        print("the tonic is: " + str(tonic_freq) + " Hz")
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



# interval of 200Hz centered at the note
def interval_to_integrate(ex, peak) :
    width_of_average = 200.
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









########################################## not used anymore

### IMPORTING DATA functions

# samplerate of the data (samples/sec.)
samplerate = 44100.

# bit depth of the data
bitdepth = np.int16

# attack (in sec.)
attack_duration = 0.4

# when a wav audio is noise
# TODO: very arbitrary
threshold_wav = 30

def import16monowav(wav_file) :
    """
    input: a .wav file 
    output: a vector of integers (one float per sample)
    """
    samplerate, data = wavfile.read(wav_file)
    #return np.fromfile(open(wav_file),bitdepth)[24:]
    return data

def import16stereo1wav(wav_file) :
    samplerate, data = wavfile.read(wav_file)
    #return np.fromfile(open(wav_file),bitdepth)[24:]
    res = [n for n,m in data]
    return res

def normalize16bits(vector) :
    """
    normalize to a float in [-1,1]
    """
    return [ sample/2**16. for sample in vector ]

def killsilence_wav(wav_16bits) :
    """
    kills first and last samples under some threshold, of a integer array
    """
    first = 0
    last = len(wav_16bits)-1
    while (abs(wav_16bits[first]) < threshold_wav) :
        first += 1
    while (abs(wav_16bits[last]) < threshold_wav) :
        last -= 1
    if first >= last :
        raise ValueError('This audio seems to be pure noise!')
    return wav_16bits[first:last]


# TODO: the silence should be removed later?
def convert_to_float(ex,channels) :
    """
    input: a wav file
    output: a vector of floats in [-1,1]
    """
    if channels == 1 :
        imported = import16monowav(ex)
    else :
        imported = import16stereo1wav(ex)
    return normalize16bits(killsilence_wav(imported))
        
    
# import wav, convert to float, take pieces and fourier transform each piece
def import_convert_transform2(ex,channels) :
    """
    input: a wav file
    output: 1) a list of transforms of windows of the audio
            2) the relative volumes wrt the attack (first window)
    """
    # import wav and convert to float
    ex_float = convert_to_float(ex,channels)
    windows = attack_and_release(ex_float)
    volumeswindows = list(map(rms,windows))
    attackvolume = volumeswindows[0]
    volumeswindows_relative = volumeswindows/attackvolume
    return list(map(transform_floats, windows)), volumeswindows_relative[1:]

 
# de comun
def import_convert_transform_release(ex,channels) :
    # import wav and convert to float
    ex_float = convert_to_float(ex,channels)
    windows = attack_and_release(ex_float)
    return transform_floats(windows[1])

def import_convert(ex,channels) :
    # import wav and convert to float
    ex_float = convert_to_float(ex,channels)
    #windows = attack_and_release(ex_float)
    return ex_float

def import_convert_release(ex,channels) :
    # import wav and convert to float
    ex_float = convert_to_float(ex,channels)
    windows = attack_and_release(ex_float)
    return windows[1]



def import_convert_transform_attack(ex,channels) :
    # import wav and convert to float
    ex_float = convert_to_float(ex,channels)
    windows = attack_and_release(ex_float)
    return transform_floats(windows[0])



def attack_and_release(vec) :
    """
    given a float vector (signal) divide it in attack and release ("the rest")
    """
    attack_samples_over2 = int(samplerate * attack_duration/2)
#    print "the number of samples in the attack is: " + str(2*attack_samples_over2)
#    print "the number of samples in the release is: " + str(4*attack_samples_over2)
    attack = vec[: 2*attack_samples_over2] 
    release_possiblytooshort = vec[attack_samples_over2 : 5*attack_samples_over2]
    tooshort = 4*attack_samples_over2 - len(release_possiblytooshort)
    release = np.pad(release_possiblytooshort,[(0,tooshort)],'constant')
    return [attack, release]
    

#parameter for tukey window
alpha_tukey = 0.05

howmuchtopad = 1000

def transform_floats(ex_float) :
    """
    given a signal of floats, apply window and transform
    """
    # apply tukey window
    ex_windowed = signal.tukey(len(ex_float),alpha=0.05) * ex_float  
    # pad with some zeros
    ex_padded = np.pad(ex_windowed,[(howmuchtopad,howmuchtopad)],'constant')
    # transform to frequency domain (since this is a real signal the transformed vector has half of the length)
    ex_transformed_complex = np.fft.rfft(ex_padded)
    # take absolute values
    ex_transformed = list(map(abs,ex_transformed_complex))
    # compute frequencies for the xlabel
    freq_label = np.fft.rfftfreq(len(ex_padded),1/samplerate)
    
    # TODO: Plancherel's theorem does not seem to be holding :(
    #print "norm of signal: " + str(rms(ex_windowed))
    #print "norm of transformed: " + str(rms(ex_transformed))
    
    return ex_transformed, freq_label
    #freq_label = [ np.arange(len(ex))/(len(ex)*2/samplerate) for ex in ex_transformed ]





### FINIDING HARMONICS functions


# TODO: this should not be sooo hardcoded 
peak_width_min = 10
peak_width_max = 20

# a semi tone
sm_tone = 2**(1/12.)
tone = 2**(1/6.)


# range of harmonics of interest
# TODO: hardcode this better
relevant_range_min = 100.
relevant_range_max = 15000.
# how many harmonics we want to take into account
num_harmonics = 9


# TODO: revise the use of ex = list(ex), doesn't make much sense to do this everywhere...

# this is just for backwards compatibility
def harmonics_energy_compatibility(ex, frq_lbl) :
    """
    just for now, delete this function later
    """
    ex = list(ex)
    tonic_pos = find_tonic(ex)
    return harmonics_energy(ex,frq_lbl,tonic_pos)


def harmonics_energy(ex, frq_lbl, tonic_pos) :
    """
    input: 1) a transformed signal 
           2) a vector with the corresponding frequencies
           3) the position of the tonic in the transformed signal
    output: 1) a vector with the frequency of the harmonics
            2) a vector with the volume of each harmonic (with respect to the volume of the tonic)
    """
    ex = list(ex)
    # compute the energy of the harmonics
    peaks_energy = [ sum(interval_to_integrate(ex,tonic_pos*n)) for n in range(1,num_harmonics+1) ]
    # and normalize with the energy of the tonic
    tonic_energy = peaks_energy[0]
    peaks_energy_nor = [ harmonic_energy/tonic_energy for harmonic_energy in peaks_energy ]
    #the frequency of the peaks
    peaks_frequency = [ (frq_lbl[tonic_pos*n] if tonic_pos*n<len(frq_lbl) else 0.) for n in range(1,num_harmonics+1) ]
    return peaks_frequency, peaks_energy_nor


def find_tonic(spec) :
    """
    input: a transformed signal
    output: the position in the array of the tonic (not the frequency of the tonic!)
    """
    spec = list(spec)
    # the total spectrum has len(spec) samples and the spectrum has 44100/2 frequencies
    max_frequency = samplerate/2
    specinitial = int(len(spec)/max_frequency * relevant_range_min)
    specfinal = int(len(spec)/max_frequency * relevant_range_max)
        # print("the minum relevant frequency is: " + str(frq_lbl[specinitial]))
        # print("this is the actual maximum: " + str(frq_lbl[len(spec)-1]))
        # print("this is the length of the specample: " + str(len(spec)))
    # find tonic
    tonic_pos = np.argmax(spec[specinitial:specfinal]) + specinitial
    print("the tonic is in: " + str(tonic_pos))
    #print("the volume is: " + str(spec[tonic_pos]))
    return tonic_pos

   
def harmonics_energy_multiwindow(exs, frq_lbls) :
    """
    input: 1) a list of transformed signals (for example: [attack, release]) we assume
              that the second vector is the release and compute the tonic based on that vector
           2) a list of the corresponding frequency for each transformed signal
    output: 1) a list of pairs (frequencies of harmonics, volume of harmonics wrt tonic)
            2) a frequency indicating the tonic of the whole signal
    """
    release_num = 1
    release = exs[release_num]
    release_frequencies = frq_lbls[release_num]
    # find tonic
    tonic_pos = find_tonic(release)
    # compute the harmonics of each window
    peaks_freqs_energies = [ harmonics_energy(ex,frqs,tonic_pos) for ex, frqs in zip(exs, frq_lbls) ]
    # the frequency of the tonic
    tonic_freq = peaks_freqs_energies[release_num][0][0]
    print("the tonic is: " + str(tonic_freq))
    return peaks_freqs_energies, tonic_freq

def harmonics_energy_multiwindow_zipped(exsfrq_lbls) :
    exs, frq_lbls = list(zip(*exsfrq_lbls))
    release_num = 1
    release = exs[release_num]
    release_frequencies = frq_lbls[release_num]
    # find tonic
    tonic_pos = find_tonic(release)
    # compute the harmonics of each window
    peaks_freqs_energies = [ harmonics_energy(ex,frqs,tonic_pos) for ex, frqs in zip(exs, frq_lbls) ]
    # the frequency of the tonic
    tonic_freq = peaks_freqs_energies[release_num][0][0]
    print("the tonic is: " + str(tonic_freq))
    return peaks_freqs_energies, tonic_freq
    
    
## NOT USED ---------------------------------


# TODO: should be called harmonics_amplitudes
def harmonics_frequency_volumes(ex, frq_lbl) :
    # filter peaks in range
    peaks_position = find_peaks_in_range(ex, frq_lbl)
    
    # compute energy of peaks (root mean square):
    peaks_energy_tmp = [ get_l2n_peak_Hz(ex,peak) for peak in peaks_position ]
    
    # find tonic
    tonic_pos = find_tonic(peaks_energy_tmp)
    #take only num_harmonics peaks
    peaks_energy = peaks_energy_tmp[tonic_pos:tonic_pos+num_harmonics]
        #print "the energy of the peaks: " + str(peaks_energy)
    # calculate volume of peaks relative to tonic
    tonic_energy = peaks_energy[0]
    peaks_energy_nor = [ harmonic_energy/tonic_energy for harmonic_energy in peaks_energy ]
        #print "the energy of the peaks normalized: " + str(peaks_energy_nor)
        ## normalize by rms
        #peaks_volume_nor = peaks_volume/rms(peaks_volume)
        
        ## put volumes in range
        ##peaks_volume_in_range = [ v-min_volume if v>min_volume else 0 for v in peaks_volume ]
    #the frequency of the peaks
    peaks_frequency = [ frq_lbl[n] for n in peaks_position[tonic_pos:tonic_pos+num_harmonics] ]
    return peaks_frequency, peaks_energy_nor

# given a signal return a list of the positions of the peaks
def find_peaks_in_range(ex, frq_lbl) :
    # find peaks (this is the position in the array)
    peaks = signal.find_peaks_cwt(ex, np.arange(peak_width_min,peak_width_max))
    peaks_in_range = filter(lambda n : frq_lbl[n]>=relevant_range_min and frq_lbl[n]<= relevant_range_max, peaks)
    return peaks_in_range


##def find_tonic2(amplitudes) :
##    maximum = max(amplitudes)
##    arbitrary = 2.
##    first = 0
##    while (amplitudes[first] < maximum/arbitrary) :
##        first += 1
##    if first >= len(amplitudes) :
##        raise ValueError('No tonic!')
##    return first


# TODO: see if this parameter is ok
tolerance_har = 0.2


def cast_harmonics(peaks_frequency, peaks_volume) :
    har_f_har_v = zip(peaks_frequency, peaks_volume)
    tonic_freq = har_f_har_v[0][0]
    peaks_volume_res = np.zeros(num_harmonics)
    peaks_frequency_new = [ tonic_freq * n for n in range(1,num_harmonics+1)]
    
    for freq, vol in har_f_har_v :
        harm_num = int(round(freq/tonic_freq))
        print("found a harmonic of number " + str(harm_num))
        print("its real frequency is " + str(freq))
        if abs(freq/tonic_freq - harm_num) > tolerance_har or harm_num > num_harmonics :
            print("i killed it :(")
            print("WARNING: strange harmonic: " + str(freq/tonic_freq))
        else :
            print("i kept it")
            peaks_volume_res[harm_num-1] = vol
            peaks_frequency_new[harm_num-1] = freq
        
    return peaks_frequency_new, peaks_volume_res
 

# rms over two tones centered at the note
def get_rms_peak_log(ex, peak) :
    windowinitial = int(peak/tone) 
    windowfinal = int(peak*tone)
    return rms(ex[windowinitial:windowfinal])

# l2n over two tones centered at the note
def get_l2n_peak_log(ex, peak) :
    windowinitial = int(peak/tone) 
    windowfinal = int(peak*tone)
    return l2n(ex[windowinitial:windowfinal])
 
# average over two tones centered at the note
def get_amplitude_peak_log(ex, peak) :
    windowinitial = int(peak/tone) 
    windowfinal = int(peak*tone)
    return np.average(ex[windowinitial:windowfinal])

# rms over a neigborhood of the peak
def get_rms_peak(ex, peak) :
    # TODO: this 10 is pretty arbitrary
    windowinitial = int(peak-peak_width_max/2) 
    windowfinal = int(peak+peak_width_max/2)
    return rms(ex[windowinitial:windowfinal])

# with a change of parameter
# TODO: do the integral with np.average
def get_amplitude_peak_log_log(ex, peak) :
    # v_n(x) = v_f(2^x),  thus \int v_n(x) = \int v_f(2^x) = \int v_f(y)/y
    windowinitial = int(peak/sm_tone) 
    windowfinal = int(peak*sm_tone)
    lenintegral = windowfinal - windowinitial
    weights = [ 1/float(y) for y in range(windowinitial,windowfinal) ]
    totalweights = sum(weights)
    averagevol = sum([v*w for v,w in zip(ex[windowinitial:windowfinal], weights)])/totalweights
    return averagevol
