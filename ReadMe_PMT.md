# PMT reconstruction integrated in CYGNO Analysis

## Running the PMT analysis code:

### In the configuration file, you will find the PMT part at the end:

1. Turn on/off PMT reconstruction mode.

    - `'pmt_mode'              : 0`

2. Describe in which channels on the digitizer are the PMTs connected to. The code assumes the PMTs are connected to the same channels in both digitizers.
    
    - `'board_pmt_channels'	: [1,2,3,4]` 

### The following options are connected to the Python library used to find the peaks in the PMT waveforms. You can find all its documentation here: [scipy.signal.find_peaks](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html )

3. _Threshold_ refers to minimum amplitude needed for the peak to be identify. Not in use.

    - `'threshold'             : 0`

4. _height\_RMS_ refers to the minimum ampltude for the peak to be identified. We use this instead of threshold.

    - `'height_RMS'			: 5,`

5. _minPeakDistance_ refers to the minimum distance in samples between two peaks. Not used.

    - `'minPeakDistance'       : 1`

6. _prominence_ refers to the minimum amplitude for the peak to be considered wrt to the neighbouring baseline. Not in use.

    - `'prominence'            : 0.1`

7. _fixed\_prom_ is a boolean that forces the prominence to be overridden with an optimized value.			

    - `'fixed_prom'			: True`

8. _width_ refers to the minimum width for the peak to be considered.

    - `'width'                 : 5`

### The following options are connected other properties of the PMT analysis such as _saving information_ and _GEM signals_

9. _resample_ refers to the number of samples used to perform a moving average action (low-pass filter). Minimum value is 1.
    
    - `'resample'              : 10`

10. _pmt\_plotpy_ lets the user save *ALL* the waveforms in '{--pdir}./waveforms' of the analysed events. Be careful, each picture has between 10-100 associated waveforms. If you run turn this variable _True_ and run on 400 events, you will tens of thousands of waveform pngs. This is useful for testing parameters and other quick tests.

    - `'pmt_plotpy'            : False`

11. _pmt\_wf\_in\_tree_ lets the user save the full raw Y array of the waveform in the corresponding tree branch. The X-axis should be manually created. This option might increase considerably the size of the output files.
    
    - `'pmt_wf_in_tree'        : False` 

12. _pmt\_verbose_ lets the user decide if and how much information to print from the waveform analysis. You can choose from *0* (no output) to *3* (full output) 

    - `'pmt_verbose'			: 0`	

13. _include\_gem_ lets the user decide if to also readout the GEM signals. Change from *0* (no) to *1* (yes). This option is only available to LNGS data at the moment. Note that the analysis is very simple and further improvements should be added into the script _waveform.py_ to properly analyse these waveforms. These are also only available in the fast digitizer

    - `'include_gem'			: 0`

14. _board\_gem\_channels_: Similarly to the PMT channels, here the user should specify to which channels of the digitizer are the GEMs connected.

    - `'board_gem_channels'	: [5,6,7]`

#### Note: Both axis of these waveforms are shown in the raw form, meaning _x_ is _samples_, and _y_ is _ADC counts_. The transformation to nanoseconds and millivolts is not automatically made.

## Digitizers relevant information:

The CYGNO DAQ system uses two digitizers, typically designated as *fast* and *slow* board/digitizer.

The intent is to be able to save signals with both small and large time extensions - fast and slow, respectively.

- Examples of *short/fast* signals are low energy (<20 keV) electron recoils (e.g. the 55Fe X-rays at 6 keV). High energy and high energy deposition particle (such as alpha particles), can also be contained in time most of the times.

- Examples of *long/slow* signals would be higher energy cosmics or electron recoils (>100 keV and/or MIPs).

### Models and settings currently used:

- **Fast digitizer:** CAEN V1742 
    - https://www.caen.it/products/v1742/
    - 12 bits
        - ADC counts are within [0,4096]
    - Dynamic range: 1Vpp
    - DC offset: Variable [-1 to +1] V
        - Chosen value for all detectors: -0.3V
            - Note: When using values different from 0, a new set of corrections should be implemented.
    - Sampling: Variable - [5 GS/s , 2.5 GS/s, 1 GS/s, 750 MS/s]
        - Chosen sampling frequency for all detectors: 750Mhz
        - X-sample duration: 1/750MS = 1.33(3) ns, for 1024 samples, thus covering 1365.3(3) ns.

- **Slow digitizer:** CAEN V1720 
    - https://www.caen.it/products/v1720/
    - 12 bits
        - ADC counts are within [0,4096]
    - Dynamic range: 2Vpp
        - (meaning it observes half the intensity in the same signal seen by the fast digitizer)
    - DC offset: [-1 to +1] V
        - No offset nor corrections are applied to this digitizer.
    - Sampling: 250Mhz
        - X-sample duration: 1/250MS = 4 ns, for 4000 samples,  thus covering 16000 ns.
            - NB: 4000 is not a typical power of 2 due to old issue with the DAQ. 
            - 4000 has been chosen to cover the longest event in the detector, which assuming an electron drift velocity of 5cm/us, is about 10 us.

#### Example: [Fast and Slow digitizers](https://imgur.com/a/g5qyqyF)

### For mechanical details and disposition of the PMTs, please check the [CYGNO official wiki git repository](https://github.com/CYGNUS-RD/WIKI-documentation/wiki/Detector-General)

### Possible updates upon request:

- Automatic conversion of GEM signal amplitude to Coulomb. Requires the full electrical scheme of the detector.
- Automatic conversion of the axis in the waveforms. Not yet implemented as it is preferable to retrieve the raw data and apply the conversion in post-processing.