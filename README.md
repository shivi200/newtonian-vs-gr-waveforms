# newtonian-vs-gr-waveforms
Comparing Newtonian and GR peak merger frequencies for BH-NS binaries using SPH data and PyCBC SEOBNRv4 - FFT analysis
"""
Newtonian vs General Relativity: Gravitational Waveform Comparison
===================================================================
Author: Shivika Lamba
Programme: BRICS Astronomy & IDIA Data Analytics Training Programme
 
Compares gravitational wave peak frequencies from black hole–neutron star
mergers modelled under two frameworks:
 
  1. Newtonian framework  — SPH simulation data from Lee & Kluzniak (1998-2000)
  2. GR framework         — PyCBC SEOBNRv4 waveform approximant
 
Key finding: GR predicts significantly higher peak merger frequencies,
confirming that relativistic effects are essential near merger.
 
Physical parameters: M_BH = 4.516 M_sun, M_NS = 1.4 M_sun
 
Requirements:
  pycbc, gwpy, lalsuite, numpy, scipy, matplotlib, seaborn
"""
 
# =============================================================================
# 1. IMPORTS
# =============================================================================
 
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from pycbc.waveform import get_td_waveform
 
# =============================================================================
# 2. LOAD NEWTONIAN SPH WAVEFORM DATA
# =============================================================================
# Data from Lee & Kluzniak (1998, 1999) and Lee (2000)
# Five configurations varying equation of state (Gamma) and spin state
#
# File format: columns = [time (s), r*h_plus (cm), r*h_cross (cm)]
# where r*h is the waveform amplitude scaled by source distance
 
DATA_FILES = [
    'data/waveform1BHNS.out',   # Gamma = 3.0, Irrotational
    'data/waveform2BHNS.out',   # Gamma = 2.0, Irrotational
    'data/waveform3BHNS.out',   # Gamma = 5/3, Irrotational
    'data/waveform4BHNS.out',   # Gamma = 3.0, Tidally Locked
    'data/waveform5BHNS.out',   # Gamma = 5/3, Tidally Locked
]
 
LABELS = [
    "Gamma=3.0, Irrotational",
    "Gamma=2.0, Irrotational",
    "Gamma=5/3, Irrotational",
    "Gamma=3.0, Tidally Locked",
    "Gamma=5/3, Tidally Locked",
]
 
# Load all five waveform files
waveforms = []
for path, label in zip(DATA_FILES, LABELS):
    data = np.loadtxt(path, skiprows=1)
    t        = data[:, 0]   # time in seconds
    h_plus   = data[:, 1]   # plus polarisation (amplitude x distance)
    h_cross  = data[:, 2]   # cross polarisation
    waveforms.append((t, h_plus, h_cross))
    print(f"Loaded {label}: {len(t)} time steps, "
          f"t = [{t[0]:.4f}, {t[-1]:.4f}] s")
 
# =============================================================================
# 3. PLOT NEWTONIAN WAVEFORMS
# =============================================================================
 
fig, axes = plt.subplots(5, 1, figsize=(14, 18))
fig.suptitle("Newtonian BH–NS Waveforms: h+ and h× Polarisations", fontsize=14)
 
for i, (ax, (t, h_plus, h_cross), label) in enumerate(
        zip(axes, waveforms, LABELS)):
    ax.plot(t, h_plus,  label='h+ (Plus)',  color='steelblue')
    ax.plot(t, h_cross, label='h× (Cross)', color='firebrick', alpha=0.8)
    ax.set_title(label, fontsize=11)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude × Distance (cm)')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.4)
 
plt.tight_layout()
plt.savefig('newtonian_waveforms.png', dpi=150)
plt.show()
 
# =============================================================================
# 4. FFT ANALYSIS — EXTRACT PEAK FREQUENCIES
# =============================================================================
 
def compute_fft(time, signal):
    """
    Compute the amplitude spectrum of a gravitational wave signal using FFT.
 
    A Hanning window is applied before the transform to reduce spectral
    leakage caused by signal discontinuities at the edges.
 
    Parameters
    ----------
    time   : array-like, time values in seconds
    signal : array-like, GW strain or amplitude (h+ or h×)
 
    Returns
    -------
    positive_freqs     : array of positive frequencies (Hz)
    amplitude_spectrum : corresponding FFT amplitude values
    peak_freq          : dominant frequency (Hz)
    """
    N = min(len(time), len(signal))
    time   = time[:N]
    signal = signal[:N]
 
    dt = time[1] - time[0]                        # sampling interval
    freqs      = fftfreq(N, dt)                   # frequency array
    fft_result = fft(signal * np.hanning(N))      # windowed FFT
 
    # Retain only positive frequencies (spectrum is symmetric for real signal)
    positive_freqs     = freqs[:N // 2]
    amplitude_spectrum = np.abs(fft_result[:N // 2])
 
    peak_idx  = np.argmax(amplitude_spectrum)
    peak_freq = positive_freqs[peak_idx]
 
    return positive_freqs, amplitude_spectrum, peak_freq
 
 
# Compute and plot FFT for each Newtonian configuration
print("\n--- Newtonian Peak Frequencies ---")
newtonian_peak_freqs = []
 
fig, axes = plt.subplots(5, 1, figsize=(12, 18))
fig.suptitle("FFT of Newtonian h+ Waveforms", fontsize=14)
 
for i, ((t, h_plus, _), label, ax) in enumerate(
        zip(waveforms, LABELS, axes)):
    freqs, amps, peak_freq = compute_fft(t, h_plus)
    newtonian_peak_freqs.append(peak_freq)
 
    ax.plot(freqs, amps, color='steelblue')
    ax.axvline(peak_freq, color='red', linestyle='--',
               label=f'Peak: {peak_freq:.2f} Hz')
    ax.set_title(label, fontsize=11)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Amplitude (a.u.)')
    ax.set_xlim(100, 150)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)
    print(f"  {label}: {peak_freq:.2f} Hz")
 
plt.tight_layout()
plt.savefig('newtonian_fft.png', dpi=150)
plt.show()
 
# =============================================================================
# 5. GR WAVEFORM — PyCBC SEOBNRv4
# =============================================================================
# Generate a relativistic waveform for the same physical system using
# the Effective One Body (EOB) approximant SEOBNRv4.
#
# Physical parameters match those used in the Newtonian simulations:
#   M_BH = 4.516 M_sun (black hole mass)
#   M_NS = 1.4   M_sun (neutron star mass)
#   f_lower = 10 Hz (starting frequency for waveform generation)
 
M_BH    = 4.516   # black hole mass in solar masses
M_NS    = 1.4     # neutron star mass in solar masses
 
print("\nGenerating GR waveform with SEOBNRv4...")
hp, hc = get_td_waveform(
    approximant='SEOBNRv4',
    mass1=M_BH,
    mass2=M_NS,
    delta_t=1.0 / 8192,    # sampling rate: 8192 Hz
    f_lower=10             # start frequency in Hz
)
 
# Plot the GR waveform
plt.figure(figsize=(14, 5))
plt.plot(hp.sample_times, hp, label='h+ (Plus)',  color='steelblue')
plt.plot(hc.sample_times, hc, label='h× (Cross)', color='darkorange',
         linestyle='--', alpha=0.8)
plt.xlabel('Time (s)')
plt.ylabel('Strain')
plt.title(f'GR Waveform — SEOBNRv4 | M_BH={M_BH} M☉, M_NS={M_NS} M☉')
plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig('gr_waveform.png', dpi=150)
plt.show()
 
# =============================================================================
# 6. GR PEAK FREQUENCY — FFT NEAR MERGER
# =============================================================================
# Focus FFT on the final 0.5 s of the waveform where the merger occurs,
# to isolate the peak merger frequency from the long inspiral signal.
 
TIME_WINDOW = 0.5   # seconds before merger to analyse
 
times  = hp.sample_times.numpy()
strain = hp.numpy()
 
mask           = times > (times[-1] - TIME_WINDOW)
times_merger   = times[mask]
strain_merger  = strain[mask]
 
gr_freqs, gr_amps, gr_peak_freq = compute_fft(times_merger, strain_merger)
 
plt.figure(figsize=(10, 5))
plt.plot(gr_freqs, gr_amps, color='steelblue', label='Amplitude Spectrum')
plt.axvline(gr_peak_freq, color='red', linestyle='--',
            label=f'Peak: {gr_peak_freq:.2f} Hz')
plt.xlim(50, 300)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (a.u.)')
plt.title('FFT of GR Waveform Near Merger (SEOBNRv4)')
plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig('gr_fft.png', dpi=150)
plt.show()
 
print(f"\nGR peak frequency near merger: {gr_peak_freq:.2f} Hz")
 
# =============================================================================
# 7. COMPARISON SUMMARY
# =============================================================================
 
print("\n" + "=" * 55)
print("NEWTONIAN vs GR — PEAK FREQUENCY COMPARISON")
print("=" * 55)
print(f"{'Configuration':<35} {'Newt. (Hz)':>10}")
print("-" * 55)
for label, freq in zip(LABELS, newtonian_peak_freqs):
    print(f"{label:<35} {freq:>10.2f}")
print("-" * 55)
print(f"{'GR (SEOBNRv4)':<35} {gr_peak_freq:>10.2f}")
print("=" * 55)
print(
    f"\nGR peak frequency is {gr_peak_freq / np.mean(newtonian_peak_freqs):.1f}x "
    f"higher than the Newtonian average — consistent with strong-field "
    f"relativistic corrections dominating near merger."
)
