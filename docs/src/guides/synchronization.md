# Synchronization

## Carrier recovery

`freqoffset` applies a known impairment and `carriersync` estimates and removes
it. Choose the phase detector from the actual constellation geometry:

- `:bpsk` for axis-aligned two-point PSK
- `:psk8` for axis-and-diagonal-aligned 8-PSK
- `:qam` for supported QAM or diagonal QPSK such as
  `pskmod(bits, 4; phase_offset=pi/4)`

`sps` changes loop gain; it does not resample. Use `sps=1` after symbol-rate
decimation. Carrier recovery has rotational ambiguity, so framed receivers need
a preamble, pilots, or differential encoding to resolve absolute phase.

Assess lock from both the phase-estimate tail slope and post-acquisition data
quality. A loop that merely improves BER is not necessarily healthy.

## Symbol timing recovery

`timingoffset` produces a controlled fractional delay. `symbolsync` accepts an
oversampled stream and returns a variable-length symbol vector plus one timing
estimate per input sample.

Detector guidance:

- Gardner and early-late are non-data-aided.
- Zero-crossing and Mueller-Muller are decision-directed.
- Start with Gardner for a pulse-shaped QPSK/QAM acquisition experiment.

The timing estimate is fractional and wraps in `[0, 1)`. Evaluate lock with a
wrapped tail spread, not a direct maximum-minus-minimum across a wrap.

## Loop tuning

Smaller `loop_bandwidth` reduces noise but slows acquisition. Larger bandwidth
tracks faster drift but passes more jitter. `damping_factor` controls transient
shape. Tune on the actual SNR, sample rate, pulse shape, and offset range; there
is no universally optimal pair.
