;;;; T1 pitch: pYIN's hidden Markov model over pitch and voicing (Mauch and
;;;; Dixon 2014, section 2.2), and the f0 contour it decodes.
;;;;
;;;; Data:
;;;;   +BINS-PER-SEMITONE+    10: pitch bins of 10 cents (librosa's resolution 0.1)
;;;;   +MAX-TRANSITION-RATE+  898/25 = 35.92 octaves a second, librosa's default:
;;;;                          pYIN's 25 bins every 256 samples at 44.1 kHz
;;;;   +SWITCH-PROBABILITY+   1/100: the probability of a voicing change
;;;;                          between two frames (pYIN's and librosa's)
;;;;   +PROBABILITY-FLOOR+    1d-300: probabilities enter the logarithms as
;;;;                          max(p, floor), so no observation is impossible
;;;;
;;;; (PITCH-BIN-COUNT pitch-range) -> n = floor(10 (high - low)) + 1 for
;;;;   PITCH-RANGE (low high), MIDI reals with 0 <= low < high < 128 (else
;;;;   TYPE-ERROR). Bin b is the pitch low + b / 10. The default range (36 84),
;;;;   C2 to C6, has 481 bins.
;;;; (PITCH-MAX-JUMP hop sample-rate) -> B = max(1, round(35.92 * 12 * 10 * HOP /
;;;;   SAMPLE-RATE)) bins, the largest pitch change between two frames, HOP and
;;;;   SAMPLE-RATE integers >= 1 (else TYPE-ERROR); ROUND is CL's (halves to
;;;;   even). 34 bins, 3.4 semitones, at hop 128 and 16 kHz.
;;;;
;;;; (DECODE-PITCH observations bin-count max-jump &key (switch-probability
;;;;               +switch-probability+)) -> a fresh SIMPLE-VECTOR with one entry
;;;;   per frame: the decoded pitch bin (an integer 0 .. n - 1) where the frame
;;;;   is voiced, NIL where it is unvoiced.
;;;;   Arguments (CHECK-TYPE, in this order): OBSERVATIONS a list or vector
;;;;   with one element per frame, each a list of (bin . probability) conses,
;;;;   BIN an integer 0 .. n - 1 and PROBABILITY a non-negative REAL; BIN-COUNT
;;;;   n an integer >= 1; MAX-JUMP B an integer >= 0; SWITCH-PROBABILITY s a
;;;;   REAL with 0 < s < 1.
;;;;   The model. 2n states: voiced (V, b) is state b and unvoiced (U, b) is
;;;;   state n + b, for b = 0 .. n - 1 (the unvoiced states keep the pitch
;;;;   through a voiceless stretch, as in pYIN).
;;;;     Observations of frame i: O(V, b) = the sum of the probabilities paired
;;;;       with bin b; O(U, b) = (1 - v_i) / n for every b, v_i = min(1, the sum
;;;;       of all the frame's probabilities).
;;;;     Transitions from (x', b') to (x, b): w(b - b') * S(x', x), with the
;;;;       triangular band w(d) = (B + 1 - |d|) / (B + 1)^2 for |d| <= B and no
;;;;       transition beyond it (not renormalized at the ends of the range), and
;;;;       S = 1 - s when the voicing stays and s when it changes.
;;;;     Initial: 1 / (2n) for every state.
;;;;   Viterbi in the log domain, each probability p taken as
;;;;   log(max(p, +PROBABILITY-FLOOR+)) (a band transition is never zero):
;;;;     delta_0(x) = log(1 / 2n) + log O_0(x)
;;;;     delta_i(x) = log O_i(x) + max_{x'} (delta_{i-1}(x') + log T(x', x))
;;;;   in DOUBLE-FLOAT. Ties go to the predecessor of smallest state index
;;;;   (voiced before unvoiced, lower bins first), and the path ends in the
;;;;   state of largest delta at the last frame (the smallest index of equal
;;;;   ones). Backtracking gives one state per frame: its bin when voiced,
;;;;   else NIL. No frames, no entries.
;;;;   The cost is O(n B) per frame: the maximum over x' factors into a
;;;;   triangular max-convolution of the voiced deltas and one of the unvoiced
;;;;   deltas, shared by both voicings of the target.
;;;;
;;;; (PYIN-CONTOUR audio &key (frame-size 1024) (hop 128) (gate-db -40)
;;;;               (pitch-range '(36 84))) -> an F0-CONTOUR of rate s (the
;;;;   audio's) and HOP. AUDIO a one-channel PCM-BUFFER; FRAME-SIZE N an even
;;;;   integer >= 4; HOP an integer >= 1; GATE-DB a REAL <= 0; PITCH-RANGE
;;;;   (low high) as for PITCH-BIN-COUNT (CHECK-TYPE, in this order). With L
;;;;   samples x_0 .. x_{L-1}, W = N / 2, n = (pitch-bin-count pitch-range):
;;;;   1. Frames, as the YIN baseline's: F = ceiling(L / HOP) (0 for L = 0).
;;;;      Frame i is the N samples from x_{c - W/2}, c = i * HOP
;;;;      ((dsp:frame-start i W hop)), as doubles, 0 outside the audio, so
;;;;      the integration window, its first W samples, is centred on frame
;;;;      i's time i * HOP / s.
;;;;   2. Gate. l_i = the RMS of frame i's first W samples. A frame with l_i
;;;;      below (db-to-gain GATE-DB) times the largest l_i has no candidates,
;;;;      and neither has any frame when that largest l_i is 0.
;;;;   3. Candidates. tau_lo = max(2, floor(s / (midi-to-hz high))) and
;;;;      tau_hi = min(W - 1, ceiling(s / (midi-to-hz low))) (15 .. 245 for
;;;;      (36 84) at 16 kHz). Each frame the gate lets in gets
;;;;      (pitch-candidates (dsp:yin-cmndf frame :window W) tau_lo tau_hi).
;;;;      A candidate (tau* . p) falls in bin round(10 ((hz-to-midi (/ s tau*))
;;;;      - low)), clamped to 0 .. n - 1.
;;;;   4. Decoding. (decode-pitch obs n (pitch-max-jump hop s)), obs_i the
;;;;      frame's (bin . p) pairs (none for a gated frame).
;;;;   5. The contour. A voiced frame of bin b gets the f0 s / tau* of its most
;;;;      probable candidate in bin b (the first in lag order of equal ones), or
;;;;      (midi-to-hz (+ low (/ b 10))) when no candidate of the frame is in b;
;;;;      an unvoiced frame gets NIL. Every frame's voicing is its VOICED
;;;;      (PITCH-CANDIDATES's second value, 0 for a gated frame) as a
;;;;      DOUBLE-FLOAT, also where the decoded state is unvoiced: pYIN's
;;;;      voicing probability, the evidence a review screen shows.

(in-package #:transcribe-backend-pyin)

(defconstant +bins-per-semitone+ 10
  "Pitch bins per semitone: 10 cents a bin.")

(defconstant +max-transition-rate+ 898/25
  "The fastest pitch change between frames, in octaves a second (librosa's 35.92).")

(defconstant +switch-probability+ 1/100
  "The probability of a change of voicing between two frames.")

(defconstant +probability-floor+ 1d-300
  "Probabilities enter the HMM's logarithms as max(p, +PROBABILITY-FLOOR+).")

(defun pitch-bin-count (pitch-range)
  "The number of 10-cent pitch bins over PITCH-RANGE (low high). See pitch.lisp."
  (declare (ignorable pitch-range))
  (%t1-pending 'pitch-bin-count))

(defun pitch-max-jump (hop sample-rate)
  "The largest pitch change between frames, in bins, at HOP and SAMPLE-RATE. See pitch.lisp."
  (declare (ignorable hop sample-rate))
  (%t1-pending 'pitch-max-jump))

(defun decode-pitch (observations bin-count max-jump &key (switch-probability +switch-probability+))
  "The Viterbi path of pYIN's pitch and voicing HMM: a bin per voiced frame,
NIL per unvoiced one. See pitch.lisp."
  (declare (ignorable observations bin-count max-jump switch-probability))
  (%t1-pending 'decode-pitch))

(defun pyin-contour (audio &key (frame-size 1024) (hop 128) (gate-db -40) (pitch-range '(36 84)))
  "pYIN's f0 contour of AUDIO, a mono PCM-BUFFER: candidates under the beta
threshold prior, a level gate, the pitch HMM's Viterbi path. See pitch.lisp."
  (declare (ignorable audio frame-size hop gate-db pitch-range))
  (%t1-pending 'pyin-contour))
