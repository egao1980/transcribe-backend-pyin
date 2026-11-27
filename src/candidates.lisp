;;;; T1 candidates: pYIN's distribution of thresholds (Mauch and Dixon 2014,
;;;; section 2.1), one frame at a time.
;;;;
;;;; YIN picks one period per frame: the first dip of d' below one threshold.
;;;; pYIN tries 100 thresholds under a prior and keeps every period some
;;;; threshold picks, with the prior mass of the thresholds that picked it.
;;;;
;;;; Data:
;;;;   +THRESHOLD-COUNT+        100: thresholds s_j = j / 100, j = 1 .. 100
;;;;   +BETA-PRIOR+             (2 18): the prior over thresholds is Beta(2, 18),
;;;;                            mean 1/10 (pYIN's middle prior; librosa's default)
;;;;   +NO-TROUGH-PROBABILITY+  1/100: the share of a threshold's mass that goes to
;;;;                            the deepest trough when no trough is below it
;;;;                            (pYIN's p_a)
;;;;
;;;; (THRESHOLD-PRIOR) -> a fresh SIMPLE-VECTOR of +THRESHOLD-COUNT+ exact
;;;;   RATIONALs: element j - 1 is P_j = I(j / 100) - I((j - 1) / 100), I the
;;;;   CDF of Beta(a, b), (a b) = +BETA-PRIOR+. For positive integers a and b,
;;;;     I(x) = sum_{k=a}^{n} C(n, k) x^k (1 - x)^(n - k),   n = a + b - 1,
;;;;   exact at a rational x; for (2 18), I(x) = 1 - (1 - x)^19 - 19 x (1 - x)^18.
;;;;   So the P_j are positive and sum to exactly 1; P_1 = I(1/100).
;;;;
;;;; (PITCH-CANDIDATES cmndf tau-lo tau-hi) -> (values candidates voiced).
;;;;   CMNDF is YIN's d'(0) .. d'(T) for one frame (audio-protocol/dsp's
;;;;   YIN-CMNDF), a VECTOR of REALs; TAU-LO and TAU-HI are INTEGERs, the lag
;;;;   range (CHECK-TYPE, in that order). Only lags with two neighbours are
;;;;   searched: max(1, TAU-LO) .. min(TAU-HI, T - 1), none when that range is
;;;;   empty.
;;;;   1. Troughs: the searched lags tau with d'(tau) < d'(tau - 1) and
;;;;      d'(tau) <= d'(tau + 1), in increasing order.
;;;;   2. Thresholds. Each threshold s_j gives its prior mass P_j
;;;;      (THRESHOLD-PRIOR) to the first trough (smallest lag) with d'(tau) <
;;;;      s_j, comparing exactly (d' a float, s_j a rational, as CL compares
;;;;      them). When no trough is below s_j, the deepest trough (the smallest
;;;;      d', the first of equal ones) gets +NO-TROUGH-PROBABILITY+ * P_j. With
;;;;      no trough at all, nothing gets anything.
;;;;   3. Periods. A trough tau that got a positive mass p is a candidate: its
;;;;      period is YIN's parabolic interpolation (step 5) of d' at tau - 1,
;;;;      tau, tau + 1 (a, b, c), tau* = tau + (a - c) / (2 (a - 2b + c)), a
;;;;      DOUBLE-FLOAT within 1/2 of tau (a > b <= c, so the denominator is
;;;;      positive).
;;;;   CANDIDATES is a fresh list of (tau* . p) conses, in increasing order of
;;;;   lag, p an exact RATIONAL > 0; VOICED is the sum of the p, an exact
;;;;   RATIONAL in 0 .. 1: pYIN's voicing probability of the frame. With the
;;;;   deepest trough's d' at v (and m = min(100, max(0, floor(100 v)))), VOICED =
;;;;   1 - I(m / 100) + I(m / 100) / 100: 1 when v < 1/100, 1/100 when no
;;;;   trough is below 1, 0 without troughs.
;;;;
;;;; A frame of a steady tone has one deep trough, at its period, and VOICED
;;;; near 1; white noise has troughs near 1 and VOICED near 1/100; silence
;;;; (d' = 1 everywhere, AU3's rule) has no troughs and VOICED 0.

(in-package #:transcribe-backend-pyin)

(defconstant +threshold-count+ 100
  "pYIN's number of thresholds: s_j = j / 100, j = 1 .. 100.")

(defparameter +beta-prior+ '(2 18)
  "The parameters (a b) of the beta prior over thresholds: mean a / (a + b) = 1/10. Data.")

(defconstant +no-trough-probability+ 1/100
  "pYIN's p_a: the share of a threshold's mass given to the deepest trough when no trough is below it.")

(defun threshold-prior ()
  "The prior mass P_1 .. P_100 of pYIN's thresholds j / 100 under Beta(2, 18),
exact rationals summing to 1. See candidates.lisp."
  (%t1-pending 'threshold-prior))

(defun pitch-candidates (cmndf tau-lo tau-hi)
  "The candidate periods of one frame from its d' CMNDF over lags TAU-LO ..
TAU-HI, with their probabilities, and the frame's voicing probability. See
candidates.lisp."
  (declare (ignorable cmndf tau-lo tau-hi))
  (%t1-pending 'pitch-candidates))
