;;;; T1 notes: from the f0 contour to note events, by a note HMM in the manner
;;;; of Tony (Mauch et al. 2015), simplified.
;;;;
;;;; Tony decodes notes from pYIN's pitch track with an HMM of attack, stable
;;;; and silent states per pitch. Here: one silent state shared by every
;;;; pitch; the attack is a chain of states that sets the minimum note length
;;;; (pruning short notes afterwards would leave a gap that splits the note
;;;; around it); a legato change goes straight from note to note; and the
;;;; level dips of ONSET-FRAMES let a note start again, at its own pitch or
;;;; another (Tony splits on a rise of the frame level). The pitch track is
;;;; smoothed for vibrato first (SMOOTH-PITCHES).
;;;;
;;;; Data:
;;;;   +NOTE-STEPS-PER-SEMITONE+   3: note states every 1/3 semitone (Tony's grid),
;;;;                               so a take sung off A440 still has a state
;;;;                               within 1/6 semitone of each note
;;;;   +SMOOTHING-BAND+            1/2: semitones, how close two frames' medians
;;;;                               are for one to smooth the other (below)
;;;;   +NOTE-SIGMA+                1/2: semitones, the spread of a note's smoothed
;;;;                               pitch around its state
;;;;   +PRIOR-PITCHED+             7/10: Tony's prior probability that a voiced
;;;;                               frame belongs to a note, given weight 1/2
;;;;   +NOTE-ON-PROBABILITY+       1/100: silent -> a note, per frame
;;;;   +NOTE-OFF-PROBABILITY+      1/100: a note -> silent, per frame
;;;;   +NOTE-CHANGE-PROBABILITY+   1/1000: a note -> a far pitch directly (legato),
;;;;                               per frame, shared among the K pitches
;;;;   +MIN-CHANGE-STEPS+          2: a far pitch is at least 2 grid steps away
;;;;                               (2/3 semitone; Tony's minimum distance is 1/2)
;;;;   +ONSET-CHANGE-PROBABILITY+  3/10: at an onset frame, a note -> a far pitch,
;;;;                               shared among the K pitches
;;;;   +REATTACK-PROBABILITY+      3/5: at an onset frame, a note -> a near pitch
;;;;                               again (the same or a grid step away), shared
;;;;                               among the 3 near pitches
;;;;   +PITCH-BAND+                3/2: semitones; a note's pitch is the mean of
;;;;                               its frames within this of their median
;;;;
;;;; (SMOOTH-PITCHES contour &key onsets (vibrato-window 1/5)) -> a fresh
;;;;   SIMPLE-VECTOR with one entry per frame of CONTOUR: the smoothed pitch y_i
;;;;   (a DOUBLE-FLOAT MIDI number) of a voiced frame, NIL for an unvoiced
;;;;   one. CONTOUR an F0-CONTOUR (rate s, HOP), ONSETS a list of integers
;;;;   (those outside 1 .. F - 1 are ignored) and VIBRATO-WINDOW a non-negative
;;;;   REAL, seconds (CHECK-TYPE, in this order). Frame i is voiced when its f0
;;;;   is a number, with pitch m_i = (hz-to-midi f0_i). A segment is a maximal
;;;;   stretch of voiced frames in which no frame but the first is an onset
;;;;   frame: unvoiced frames and onset frames both cut. R = floor(VIBRATO-WINDOW
;;;;   * s / (2 HOP)) frames (12 at the defaults: windows of 25 frames, 200 ms,
;;;;   one period of a 5 Hz vibrato). For a voiced frame i, J_i is the frames j
;;;;   of its own segment with |j - i| <= R, and
;;;;        med_i = the median of m_j over J_i (the mean of the two middle
;;;;                values for an even count), then
;;;;        y_i   = the mean of m_j over the j in J_i with |med_j - med_i| <=
;;;;                +SMOOTHING-BAND+ and |m_j - med_i| <= +PITCH-BAND+ (med_i
;;;;                itself when there is no such j),
;;;;      in DOUBLE-FLOAT. The median keeps a step between two notes sharp; the
;;;;      mean over the frames on the same side of it averages the vibrato out
;;;;      (a 200 ms mean leaves 9 % of a 5.5 Hz vibrato), and the second band
;;;;      keeps octave errors out of it. Cutting at onset frames keeps the
;;;;      vibrato of one note out of the next one's pitch.
;;;;
;;;; (CONTOUR-NOTES contour &key onsets (pitch-range '(36 84)) (min-duration
;;;;                1/10) (vibrato-window 1/5)) -> a fresh SIMPLE-VECTOR of
;;;;   NOTE-EVENTs, one line in time order. CONTOUR an F0-CONTOUR (rate s, HOP,
;;;;   F frames); ONSETS a list of integers, frame indices (ONSET-FRAMES; those
;;;;   outside 1 .. F - 1 are ignored); PITCH-RANGE (low high) as for
;;;;   PITCH-BIN-COUNT; MIN-DURATION and VIBRATO-WINDOW non-negative REALs,
;;;;   seconds (CHECK-TYPE, in this order). No frames, no notes.
;;;;   1. Frames. Frame i is voiced when its f0 is a number: pitch m_i =
;;;;      (hz-to-midi f0_i), voicing w_i (the contour's).
;;;;   2. Vibrato smoothing: y = (smooth-pitches contour :onsets ONSETS
;;;;      :vibrato-window VIBRATO-WINDOW), above.
;;;;   3. The note HMM. Pitches mu_k = low + k / 3, k = 0 .. K - 1, K =
;;;;      floor(3 (high - low)) + 1 (145 for (36 84)). Notes last at least
;;;;      M_min = max(2, ceiling(MIN-DURATION * s / HOP)) frames (13 at the
;;;;      defaults: 104 ms), so each pitch has a chain of M = M_min - 1 attack
;;;;      states A_k^1 .. A_k^M and a stable state S_k; one silent state Q.
;;;;      State order: Q, then for each k in turn A_k^1 .. A_k^M, S_k.
;;;;      Observations at frame i:
;;;;        unvoiced: Q 1; every note state 0 (a note never covers an unvoiced
;;;;          frame);
;;;;        voiced: with pi_i = (w_i + 7/10) / 2, Q has 1 - pi_i, and A_k^j and S_k
;;;;          have pi_i * exp(-(y_i - mu_k)^2 / (2 sigma^2)), sigma = +NOTE-SIGMA+.
;;;;      Transitions into frame i (on, off, change, onset-change, re: the
;;;;      probabilities above; far: |k' - k| >= +MIN-CHANGE-STEPS+, near: the
;;;;      others):
;;;;        Q -> Q: 1 - on;  Q -> A_k^1: on / K
;;;;        A_k^j -> A_k^(j+1) (j < M) and A_k^M -> S_k: 1
;;;;        S_k -> Q: off
;;;;        frame i not an onset frame:
;;;;          S_k -> A_k'^1, k' far: change / K each
;;;;          S_k -> S_k: 1 - off - (far count) change / K
;;;;        frame i an onset frame:
;;;;          S_k -> A_k'^1, k' far: onset-change / K each
;;;;          S_k -> A_k'^1, k' near: re / 3 each (3 = 2 +MIN-CHANGE-STEPS+ - 1)
;;;;          S_k -> S_k: 1 - off - (far count) onset-change / K
;;;;                      - (near count) re / 3
;;;;      Initial (frame 0): Q 1 - on, each A_k^1 on / K, the rest 0.
;;;;      Viterbi in the log domain (log 0 is impossible), in DOUBLE-FLOAT. Ties
;;;;      go to the predecessor of smallest state index, and the path ends in the
;;;;      state of largest delta at the last frame (the smallest index of equal
;;;;      ones). The path never takes an impossible step: staying in Q is
;;;;      always possible.
;;;;   4. Notes. Each frame on the path whose state is some A_k^1 starts a note;
;;;;      the note runs over the frames a .. b up to the next start or the next
;;;;      Q frame, all of them voiced. It becomes a NOTE-EVENT with onset a * HOP
;;;;      / s, offset (b + 1) * HOP / s, pitch the mean of the raw pitches m_a
;;;;      .. m_b that lie within +PITCH-BAND+ of their median (the mean, not the
;;;;      median, because a vibrato's median over a part of its period is off
;;;;      centre; the band leaves out a leap's first frames and octave
;;;;      errors), and confidence the mean of w_a .. w_b (at most 1).
;;;;   So: vibrato of up to a semitone's depth stays one note (the smoothing
;;;;   averages it, and a new note needs a change of at least 2/3 semitone held
;;;;   long enough to pay for the transition); at an onset frame a note starts
;;;;   again, at the same pitch or another, when both parts last M_min frames;
;;;;   nothing shorter than M_min frames is a note; and voiced frames of low
;;;;   voicing (noise: w near 1/100, pi below 1/2) stay silent.

(in-package #:transcribe-backend-pyin)

(defconstant +note-steps-per-semitone+ 3
  "Note states per semitone.")

(defconstant +smoothing-band+ 1/2
  "Semitones: frames whose medians differ by more don't smooth each other.")

(defconstant +note-sigma+ 1/2
  "The spread, in semitones, of a note's smoothed pitch around its state.")

(defconstant +prior-pitched+ 7/10
  "Tony's prior probability that a voiced frame belongs to a note (weight 1/2).")

(defconstant +note-on-probability+ 1/100
  "The probability per frame that a note starts after silence.")

(defconstant +note-off-probability+ 1/100
  "The probability per frame that a note ends in silence.")

(defconstant +note-change-probability+ 1/1000
  "The probability per frame that a note moves straight to a far pitch (shared among the pitches).")

(defconstant +min-change-steps+ 2
  "Far pitches are at least this many grid steps (of 1/3 semitone) away.")

(defconstant +onset-change-probability+ 3/10
  "At an onset frame, the probability that a note moves to a far pitch (shared among the pitches).")

(defconstant +reattack-probability+ 3/5
  "At an onset frame, the probability that a note starts again at a near pitch (shared among 3).")

(defconstant +pitch-band+ 3/2
  "Semitones: a note's pitch is the mean of its frames within this of their median.")

(defun smooth-pitches (contour &key onsets (vibrato-window 1/5))
  "The vibrato-smoothed pitch y_i of each voiced frame of CONTOUR (NIL for an
unvoiced one): a median, then a mean over the frames on the same side of any
step, within segments cut at ONSETS. See notes.lisp."
  (declare (ignorable contour onsets vibrato-window))
  (%t1-pending 'smooth-pitches))

(defun contour-notes (contour &key onsets (pitch-range '(36 84)) (min-duration 1/10)
                                   (vibrato-window 1/5))
  "Note events from CONTOUR: vibrato smoothing, then a note HMM (attack and
stable states per 1/3 semitone, one silent state) that may start a note again
at ONSETS. See notes.lisp."
  (declare (ignorable contour onsets pitch-range min-duration vibrato-window))
  (%t1-pending 'contour-notes))
