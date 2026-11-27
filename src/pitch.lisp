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

(deftype %doubles ()
  "The typed vectors of the HMMs' inner loops."
  '(simple-array double-float (*)))

(defun pitch-bin-count (pitch-range)
  "The number of 10-cent pitch bins over PITCH-RANGE (low high). See pitch.lisp."
  (check-type pitch-range (satisfies %pitch-range-p) "a list (low high) of reals, 0 <= low < high < 128")
  (destructuring-bind (low high) pitch-range
    (1+ (floor (* +bins-per-semitone+ (- (rational high) (rational low)))))))

(defun pitch-max-jump (hop sample-rate)
  "The largest pitch change between frames, in bins, at HOP and SAMPLE-RATE. See pitch.lisp."
  (check-type hop (integer 1))
  (check-type sample-rate (integer 1))
  (max 1 (round (* +max-transition-rate+ 12 +bins-per-semitone+ hop) sample-rate)))

;;; DECODE-PITCH's arguments. The bound on the bins needs BIN-COUNT, which is
;;; checked after OBSERVATIONS: OBSERVATIONS is checked for its shape first,
;;; then once more against the bound, both times by CHECK-TYPE.

(defvar *%bin-count* nil
  "While DECODE-PITCH checks its observations against its bin count n: n, the
bound every bin stays below. Otherwise NIL.")

(defun %pitch-observations-p (x)
  "True for a list or vector of frames, each a list of (bin . probability)
conses, BIN an integer >= 0 (and below *%BIN-COUNT* when that is bound) and
PROBABILITY a non-negative real."
  (and (or (vectorp x) (%proper-list-p x))
       (every (lambda (frame)
                (and (%proper-list-p frame)
                     (every (lambda (pair)
                              (and (consp pair)
                                   (typep (car pair) '(integer 0))
                                   (or (null *%bin-count*) (< (car pair) *%bin-count*))
                                   (typep (cdr pair) '(real 0))))
                            frame)))
             x)))

(defun %log-probability (p)
  "log(max(P, +PROBABILITY-FLOOR+)) for a non-negative real P, a DOUBLE-FLOAT.
P is compared with the floor exactly, so a P too small for a double meets it."
  (if (< p +probability-floor+)
      (log +probability-floor+)
      (log (float p 1d0))))

;;; The band step: frame i - 1's deltas PV (voiced) and PU (unvoiced) to frame
;;; i's CV and CU. For each target bin b, one pass over the band b' = b - B' ..
;;; b + B' (B' = min(B, n - 1), clipped to the range) gives the band maxima
;;; MV(b) = max delta_V(b') + log w(b - b') and MU(b) likewise, scanned upward
;;; with a strict > (the smallest b' of equal ones). Then
;;;   delta_V(b) = log O_V(b) + max(MV(b) + log(1 - s), MU(b) + log s)
;;;   delta_U(b) = log O_U    + max(MV(b) + log s, MU(b) + log(1 - s)),
;;; the voiced predecessor winning a tie (its index is the smaller). The
;;; backpointer of target x (V b = b, U b = n + b) goes to BP[BASE + x]. Two
;;; copies: backpointers of 16 bits, and of 32 bits when 2n > 65536.

(macrolet ((define-band-step (name element-type)
             `(defun ,name (pv pu cv cu lw logov logou lstay lswitch n band bp base)
                "One frame of the pitch Viterbi: see the band step above."
                (declare (type %doubles pv pu cv cu lw logov)
                         (type double-float logou lstay lswitch)
                         (type (integer 1 #.(floor most-positive-fixnum 4)) n)
                         (type (integer 0 #.(floor most-positive-fixnum 4)) band)
                         (type (integer 0 #.(floor most-positive-fixnum 2)) base)
                         (type (simple-array ,element-type (*)) bp)
                         (optimize speed (safety 0)))
                (dotimes (b n)
                  (let* ((lo (max 0 (- b band)))
                         (hi (min (1- n) (+ b band)))
                         (l (aref lw (- b lo)))
                         (mv (+ (aref pv lo) l))
                         (av lo)
                         (mu (+ (aref pu lo) l))
                         (au lo))
                    (declare (type fixnum lo hi av au) (type double-float l mv mu))
                    (loop for b2 of-type fixnum from (1+ lo) to b
                          do (let* ((l (aref lw (- b b2)))
                                    (x (+ (aref pv b2) l))
                                    (y (+ (aref pu b2) l)))
                               (when (> x mv) (setf mv x av b2))
                               (when (> y mu) (setf mu y au b2))))
                    (loop for b2 of-type fixnum from (1+ b) to hi
                          do (let* ((l (aref lw (- b2 b)))
                                    (x (+ (aref pv b2) l))
                                    (y (+ (aref pu b2) l)))
                               (when (> x mv) (setf mv x av b2))
                               (when (> y mu) (setf mu y au b2))))
                    (let ((vv (+ mv lstay))
                          (uv (+ mu lswitch)))
                      (if (>= vv uv)
                          (setf (aref cv b) (+ (aref logov b) vv)
                                (aref bp (+ base b)) av)
                          (setf (aref cv b) (+ (aref logov b) uv)
                                (aref bp (+ base b)) (+ n au))))
                    (let ((vu (+ mv lswitch))
                          (uu (+ mu lstay)))
                      (if (>= vu uu)
                          (setf (aref cu b) (+ logou vu)
                                (aref bp (+ base n b)) av)
                          (setf (aref cu b) (+ logou uu)
                                (aref bp (+ base n b)) (+ n au)))))))))
  (define-band-step %pitch-step/16 (unsigned-byte 16))
  (define-band-step %pitch-step/32 (unsigned-byte 32)))

(defun %decode-pitch (frames n jump switch)
  "DECODE-PITCH on checked arguments, FRAMES a SIMPLE-VECTOR of frames."
  (let ((count (length frames)))
    (when (zerop count)
      (return-from %decode-pitch (vector)))
    (let* ((states (* 2 n))
           (band (min jump (1- n)))
           (wide (> states 65536))
           (bp (make-array (* (1- count) states)
                           :element-type (if wide '(unsigned-byte 32) '(unsigned-byte 16))))
           (lw (make-array (1+ band) :element-type 'double-float))
           (lfloor (%log-probability 0))
           (lstay (%log-probability (- 1 switch)))
           (lswitch (%log-probability switch))
           (linit (%log-probability (/ 1 states)))
           (pv (make-array n :element-type 'double-float))
           (pu (make-array n :element-type 'double-float))
           (cv (make-array n :element-type 'double-float))
           (cu (make-array n :element-type 'double-float))
           ;; log O_V of the current frame: the floor but where it has
           ;; candidates. SUMS and SEEN are reset after each frame.
           (logov (make-array n :element-type 'double-float :initial-element lfloor))
           (sums (make-array n :initial-element 0))
           (seen (make-array n :element-type 'bit :initial-element 0))
           (touched '())
           (logou 0d0))
      (declare (type %doubles pv pu cv cu logov lw) (type double-float logou))
      ;; w(d) = (B + 1 - d) / (B + 1)^2 for the true B, d = 0 .. B'.
      (dotimes (d (1+ band))
        (setf (aref lw d) (%log-probability (/ (- (1+ jump) d) (expt (1+ jump) 2)))))
      (flet ((observe (frame)
               ;; O_V(b): the sum of FRAME's probabilities of bin b; O_U =
               ;; (1 - min(1, the sum of them all)) / n.
               (let ((total 0))
                 (dolist (pair frame)
                   (let ((bin (car pair)))
                     (when (zerop (sbit seen bin))
                       (setf (sbit seen bin) 1)
                       (push bin touched))
                     (setf (svref sums bin) (+ (svref sums bin) (cdr pair))
                           total (+ total (cdr pair)))))
                 (dolist (bin touched)
                   (setf (aref logov bin) (%log-probability (svref sums bin))))
                 (setf logou (%log-probability (/ (- 1 (min 1 total)) n)))))
             (forget ()
               (dolist (bin touched)
                 (setf (aref logov bin) lfloor
                       (svref sums bin) 0
                       (sbit seen bin) 0))
               (setf touched '())))
        (observe (svref frames 0))
        (dotimes (b n)
          (setf (aref cv b) (+ linit (aref logov b))
                (aref cu b) (+ linit logou)))
        (forget)
        (loop for i from 1 below count
              do (rotatef pv cv)
                 (rotatef pu cu)
                 (observe (svref frames i))
                 (if wide
                     (%pitch-step/32 pv pu cv cu lw logov logou lstay lswitch n band bp
                                     (* (1- i) states))
                     (%pitch-step/16 pv pu cv cu lw logov logou lstay lswitch n band bp
                                     (* (1- i) states)))
                 (forget)))
      ;; The end state: the largest delta, V before U and lower bins first.
      (let ((end 0)
            (best (aref cv 0)))
        (declare (type double-float best))
        (dotimes (b n)
          (when (> (aref cv b) best) (setf best (aref cv b) end b)))
        (dotimes (b n)
          (when (> (aref cu b) best) (setf best (aref cu b) end (+ n b))))
        (let ((path (make-array count))
              (x end))
          (loop for i from (1- count) downto 0
                do (setf (svref path i) (and (< x n) x))
                   (when (plusp i)
                     (setf x (aref bp (+ (* (1- i) states) x)))))
          path)))))

(defun decode-pitch (observations bin-count max-jump &key (switch-probability +switch-probability+))
  "The Viterbi path of pYIN's pitch and voicing HMM: a bin per voiced frame,
NIL per unvoiced one. See pitch.lisp."
  (check-type observations (satisfies %pitch-observations-p)
              "a list or vector of frames, each a list of (bin . probability) conses")
  (check-type bin-count (integer 1))
  (check-type max-jump (integer 0))
  (check-type switch-probability (real (0) (1)))
  (let ((*%bin-count* bin-count))
    (check-type observations (satisfies %pitch-observations-p)
                "a list or vector of frames of (bin . probability) conses, each bin below BIN-COUNT"))
  (%decode-pitch (coerce observations 'simple-vector) bin-count max-jump switch-probability))

;;; PYIN-CONTOUR

(defun %frame-samples (frame samples start)
  "FRAME, a DOUBLE-FLOAT vector, filled with SAMPLES from index START on (as
doubles), 0 outside SAMPLES."
  (declare (type %doubles frame) (type audio:sample-vector samples) (type fixnum start)
           (optimize speed))
  (let ((size (length frame))
        (length (length samples)))
    (fill frame 0d0)
    (loop for j of-type fixnum from (max 0 (- start)) below (min size (- length start))
          do (setf (aref frame j) (float (aref samples (+ start j)) 1d0)))
    frame))

(defun %candidate-bin (f0 low n)
  "The pitch bin of F0 Hz: round(10 ((hz-to-midi f0) - LOW)), clamped to 0 .. N - 1."
  (max 0 (min (1- n) (round (* +bins-per-semitone+ (- (tp:hz-to-midi f0) low))))))

(defun pyin-contour (audio &key (frame-size 1024) (hop 128) (gate-db -40) (pitch-range '(36 84)))
  "pYIN's f0 contour of AUDIO, a mono PCM-BUFFER: candidates under the beta
threshold prior, a level gate, the pitch HMM's Viterbi path. See pitch.lisp."
  (check-type audio (satisfies %mono-buffer-p) "a one-channel PCM-BUFFER")
  (check-type frame-size (and (integer 4) (satisfies evenp)) "an even integer >= 4")
  (check-type hop (integer 1))
  (check-type gate-db (real * 0))
  (check-type pitch-range (satisfies %pitch-range-p) "a list (low high) of reals, 0 <= low < high < 128")
  (let* ((s (audio:pcm-buffer-sample-rate audio))
         (samples (audio:pcm-buffer-samples audio))
         (w (floor frame-size 2))
         (low (rational (first pitch-range)))
         (high (rational (second pitch-range)))
         (n (pitch-bin-count pitch-range))
         (count (dsp:frame-count (length samples) w hop))
         ;; l_i: the RMS of frame i's integration window, the W samples from
         ;; (frame-start i w hop), zeros outside the audio.
         (levels (dsp:frame-rms audio w hop))
         (loudest (reduce #'max levels :initial-value 0d0))
         (gate (* (audio:db-to-gain gate-db) loudest))
         (tau-lo (max 2 (floor s (tp:midi-to-hz high))))
         (tau-hi (min (1- w) (ceiling s (tp:midi-to-hz low))))
         (frame (make-array frame-size :element-type 'double-float :initial-element 0d0))
         ;; Per frame: its (bin tau* . p) candidates, and its (bin . p) pairs.
         (binned (make-array count :initial-element '()))
         (observations (make-array count :initial-element '()))
         (voicing (make-array count :element-type 'double-float :initial-element 0d0)))
    (declare (type %doubles levels))
    (when (plusp loudest)
      (dotimes (i count)
        (unless (< (aref levels i) gate)
          (%frame-samples frame samples (dsp:frame-start i w hop))
          (multiple-value-bind (candidates voiced)
              (pitch-candidates (dsp:yin-cmndf frame :window w) tau-lo tau-hi)
            (let ((these (mapcar (lambda (c) (list* (%candidate-bin (/ s (car c)) low n) c))
                                 candidates)))
              (setf (svref binned i) these
                    (svref observations i) (mapcar (lambda (e) (cons (first e) (cddr e))) these)
                    (aref voicing i) (float voiced 1d0)))))))
    (let ((path (decode-pitch observations n (pitch-max-jump hop s)))
          (f0 (make-array count :initial-element nil)))
      (dotimes (i count)
        (let ((b (svref path i)))
          (when b
            ;; The most probable candidate of the frame in bin B, the first in
            ;; lag order of equal ones; else the bin's own pitch.
            (let ((best nil))
              (dolist (e (svref binned i))
                (when (and (= (first e) b)
                           (or (null best) (> (cddr e) (cddr best))))
                  (setf best e)))
              (setf (svref f0 i)
                    (if best
                        (/ s (second best))
                        (tp:midi-to-hz (+ low (/ b +bins-per-semitone+)))))))))
      (tp:make-f0-contour :sample-rate s :hop hop :f0 f0 :voicing voicing))))
