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

;;; Shared by SMOOTH-PITCHES and CONTOUR-NOTES.

(defun %integer-list-p (x)
  "True for a proper list of integers."
  (and (%proper-list-p x) (every #'integerp x)))

(defun %onset-bits (onsets count)
  "A fresh SIMPLE-BIT-VECTOR of COUNT bits: 1 at each frame of ONSETS within 1
.. COUNT - 1 (the onset frames), 0 elsewhere."
  (let ((bits (make-array count :element-type 'bit :initial-element 0)))
    (dolist (o onsets bits)
      (when (< 0 o count)
        (setf (sbit bits o) 1)))))

(defun %frame-pitches (contour)
  "(values pitches voiced) of CONTOUR: PITCHES the m_i = (hz-to-midi f0_i), a
DOUBLE-FLOAT vector with 0d0 at unvoiced frames, and VOICED a SIMPLE-BIT-VECTOR,
1 where f0 is a number."
  (let* ((f0 (tp:f0-contour-f0 contour))
         (count (length f0))
         (pitches (make-array count :element-type 'double-float :initial-element 0d0))
         (voiced (make-array count :element-type 'bit :initial-element 0)))
    (dotimes (i count (values pitches voiced))
      (let ((f (svref f0 i)))
        (when f
          (setf (aref pitches i) (tp:hz-to-midi f)
                (sbit voiced i) 1))))))

(defun %median (buffer count)
  "The median of BUFFER[0, COUNT), COUNT >= 1, sorting that part in place: the
middle value, or the mean of the two middle values for an even COUNT."
  (declare (type %doubles buffer) (type (integer 1 #.array-dimension-limit) count))
  (if (<= count 32)
      ;; Insertion sort: the windows are short (25 frames at the defaults).
      (loop for i of-type fixnum from 1 below count
            do (let ((v (aref buffer i))
                     (j (1- i)))
                 (declare (type double-float v) (type fixnum j))
                 (loop while (and (>= j 0) (> (aref buffer j) v))
                       do (setf (aref buffer (1+ j)) (aref buffer j))
                          (decf j))
                 (setf (aref buffer (1+ j)) v)))
      (replace buffer (sort (subseq buffer 0 count) #'<)))
  (multiple-value-bind (half odd) (floor count 2)
    (if (= odd 1)
        (aref buffer half)
        (/ (+ (aref buffer (1- half)) (aref buffer half)) 2d0))))

;;; SMOOTH-PITCHES

(defun %smooth (pitches voiced onset-bits radius)
  "The smoothed pitches y of the frame PITCHES (VOICED and ONSET-BITS bit
vectors, windows of RADIUS frames on each side): a fresh SIMPLE-VECTOR, y_i or
NIL. See SMOOTH-PITCHES."
  (declare (type %doubles pitches) (type simple-bit-vector voiced onset-bits)
           (type fixnum radius))
  (let* ((count (length pitches))
         ;; Each voiced frame's segment, first .. last frame.
         (starts (make-array count :element-type 'fixnum :initial-element -1))
         (ends (make-array count :element-type 'fixnum :initial-element -1))
         (medians (make-array count :element-type 'double-float :initial-element 0d0))
         (buffer (make-array (max 1 (min count (1+ (* 2 radius))))
                             :element-type 'double-float :initial-element 0d0))
         ;; Both bands are dyadic, so these comparisons are exact.
         (band (float +smoothing-band+ 1d0))
         (pitch-band (float +pitch-band+ 1d0))
         (y (make-array count :initial-element nil)))
    ;; Segments: an unvoiced frame ends one, an onset frame starts one.
    (let ((start -1))
      (declare (type fixnum start))
      (dotimes (i count)
        (cond ((zerop (sbit voiced i))
               (setf start -1))
              (t
               (when (or (minusp start) (= 1 (sbit onset-bits i)))
                 (setf start i))
               (setf (aref starts i) start)))))
    (let ((end -1))
      (declare (type fixnum end))
      (loop for i from (1- count) downto 0
            do (cond ((zerop (sbit voiced i))
                      (setf end -1))
                     (t
                      (when (minusp end)
                        (setf end i))
                      (setf (aref ends i) end)
                      (when (= (aref starts i) i)
                        (setf end -1))))))
    (flet ((lo (i) (max (aref starts i) (- i radius)))
           (hi (i) (min (aref ends i) (+ i radius))))
      ;; Pass 1: med_i, the median of m_j over J_i.
      (dotimes (i count)
        (when (= 1 (sbit voiced i))
          (let ((lo (lo i))
                (hi (hi i)))
            (replace buffer pitches :start2 lo :end2 (1+ hi))
            (setf (aref medians i) (%median buffer (1+ (- hi lo)))))))
      ;; Pass 2: y_i, the mean of the m_j of J_i on the same side as med_i.
      (dotimes (i count)
        (when (= 1 (sbit voiced i))
          (let ((c (aref medians i))
                (sum 0d0)
                (k 0))
            (declare (type double-float c sum) (type fixnum k))
            (loop for j from (lo i) to (hi i)
                  do (when (and (<= (abs (- (aref medians j) c)) band)
                                (<= (abs (- (aref pitches j) c)) pitch-band))
                       (incf sum (aref pitches j))
                       (incf k)))
            (setf (svref y i) (if (zerop k) c (/ sum k)))))))
    y))

(defun smooth-pitches (contour &key onsets (vibrato-window 1/5))
  "The vibrato-smoothed pitch y_i of each voiced frame of CONTOUR (NIL for an
unvoiced one): a median, then a mean over the frames on the same side of any
step, within segments cut at ONSETS. See notes.lisp."
  (check-type contour tp:f0-contour)
  (check-type onsets (satisfies %integer-list-p) "a list of integers")
  (check-type vibrato-window (real 0))
  (multiple-value-bind (pitches voiced) (%frame-pitches contour)
    (let ((count (length pitches)))
      (%smooth pitches voiced (%onset-bits onsets count)
               ;; R = floor(VIBRATO-WINDOW s / (2 HOP)); no window reaches
               ;; further than the contour.
               (min count (floor (* vibrato-window (tp:f0-contour-sample-rate contour))
                                 (* 2 (tp:f0-contour-hop contour))))))))

;;; CONTOUR-NOTES: the note HMM.
;;;
;;; Which states are possible depends on one count only. With r_i the number
;;; of voiced frames in a row that end at frame i (0 at an unvoiced one), Q is
;;; possible at every frame, A_k^j iff j <= r_i, and S_k iff r_i >= M + 1: an
;;; unvoiced frame makes every note state impossible, A_k^j has the one
;;; predecessor A_k^(j-1), and S_k is reached only through A_k^M. Every
;;; transition the model allows is positive (S_k -> S_k at an onset frame is
;;; at least 9/100) and every emission is finite (1 - pi_i >= 3/20), so the
;;; deltas of the possible states are finite and the impossible ones are
;;; never read: no -infinity and no sentinel. Frame 0 follows a virtual Q of
;;; delta 0 with r = 0: the header's initial distribution is exactly Q's row
;;; of transitions, and 0d0 + x = x.

(defun %log-positive (p)
  "log P for a positive real P, a DOUBLE-FLOAT."
  (log (float p 1d0)))

(defun %note-path (y w voiced onset-bits k-count attacks low)
  "The Viterbi path of the note HMM over the smoothed pitches Y (a
DOUBLE-FLOAT vector, read where VOICED is 1), the voicing W and the onset
frames ONSET-BITS, with K-COUNT pitches from LOW and ATTACKS (M) attack states
each: (values path starts), PATH a FIXNUM vector with each frame's pitch k (-1
for Q) and STARTS a bit vector, 1 where the frame's state is an A_k^1."
  (declare (type %doubles y w) (type simple-bit-vector voiced onset-bits)
           (type (integer 1 1000) k-count) (type (integer 1 #.array-dimension-limit) attacks))
  (let* ((count (length y))
         (m attacks)
         (h (1- +min-change-steps+))            ; near: |k' - k| <= h
         (near-share (1- (* 2 +min-change-steps+)))
         (on +note-on-probability+)
         (off +note-off-probability+)
         (l-qq (%log-positive (- 1 on)))
         (l-on (%log-positive (/ on k-count)))
         (l-off (%log-positive off))
         (l-far (%log-positive (/ +note-change-probability+ k-count)))
         (l-far-onset (%log-positive (/ +onset-change-probability+ k-count)))
         (l-near (%log-positive (/ +reattack-probability+ near-share)))
         (l-stay (make-array k-count :element-type 'double-float))
         (l-stay-onset (make-array k-count :element-type 'double-float))
         (mu (make-array k-count :element-type 'double-float))
         (inverse (float (/ 1 (* 2 (expt +note-sigma+ 2))) 1d0))   ; 1 / (2 sigma^2)
         (prior (float +prior-pitched+ 1d0))
         ;; Frame i - 1: delta(Q), r, delta(A_k^j) at (j - 1) K + k, delta(S_k).
         (q 0d0)
         (run 0)
         (chain (make-array (* m k-count) :element-type 'double-float :initial-element 0d0))
         (s-prev (make-array k-count :element-type 'double-float :initial-element 0d0))
         (s-cur (make-array k-count :element-type 'double-float :initial-element 0d0))
         (e (make-array k-count :element-type 'double-float :initial-element 0d0))
         (pre-v (make-array k-count :element-type 'double-float :initial-element 0d0))
         (pre-a (make-array k-count :element-type 'fixnum :initial-element 0))
         (suf-v (make-array k-count :element-type 'double-float :initial-element 0d0))
         (suf-a (make-array k-count :element-type 'fixnum :initial-element 0))
         ;; Backpointers. Q's and A_k^1's: 0 for Q, 1 + k' for S_k'. S_k's: 0
         ;; for A_k^M, 1 for S_k. A_k^j (j > 1) has one predecessor.
         (bp-q (make-array count :element-type '(unsigned-byte 16) :initial-element 0))
         (bp-a (make-array (* count k-count) :element-type '(unsigned-byte 16) :initial-element 0))
         (bp-s (make-array (* count k-count) :element-type 'bit :initial-element 0)))
    (declare (type double-float l-qq l-on l-off l-far l-far-onset l-near inverse prior q)
             (type fixnum run)
             (type %doubles l-stay l-stay-onset mu chain s-prev s-cur e pre-v suf-v)
             (type (simple-array fixnum (*)) pre-a suf-a)
             (type (simple-array (unsigned-byte 16) (*)) bp-q bp-a)
             (type simple-bit-vector bp-s))
    (dotimes (k k-count)
      (let* ((near (1+ (- (min (1- k-count) (+ k h)) (max 0 (- k h)))))
             (far (- k-count near)))
        (setf (aref l-stay k)
              (%log-positive (- 1 off (* far (/ +note-change-probability+ k-count))))
              (aref l-stay-onset k)
              (%log-positive (- 1 off
                                (* far (/ +onset-change-probability+ k-count))
                                (* near (/ +reattack-probability+ near-share))))
              (aref mu k) (float (+ low (/ k +note-steps-per-semitone+)) 1d0))))
    (dotimes (i count)
      (let* ((rp run)
             (s-live (> rp m))                  ; S possible at i - 1
             (voiced-p (= 1 (sbit voiced i)))
             (onset-p (= 1 (sbit onset-bits i)))
             (pi-i (if voiced-p (* 0.5d0 (+ (aref w i) prior)) 0d0))
             (q-new 0d0))
        (declare (type fixnum rp) (type double-float pi-i q-new))
        ;; Q, from Q, then S_0 .. S_{K-1}.
        (let ((best (+ q l-qq))
              (arg 0))
          (declare (type double-float best) (type fixnum arg))
          (when s-live
            (dotimes (k k-count)
              (let ((v (+ (aref s-prev k) l-off)))
                (when (> v best)
                  (setf best v arg (1+ k))))))
          (setf q-new (+ (if voiced-p (log (- 1d0 pi-i)) 0d0) best)
                (aref bp-q i) arg))
        (cond
          ((not voiced-p)
           (setf run 0))
          (t
           (let ((lpi (log pi-i))
                 (yi (aref y i)))
             (declare (type double-float lpi yi))
             (dotimes (k k-count)
               (let ((d (- yi (aref mu k))))
                 (setf (aref e k) (- lpi (* inverse d d))))))
           ;; S_k, from A_k^M (the smaller index, so it wins a tie), then S_k.
           (when (>= rp m)
             (let ((stay (if onset-p l-stay-onset l-stay))
                   (row (* (1- m) k-count)))
               (declare (type %doubles stay) (type fixnum row))
               (dotimes (k k-count)
                 (let ((best (aref chain (+ row k)))
                       (arg 0))
                   (declare (type double-float best) (type bit arg))
                   (when s-live
                     (let ((v (+ (aref s-prev k) (aref stay k))))
                       (when (> v best)
                         (setf best v arg 1))))
                   (setf (aref s-cur k) (+ (aref e k) best)
                         (sbit bp-s (+ (* i k-count) k)) arg)))))
           ;; The far candidates delta(S_k') + log T: prefix maxima over k' <=
           ;; t (upward, strict >) and suffix maxima over k' >= t (downward,
           ;; >=), each keeping the smallest k' of equal values.
           (when s-live
             (let ((c (if onset-p l-far-onset l-far)))
               (let ((bv 0d0) (ba 0))
                 (declare (type double-float bv) (type fixnum ba))
                 (dotimes (k k-count)
                   (let ((v (+ (aref s-prev k) c)))
                     (when (or (zerop k) (> v bv))
                       (setf bv v ba k)))
                   (setf (aref pre-v k) bv (aref pre-a k) ba)))
               (let ((bv 0d0) (ba 0))
                 (declare (type double-float bv) (type fixnum ba))
                 (loop for k of-type fixnum from (1- k-count) downto 0
                       do (let ((v (+ (aref s-prev k) c)))
                            (when (or (= k (1- k-count)) (>= v bv))
                              (setf bv v ba k)))
                          (setf (aref suf-v k) bv (aref suf-a k) ba)))))
           ;; The chains, top down: A_k^j = e_k + A_k^(j-1) of frame i - 1.
           (loop for j of-type fixnum from (min (1+ rp) m) downto 2
                 do (let ((to (* (1- j) k-count))
                          (from (* (- j 2) k-count)))
                      (declare (type fixnum to from))
                      (dotimes (k k-count)
                        (setf (aref chain (+ to k)) (+ (aref e k) (aref chain (+ from k)))))))
           ;; A_k^1, its predecessors in index order: Q; the far S_k' <= k - 2;
           ;; at an onset frame the near S_{k-1}, S_k, S_{k+1}; the far S_k' >= k + 2.
           (dotimes (k k-count)
             (let ((best (+ q l-on))
                   (arg 0))
               (declare (type double-float best) (type fixnum arg))
               (when s-live
                 (let ((below (- k +min-change-steps+))
                       (above (+ k +min-change-steps+)))
                   (declare (type fixnum below above))
                   (when (>= below 0)
                     (let ((v (aref pre-v below)))
                       (when (> v best)
                         (setf best v arg (1+ (aref pre-a below))))))
                   (when onset-p
                     (loop for k2 of-type fixnum from (max 0 (- k h)) to (min (1- k-count) (+ k h))
                           do (let ((v (+ (aref s-prev k2) l-near)))
                                (when (> v best)
                                  (setf best v arg (1+ k2))))))
                   (when (< above k-count)
                     (let ((v (aref suf-v above)))
                       (when (> v best)
                         (setf best v arg (1+ (aref suf-a above))))))))
               (setf (aref chain k) (+ (aref e k) best)
                     (aref bp-a (+ (* i k-count) k)) arg)))
           (setf run (1+ rp))))
        (rotatef s-prev s-cur)
        (setf q q-new)))
    ;; The end state: the largest delta of the possible states, in state
    ;; order (Q; then per k, A_k^1 .. A_k^M, S_k), the first of equal ones.
    (let ((best q)
          (kind :q)
          (bk 0)
          (bj 0))
      (declare (type double-float best) (type fixnum bk bj))
      (dotimes (k k-count)
        (loop for j of-type fixnum from 1 to (min run m)
              do (let ((v (aref chain (+ (* (1- j) k-count) k))))
                   (when (> v best)
                     (setf best v kind :a bk k bj j))))
        (when (> run m)
          (let ((v (aref s-prev k)))
            (when (> v best)
              (setf best v kind :s bk k)))))
      ;; Backtrack.
      (let ((path (make-array count :element-type 'fixnum :initial-element -1))
            (starts (make-array count :element-type 'bit :initial-element 0)))
        (flet ((enter (p)
                 ;; P from a backpointer of Q or A^1: 0 for Q, 1 + k' for S_k'.
                 (if (zerop p)
                     (setf kind :q)
                     (setf kind :s bk (1- p)))))
          (loop for i from (1- count) downto 0
                do (ecase kind
                     (:q
                      (setf (aref path i) -1)
                      (enter (aref bp-q i)))
                     (:a
                      (setf (aref path i) bk)
                      (cond ((= bj 1)
                             (setf (sbit starts i) 1)
                             (enter (aref bp-a (+ (* i k-count) bk))))
                            (t
                             (decf bj))))
                     (:s
                      (setf (aref path i) bk)
                      (when (zerop (sbit bp-s (+ (* i k-count) bk)))
                        (setf kind :a bj m))))))
        (values path starts)))))

(defun %note-event (a b pitches w hop s)
  "The NOTE-EVENT of frames A .. B: onset A HOP / S, offset (B + 1) HOP / S, the
mean pitch of those within +PITCH-BAND+ of their median, the mean voicing."
  (declare (type %doubles pitches w))
  (let* ((count (1+ (- b a)))
         (median (%median (subseq pitches a (1+ b)) count))
         (band (float +pitch-band+ 1d0))
         (sum 0d0)
         (k 0)
         (voicing 0d0))
    (declare (type double-float median sum voicing) (type fixnum k))
    (loop for i from a to b
          do (let ((x (aref pitches i)))
               (when (<= (abs (- x median)) band)
                 (incf sum x)
                 (incf k))
               (incf voicing (aref w i))))
    (tp:make-note-event :onset (/ (* a hop) s)
                        :offset (/ (* (1+ b) hop) s)
                        :pitch (if (zerop k) median (/ sum k))
                        :confidence (min 1d0 (/ voicing count)))))

(defun %path-notes (path starts pitches w hop s)
  "The notes of a note-HMM path: each A_k^1 frame starts a note, which runs
up to the next start or the next Q frame. A fresh SIMPLE-VECTOR."
  (let ((notes '())
        (open nil))
    (flet ((finish (b)
             (when open
               (push (%note-event open b pitches w hop s) notes)
               (setf open nil))))
      (dotimes (i (length path))
        (cond ((minusp (aref path i))
               (finish (1- i)))
              ((= 1 (sbit starts i))
               (finish (1- i))
               (setf open i))))
      (finish (1- (length path))))
    (coerce (nreverse notes) 'simple-vector)))

(defun contour-notes (contour &key onsets (pitch-range '(36 84)) (min-duration 1/10)
                                   (vibrato-window 1/5))
  "Note events from CONTOUR: vibrato smoothing, then a note HMM (attack and
stable states per 1/3 semitone, one silent state) that may start a note again
at ONSETS. See notes.lisp."
  (check-type contour tp:f0-contour)
  (check-type onsets (satisfies %integer-list-p) "a list of integers")
  (check-type pitch-range (satisfies %pitch-range-p) "a list (low high) of reals, 0 <= low < high < 128")
  (check-type min-duration (real 0))
  (check-type vibrato-window (real 0))
  (let ((count (tp:contour-length contour)))
    (if (zerop count)
        (vector)
        (let* ((s (tp:f0-contour-sample-rate contour))
               (hop (tp:f0-contour-hop contour))
               (w (tp:f0-contour-voicing contour))
               (low (rational (first pitch-range)))
               (high (rational (second pitch-range)))
               (k-count (1+ (floor (* +note-steps-per-semitone+ (- high low)))))
               (attacks (1- (max 2 (ceiling (* min-duration s) hop))))
               (smoothed (smooth-pitches contour :onsets onsets :vibrato-window vibrato-window))
               (y (make-array count :element-type 'double-float :initial-element 0d0)))
          (dotimes (i count)
            (let ((v (svref smoothed i)))
              (when v
                (setf (aref y i) v))))
          (multiple-value-bind (pitches voiced) (%frame-pitches contour)
            (multiple-value-bind (path starts)
                (%note-path y w voiced (%onset-bits onsets count) k-count attacks low)
              (%path-notes path starts pitches w hop s)))))))
