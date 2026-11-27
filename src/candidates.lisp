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

;;; The prior. I is kept as a table of its values at the thresholds, C_m =
;;; I(m / 100) for m = 0 .. 100 (C_0 = 0, C_100 = 1), so that the mass of the
;;; thresholds j = m + 1 .. m' is C_m' - C_m. The table is computed once per
;;; value of +BETA-PRIOR+ and never handed out.

(defvar *%threshold-cdf* nil
  "The cache of %THRESHOLD-CDF: (prior . cdf), CDF the SIMPLE-VECTOR of I(m /
100), m = 0 .. 100, under the beta prior PRIOR. Internal.")

(defun %binomial (n k)
  "C(N, K) for integers 0 <= K <= N, exact."
  (let ((c 1))
    ;; After step i, C is C(N - K + i, i), an integer.
    (loop for i from 1 to k
          do (setf c (/ (* c (+ (- n k) i)) i)))
    c))

(defun %beta-cdf (x a b)
  "I(X), the CDF of Beta(A, B) for positive integers A and B, exact at a
rational X: the sum over k = A .. n of C(n, k) X^k (1 - X)^(n - k), n = A + B - 1."
  (let ((n (+ a b -1)))
    (loop for k from a to n
          sum (* (%binomial n k) (expt x k) (expt (- 1 x) (- n k))))))

(defun %threshold-cdf ()
  "The SIMPLE-VECTOR C_0 .. C_100, C_m = I(m / 100) under +BETA-PRIOR+, exact
rationals. Shared: callers never modify it or return it."
  (let ((cache *%threshold-cdf*))
    (if (and cache (equal (car cache) +beta-prior+))
        (cdr cache)
        (destructuring-bind (a b) +beta-prior+
          (let ((cdf (make-array (1+ +threshold-count+))))
            (dotimes (m (1+ +threshold-count+))
              (setf (svref cdf m) (%beta-cdf (/ m +threshold-count+) a b)))
            (setf *%threshold-cdf* (cons (copy-list +beta-prior+) cdf))
            cdf)))))

(defun threshold-prior ()
  "The prior mass P_1 .. P_100 of pYIN's thresholds j / 100 under Beta(2, 18),
exact rationals summing to 1. See candidates.lisp."
  (let ((cdf (%threshold-cdf))
        (p (make-array +threshold-count+)))
    (dotimes (j +threshold-count+ p)
      (setf (svref p j) (- (svref cdf (1+ j)) (svref cdf j))))))

;;; Candidates. Threshold s_j goes to the first trough t with d'_t < j / 100,
;;; so trough t receives exactly the j with d'_t < j / 100 <= u_t, u_t the
;;; smallest d' of the troughs before it (1 for the first): only the troughs
;;; deeper than every earlier one receive anything. With m(x) = min(100,
;;; max(0, floor(100 x))) on the exact rational of x, those j are m(d'_t) + 1
;;; .. m(u_t), a mass of C_m(u_t) - C_m(d'_t). The thresholds below no trough
;;; are j = 1 .. m(v), v the deepest trough's d': p_a C_m(v) more for it.

(defun %real-vector-p (x)
  "True for a VECTOR of REALs."
  (and (vectorp x)
       (or (typep x '(simple-array double-float (*)))
           (every #'realp x))))

(defun %threshold-index (x)
  "m(X) = min(100, max(0, floor(100 X))), exact: X is taken as the rational it
is, so a float near j / 100 lands on the side it lies on."
  (min +threshold-count+ (max 0 (floor (* +threshold-count+ (rational x))))))

(defun %deepening-troughs (d lo hi)
  "The troughs of D among lags LO .. HI (1 <= LO <= HI <= (length D) - 2, so
each lag has both neighbours) that are deeper than every trough before them,
as a fresh list of (tau . d'(tau)), the last (the deepest, the first of equal
ones) first."
  (declare (type fixnum lo hi))
  (macrolet ((scan (&optional type)
               `(let ((d d)
                      (best nil)
                      (troughs '()))
                  ,@(when type `((declare (type ,type d))))
                  (loop for tau of-type fixnum from lo to hi
                        do (let ((v (aref d tau)))
                             (when (and (< v (aref d (1- tau)))
                                        (<= v (aref d (1+ tau)))
                                        (or (null best) (< v best)))
                               (push (cons tau v) troughs)
                               (setf best v))))
                  troughs)))
    (if (typep d '(simple-array double-float (*)))
        (scan (simple-array double-float (*)))
        (scan))))

(defun %parabolic-period (d tau)
  "YIN's parabolic interpolation of D at TAU - 1, TAU, TAU + 1 (a, b, c): tau +
(a - c) / (2 (a - 2b + c)), computed exactly and rounded once to a
DOUBLE-FLOAT. At a trough a > b <= c, so the denominator is positive."
  (let ((a (rational (aref d (1- tau))))
        (b (rational (aref d tau)))
        (c (rational (aref d (1+ tau)))))
    (float (+ tau (/ (- a c) (* 2 (+ (- a (* 2 b)) c)))) 1d0)))

(defun pitch-candidates (cmndf tau-lo tau-hi)
  "The candidate periods of one frame from its d' CMNDF over lags TAU-LO ..
TAU-HI, with their probabilities, and the frame's voicing probability. See
candidates.lisp."
  (check-type cmndf (satisfies %real-vector-p) "a vector of reals")
  (check-type tau-lo integer)
  (check-type tau-hi integer)
  (let* ((lo (max 1 tau-lo))
         (hi (min tau-hi (- (length cmndf) 2)))
         (troughs (and (<= lo hi) (%deepening-troughs cmndf lo hi))))
    (if (null troughs)
        (values '() 0)
        (let ((cdf (%threshold-cdf))
              (upper +threshold-count+)
              (masses '()))
          ;; In lag order: each trough's threshold mass, MASSES ending up
          ;; deepest first, parallel to TROUGHS. UPPER is m(u_t), and after
          ;; the last trough m(v).
          (dolist (trough (reverse troughs))
            (let ((m (%threshold-index (cdr trough))))
              (push (if (< m upper) (- (svref cdf upper) (svref cdf m)) 0) masses)
              (setf upper (min upper m))))
          (incf (first masses) (* +no-trough-probability+ (svref cdf upper)))
          (let ((candidates '())
                (voiced 0))
            ;; Deepest (last in lag order) first, so CANDIDATES comes out in
            ;; lag order.
            (loop for (tau . nil) in troughs
                  for p in masses
                  when (plusp p)
                    do (incf voiced p)
                       (push (cons (%parabolic-period cmndf tau) p) candidates))
            (values candidates voiced))))))
