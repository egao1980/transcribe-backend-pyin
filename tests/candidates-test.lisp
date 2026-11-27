;;;; THRESHOLD-PRIOR and PITCH-CANDIDATES: pYIN's beta prior over thresholds and
;;;; one frame's candidate periods, exact.

(in-package #:transcribe-backend-pyin/tests)

(defun beta-2-18 (x)
  "The CDF of Beta(2, 18) in closed form: 1 - (1 - x)^19 - 19 x (1 - x)^18."
  (- 1 (expt (- 1 x) 19) (* 19 x (expt (- 1 x) 18))))

(defun cmndf (&rest pairs)
  "A d' vector of 24 lags, 1d0 everywhere but at the (lag value) PAIRS."
  (let ((d (make-array 24 :element-type 'double-float :initial-element 1d0)))
    (loop for (lag value) on pairs by #'cddr
          do (setf (aref d lag) (float value 1d0)))
    d))

(deftest threshold-prior-values
  (let ((p (threshold-prior)))
    (ok (and (simple-vector-p p) (= 100 (length p))))
    (ok (every (lambda (x) (and (typep x 'rational) (plusp x))) p) "positive rationals")
    (ok (= 1 (reduce #'+ p)) "they sum to exactly 1")
    (ok (loop for j from 1 to 100
              always (= (svref p (1- j)) (- (beta-2-18 (/ j 100)) (beta-2-18 (/ (1- j) 100)))))
        "P_j = I(j/100) - I((j-1)/100), Beta(2, 18)")
    (ok (= (svref p 0) (beta-2-18 1/100)))
    (ok (< 1/10 (loop for j from 1 to 100 sum (* (/ j 100) (svref p (1- j)))) 11/100)
        "mean threshold about 1/10")
    (ok (= 6 (1+ (position (reduce #'max p) p))) "the largest mass is on (5/100, 6/100], around the mode 1/18")
    (setf (svref p 0) 0)
    (ok (plusp (svref (threshold-prior) 0)) "a fresh vector each call")))

(deftest pitch-candidates-one-trough
  ;; d'(7) = 1/2, d'(8) = 1/16, d'(9) = 1/4: thresholds 7/100 .. 1 see it; 1/100 ..
  ;; 6/100 don't, and give it 1/100 of their mass, as the deepest trough.
  (multiple-value-bind (candidates voiced)
      (pitch-candidates (cmndf 7 1/2 8 1/16 9 1/4) 2 20)
    (let ((expected (+ (- 1 (beta-2-18 6/100)) (/ (beta-2-18 6/100) 100))))
      (ok (= 1 (length candidates)))
      (ok (close-to (car (first candidates)) 8.2d0 1d-12) "tau* = 8 + (a - c) / (2 (a - 2b + c))")
      (ok (typep (car (first candidates)) 'double-float))
      (ok (eql expected (cdr (first candidates))) "an exact rational")
      (ok (eql expected voiced)))))

(deftest pitch-candidates-two-troughs
  ;; The first trough (lag 5, d' 3/8) takes thresholds 38/100 .. 1; the second
  ;; (lag 12, d' 1/16) takes 7/100 .. 37/100 and, as the deepest, 1/100 of the rest.
  (multiple-value-bind (candidates voiced)
      (pitch-candidates (cmndf 4 1 5 3/8 6 3/4 7 4/5 8 17/20 9 7/8 10 9/10 11 1/2 12 1/16 13 1/2)
                        2 20)
    (ok (= 2 (length candidates)))
    (ok (close-to (car (first candidates)) 5.125d0 1d-12))
    (ok (eql (- 1 (beta-2-18 37/100)) (cdr (first candidates))))
    (ok (close-to (car (second candidates)) 12d0 1d-12))
    (ok (eql (+ (- (beta-2-18 37/100) (beta-2-18 6/100)) (/ (beta-2-18 6/100) 100))
             (cdr (second candidates))))
    (ok (eql (+ (- 1 (beta-2-18 6/100)) (/ (beta-2-18 6/100) 100)) voiced))))

(deftest pitch-candidates-rules
  (testing "a later, shallower trough gets nothing"
    (let ((cs (pitch-candidates (cmndf 5 1/16 6 3/4 11 1/2 12 3/8 13 1/2) 2 20)))
      (ok (= 1 (length cs)))
      (ok (close-to (car (first cs)) 5 1/2))))
  (testing "a trough exactly at a threshold is not below it"
    (multiple-value-bind (cs voiced) (pitch-candidates (cmndf 7 3/4 8 1/2 9 3/4) 2 20)
      (ok (= 1 (length cs)))
      (ok (eql (+ (- 1 (beta-2-18 50/100)) (/ (beta-2-18 50/100) 100)) voiced))))
  (testing "of equal troughs the first is the deepest"
    (let ((cs (pitch-candidates (cmndf 5 1/4 6 1/2 11 1/2 12 1/4 13 1/2) 2 20)))
      (ok (= 1 (length cs)))
      (ok (close-to (car (first cs)) 5 1/2))))
  (testing "no trough below 1: the deepest gets 1/100"
    (multiple-value-bind (cs voiced)
        (pitch-candidates (map 'vector (lambda (x) (+ 3/2 (* 1/100 (abs (- x 9))))) (loop for i below 24 collect i))
                          2 20)
      (ok (= 1 (length cs)))
      (ok (eql 1/100 voiced))
      (ok (close-to (car (first cs)) 9 1/2))))
  (testing "no troughs"
    (ok (equal '(nil 0) (multiple-value-list (pitch-candidates (cmndf) 2 20))) "flat d', as silence")
    (ok (equal '(nil 0) (multiple-value-list
                         (pitch-candidates (map 'vector (lambda (i) (/ (1+ i) 30)) (loop for i below 24 collect i))
                                           2 20)))
        "rising d'"))
  (testing "the lag range"
    (ok (null (pitch-candidates (cmndf 7 1/2 8 1/16 9 1/4) 9 20)) "below TAU-LO")
    (ok (null (pitch-candidates (cmndf 7 1/2 8 1/16 9 1/4) 2 7)) "above TAU-HI")
    (ok (= 1 (length (pitch-candidates (cmndf 7 1/2 8 1/16 9 1/4) -5 500)))
        "clipped to the lags with two neighbours")
    (ok (null (pitch-candidates (cmndf 22 1/2 23 1/16) 2 500)) "the last lag has no neighbour")
    (ok (null (pitch-candidates (cmndf 7 1/2 8 1/16 9 1/4) 12 8)) "an empty range"))
  (testing "argument checks"
    (ok (signals (pitch-candidates (list 1 1/2 1) 1 1) 'type-error))
    (ok (signals (pitch-candidates (cmndf) 1.5 10) 'type-error))
    (ok (signals (pitch-candidates (cmndf) 1 nil) 'type-error))))

(deftest pitch-candidates-real-frames
  (testing "a 200 Hz sine: one sure candidate at its period, 80 samples"
    (let ((frame (sine-samples 200 1024/16000)))
      (multiple-value-bind (cs voiced)
          (pitch-candidates (dsp:yin-cmndf frame :window 512) 15 245)
        (let ((best (reduce (lambda (a b) (if (>= (cdr a) (cdr b)) a b)) cs)))
          (ok (close-to (car best) 80 0.05d0))
          (ok (> voiced 99/100) (format nil "voiced ~,4F" (float voiced)))))))
  (testing "white noise: about 1/100"
    (let* ((random (make-park-miller 9))
           (frame (map 'vector (lambda (x) (declare (ignore x)) (- (* 2 (funcall random)) 1))
                       (make-list 1024))))
      (multiple-value-bind (cs voiced)
          (pitch-candidates (dsp:yin-cmndf frame :window 512) 15 245)
        (declare (ignore cs))
        (ok (< voiced 1/20) (format nil "voiced ~,4F" (float voiced)))))))
