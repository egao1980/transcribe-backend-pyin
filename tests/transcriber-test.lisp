;;;; PYIN-TRANSCRIBER: the class, its initargs and their checks, and the data.
;;;; Given in the contract: these pass on the stubs.

(in-package #:transcribe-backend-pyin/tests)

(deftest transcriber-class
  (let ((p (pyin)))
    (ok (typep p 'transcriber) "a transcribe-protocol TRANSCRIBER")
    (ok (= 1024 (pyin-transcriber-frame-size p)))
    (ok (= 128 (pyin-transcriber-hop p)))
    (ok (= -40 (pyin-transcriber-gate-db p)))
    (ok (eql 1/10 (pyin-transcriber-min-duration p)))
    (ok (eql 1/5 (pyin-transcriber-vibrato-window p)))
    (ok (= 16000 (transcriber-sample-rate p)) "the protocol's default rate")
    (ok (equal '(:name "pyin-transcriber" :frame-size 1024 :hop 128 :gate-db -40
                 :min-duration 1/10 :vibrato-window 1/5)
               (transcriber-info p))))
  (let ((p (pyin :frame-size 2048 :hop 256 :gate-db -60 :min-duration 0 :vibrato-window 1/4)))
    (ok (equal '(:name "pyin-transcriber" :frame-size 2048 :hop 256 :gate-db -60
                 :min-duration 0 :vibrato-window 1/4)
               (transcriber-info p)))
    (ok (not (eq (transcriber-info p) (transcriber-info p))) "a fresh list")))

(deftest-parametrize transcriber-initarg-checks
    ((initarg value)
     :ids ("odd-frame" "small-frame" "float-frame" "hop-0" "float-hop" "gate-positive"
           "gate-not-real" "min-duration-negative" "min-duration-not-real"
           "vibrato-window-negative")
     (:frame-size 1023)
     (:frame-size 2)
     (:frame-size 1024.0)
     (:hop 0)
     (:hop 12.5)
     (:gate-db 1)
     (:gate-db "-40")
     (:min-duration -1/10)
     (:min-duration :short)
     (:vibrato-window -1/5))
  (ok (signals (pyin initarg value) 'type-error)))

(deftest transcriber-data
  (ok (= 100 +threshold-count+))
  (ok (equal '(2 18) +beta-prior+))
  (ok (eql 1/100 +no-trough-probability+))
  (ok (= 10 +bins-per-semitone+))
  (ok (eql 898/25 +max-transition-rate+) "35.92 octaves a second")
  (ok (eql 1/100 +switch-probability+))
  (ok (eql 1d-300 +probability-floor+))
  (ok (= 6 +dip-db+))
  (ok (eql 1/20 +dip-span+))
  (ok (= 3 +note-steps-per-semitone+))
  (ok (eql 1/2 +smoothing-band+))
  (ok (eql 1/2 +note-sigma+))
  (ok (eql 7/10 +prior-pitched+))
  (ok (eql 1/100 +note-on-probability+))
  (ok (eql 1/100 +note-off-probability+))
  (ok (eql 1/1000 +note-change-probability+))
  (ok (= 2 +min-change-steps+))
  (ok (eql 3/10 +onset-change-probability+))
  (ok (eql 3/5 +reattack-probability+))
  (ok (eql 3/2 +pitch-band+)))
