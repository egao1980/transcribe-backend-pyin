;;;; CONTOUR-NOTES: vibrato smoothing and the note HMM, on contours made by hand
;;;; (125 frames a second: hop 128 at 16 kHz).

(in-package #:transcribe-backend-pyin/tests)

(defun spans (notes)
  "Each note as (onset-frame offset-frame), frames of 8 ms."
  (map 'list (lambda (n) (list (* 125 (note-event-onset n)) (* 125 (note-event-offset n)))) notes))

(defun vibrato (centre depth rate)
  "A pitch function: CENTRE with a vibrato of DEPTH cents at RATE Hz."
  (lambda (time) (+ centre (* depth 1/100 (sin (* 2 pi rate time))))))

(deftest contour-notes-steady
  ;; Frames 25 .. 149 at MIDI 60 with voicing 0.9, silence around.
  (let* ((c (contour-of (steps '(nil 1/5) '(60 1) '(nil 1/5)) 175 :voicing 9/10))
         (ns (contour-notes c)))
    (ok (= 1 (length ns)))
    (ok (equal '((25 150)) (spans ns)) "onset at the first voiced frame, offset after the last")
    (ok (eql 1/5 (note-event-onset (svref ns 0))))
    (ok (eql 6/5 (note-event-offset (svref ns 0))))
    (ok (close-to (note-event-pitch (svref ns 0)) 60 1d-9))
    (ok (close-to (note-event-confidence (svref ns 0)) 9/10 1d-9) "the mean voicing")
    (ok (simple-vector-p ns))))

(deftest contour-notes-steps
  (testing "a whole tone and a semitone, legato"
    (let ((ns (contour-notes (contour-of (steps '(nil 1/5) '(60 1/2) '(62 1/2) '(63 1/2) '(nil 1/5)) 225))))
      (ok (equal '((25 88) (88 150) (150 213)) (spans ns)))
      (ok (equal '(60 62 63) (mapcar #'round (pitches ns))))))
  (testing "a neighbour figure in eighths: 71 72 71"
    (let ((ns (contour-notes (contour-of (steps '(nil 1/5) '(71 1/4) '(72 1/4) '(71 1/4) '(nil 1/5)) 150))))
      (ok (equal '(71 72 71) (mapcar #'round (pitches ns))))))
  (testing "a leap of a twelfth"
    (let ((ns (contour-notes (contour-of (steps '(nil 1/5) '(48 1/2) '(67 1/2) '(nil 1/5)) 175))))
      (ok (equal '((25 88) (88 150)) (spans ns)))))
  (testing "a drift of 30 cents within a note is no step"
    (let ((ns (contour-notes (contour-of (lambda (time) (and (< 1/5 time 2) (+ 60 (* 1/6 time)))) 275))))
      (ok (= 1 (length ns))))))

(deftest-parametrize contour-notes-vibrato
    ((depth rate)
     :ids ("60c-5.5Hz" "100c-5.5Hz" "100c-5Hz" "100c-7Hz")
     (60 11/2) (100 11/2) (100 5) (100 7))
  ;; Two seconds of one note with vibrato: one note, its pitch the centre.
  (let* ((pitch (vibrato 64 depth rate))
         (ns (contour-notes (contour-of (lambda (time) (and (< 1/5 time 11/5) (funcall pitch time))) 300))))
    (ok (= 1 (length ns)) (format nil "~D notes" (length ns)))
    (ok (and (= 1 (length ns)) (close-to (note-event-pitch (svref ns 0)) 64 1/10)) "within 10 cents")))

(deftest contour-notes-onsets
  (let ((c (contour-of (steps '(nil 1/5) '(60 1) '(nil 1/5)) 175)))
    (testing "a repeated pitch splits at an onset frame"
      (ok (equal '((25 88) (88 150)) (spans (contour-notes c :onsets '(88))))))
    (testing "only where both parts last the minimum (13 frames)"
      (ok (equal '((25 150)) (spans (contour-notes c :onsets '(30)))) "5 frames after the start")
      (ok (equal '((25 150)) (spans (contour-notes c :onsets '(140)))) "10 frames before the end")
      (ok (equal '((25 38) (38 150)) (spans (contour-notes c :onsets '(38)))) "13 frames: enough"))
    (testing "several onsets, and onsets outside the contour"
      (ok (equal '((25 63) (63 100) (100 150)) (spans (contour-notes c :onsets '(-4 0 63 100 175 900))))))
    (testing "an onset at a pitch change keeps the change there"
      (let ((c2 (contour-of (steps '(nil 1/5) '(60 1/2) '(61 1/2) '(nil 1/5)) 175)))
        (ok (equal '((25 88) (88 150)) (spans (contour-notes c2 :onsets '(88)))))))))

(deftest contour-notes-vibrato-across-an-onset
  ;; E then F, each with a 5 Hz vibrato of 80 cents from its own onset at
  ;; frame 25 and 88: the onset places the boundary, not the vibrato.
  (let* ((c (contour-of (lambda (time)
                          (cond ((< time 1/5) nil)
                                ((< time 7/10) (+ 64 (* 4/5 (sin (* 2 pi 5 (- time 1/5))))))
                                ((< time 6/5) (+ 65 (* 4/5 (sin (* 2 pi 5 (- time 7/10))))))))
                        175))
         (ns (contour-notes c :onsets '(88))))
    (ok (equal '((25 88) (88 150)) (spans ns)))
    (ok (equal '(64 65) (mapcar #'round (pitches ns))))))

(deftest contour-notes-minimum-length
  (flet ((runs (frames &rest keys)
           (apply #'contour-notes (contour-of (lambda (time) (and (<= 1/5 time) (< time (+ 1/5 (/ frames 125))) 60))
                                              100)
                  keys)))
    (ok (= 0 (length (runs 12))) "12 frames (96 ms): too short")
    (ok (= 1 (length (runs 13))) "13 frames: M_min = ceiling(0.1 * 125)")
    (ok (= 0 (length (runs 20 :min-duration 1/5))) "M_min 25")
    (ok (= 1 (length (runs 25 :min-duration 1/5))))
    (ok (= 0 (length (runs 10))) "10 frames: too short by default")
    (ok (= 1 (length (runs 10 :min-duration 0))) "but not when M_min is 2")
    (ok (= 0 (length (runs 1 :min-duration 0))) "M_min is at least 2")))

(deftest contour-notes-voicing-and-gaps
  (testing "noise-like voicing stays silent"
    (ok (= 0 (length (contour-notes (contour-of (steps '(nil 1/5) '(60 1) '(nil 1/5)) 175 :voicing 1/100)))))
    (ok (= 1 (length (contour-notes (contour-of (steps '(nil 1/5) '(60 1) '(nil 1/5)) 175 :voicing 1/2))))))
  (testing "an unvoiced frame ends a note"
    (let ((ns (contour-notes (contour-of (steps '(nil 1/5) '(60 1/2) '(nil 1/125) '(60 1/2) '(nil 1/5)) 176))))
      (ok (equal '((25 88) (89 151)) (spans ns)))))
  (testing "pitches outside the range are no notes"
    (ok (= 0 (length (contour-notes (contour-of (steps '(nil 1/5) '(30 1) '(nil 1/5)) 175)))))
    (ok (= 1 (length (contour-notes (contour-of (steps '(nil 1/5) '(30 1) '(nil 1/5)) 175)
                                    :pitch-range '(24 48)))))))

(deftest contour-notes-pitch
  ;; A note at 62 whose first three frames are a leap's glide (55 .. 60) and
  ;; with two frames an octave up: the band leaves them out.
  (let* ((c (contour-of (lambda (time)
                          (let ((i (round (* 125 time))))
                            (cond ((< i 25) nil)
                                  ((< i 28) (+ 55 (* 5/3 (- i 25))))
                                  ((< 60 i 63) 74)
                                  ((< i 100) 62))))
                        110))
         (ns (contour-notes c)))
    (ok (= 1 (length ns)))
    (ok (and (= 1 (length ns)) (close-to (note-event-pitch (svref ns 0)) 62 1d-9))))
  ;; 20 frames at 62 then 30 at 62.4: one note, its pitch the mean 62.24, not
  ;; the median 62.4.
  (let ((ns (contour-notes (contour-of (steps '(nil 1/5) '(62 4/25) '(62.4d0 6/25) '(nil 1/5)) 125))))
    (ok (= 1 (length ns)))
    (ok (and (= 1 (length ns)) (close-to (note-event-pitch (svref ns 0)) 62.24d0 1d-9)) "the mean")))

(deftest contour-notes-arguments
  (let ((c (contour-of (steps '(60 1)) 125)))
    (ok (equalp #() (contour-notes (make-f0-contour))) "no frames")
    (ok (signals (contour-notes '(1 2 3)) 'type-error))
    (ok (signals (contour-notes c :onsets 5) 'type-error))
    (ok (signals (contour-notes c :onsets '(1.5)) 'type-error))
    (ok (signals (contour-notes c :pitch-range '(60 50)) 'type-error))
    (ok (signals (contour-notes c :min-duration -1) 'type-error))
    (ok (signals (contour-notes c :vibrato-window -1) 'type-error))))

(defun smoothed (pitch-at frames &rest keys)
  (apply #'smooth-pitches (contour-of pitch-at frames) keys))

(deftest smooth-pitches-vibrato
  ;; A 5 Hz vibrato of a semitone has a period of 25 frames, the window: away
  ;; from the ends of the run each mean covers one period exactly, while the
  ;; median of a period sampled off its zero crossings is 0.025 off.
  (let ((y (smoothed (lambda (time) (and (< 1/5 time 11/5) (+ 64 (sin (+ (* 2 pi 5 time) 1/10))))) 300)))
    (ok (= 300 (length y)))
    (ok (every #'null (subseq y 0 26)) "unvoiced frames have none")
    (ok (every (lambda (v) (close-to v 64 1/100)) (subseq y 50 250)) "the mean, not the median")
    (ok (every (lambda (v) (typep v 'double-float)) (subseq y 50 250))))
  (let ((y (smoothed (lambda (time) (and (< 1/5 time 11/5) (+ 64 (sin (+ (* 2 pi 5 time) 1/10))))) 300
                     :vibrato-window 0)))
    (ok (close-to (svref y 100) (+ 64 (sin (+ (* 2 pi 5 100/125) 1/10))) 1d-9) "a window of one frame")))

(deftest smooth-pitches-steps-and-outliers
  (testing "a step stays sharp"
    (let ((y (smoothed (steps '(60 2/5) '(62 2/5)) 100)))
      (ok (every (lambda (v) (close-to v 60 1d-9)) (subseq y 0 50)))
      (ok (every (lambda (v) (close-to v 62 1d-9)) (subseq y 50 100)))))
  (testing "an octave error is left out"
    (let ((y (smoothed (lambda (time) (if (< 48/125 time 51/125) 72 60)) 100)))
      (ok (every (lambda (v) (close-to v 60 1d-9)) (concatenate 'list (subseq y 0 49) (subseq y 51 100))))))
  (testing "a step of 0.4 semitone is smoothed, but not across an onset"
    (let ((y (smoothed (steps '(60 2/5) '(60.4d0 2/5)) 100)))
      (ok (< 60.01 (svref y 49) 60.39) "mixed without an onset"))
    (let ((y (smoothed (steps '(60 2/5) '(60.4d0 2/5)) 100 :onsets '(50))))
      (ok (every (lambda (v) (close-to v 60 1d-9)) (subseq y 0 50)))
      (ok (every (lambda (v) (close-to v 60.4d0 1d-9)) (subseq y 50 100)))))
  (testing "argument checks"
    (ok (equalp #() (smooth-pitches (make-f0-contour))))
    (ok (signals (smooth-pitches #(60)) 'type-error))
    (ok (signals (smooth-pitches (contour-of (steps '(60 1)) 10) :onsets '(a)) 'type-error))
    (ok (signals (smooth-pitches (contour-of (steps '(60 1)) 10) :vibrato-window -1) 'type-error))))
