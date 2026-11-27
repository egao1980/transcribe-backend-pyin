;;;; DECODE-PITCH (the pitch and voicing HMM) and PYIN-CONTOUR.

(in-package #:transcribe-backend-pyin/tests)

(defun frames (count &rest pairs)
  "COUNT frames, each with the (bin . probability) PAIRS."
  (make-list count :initial-element (loop for (b p) on pairs by #'cddr collect (cons b p))))

(deftest pitch-bins-and-jumps
  (ok (= 481 (pitch-bin-count '(36 84))) "C2 .. C6 in 10-cent bins")
  (ok (= 321 (pitch-bin-count '(28 60))))
  (ok (= 11 (pitch-bin-count '(69 70))))
  (ok (= 6 (pitch-bin-count '(60 60.55))) "floor(5.5) + 1")
  (ok (= 34 (pitch-max-jump 128 16000)) "3.4 semitones a frame at 8 ms")
  (ok (= 25 (pitch-max-jump 256 44100)) "pYIN's own 25 bins at 256 / 44.1 kHz")
  (ok (= 1 (pitch-max-jump 1 16000)) "at least one bin")
  (dolist (bad '((84 36) (36 36) (-1 60) (36 128) (36) (36 84 90) "36 84" nil))
    (ok (signals (pitch-bin-count bad) 'type-error) (format nil "range ~S" bad)))
  (ok (signals (pitch-max-jump 0 16000) 'type-error))
  (ok (signals (pitch-max-jump 128 16000.0) 'type-error)))

(deftest decode-pitch-paths
  (testing "steady evidence"
    (ok (equalp (make-array 6 :initial-element 10) (decode-pitch (frames 6 10 9/10) 20 3)))
    (ok (equalp (make-array 4 :initial-element nil) (decode-pitch (frames 4) 20 3)) "no candidates")
    (ok (equalp #() (decode-pitch '() 20 3)) "no frames"))
  (testing "a jump beyond the band is no pitch"
    (let ((obs (append (frames 5 10 9/10) (frames 1 18 9/10) (frames 5 10 9/10))))
      (ok (equalp (concatenate 'vector (make-array 5 :initial-element 10) #(nil)
                               (make-array 5 :initial-element 10))
                  (decode-pitch obs 20 3)))
      (ok (eql 18 (svref (decode-pitch obs 20 8) 5)) "within a band of 8 it is followed")))
  (testing "continuity picks among candidates"
    ;; Two candidates a frame, the stronger one at bin 4 but for one frame.
    (let ((obs (append (frames 5 4 1/2 15 2/5) (frames 1 4 2/5 15 1/2) (frames 5 4 1/2 15 2/5))))
      (ok (every (lambda (b) (eql b 4)) (decode-pitch obs 20 3)))))
  (testing "voicing"
    (let ((obs (append (frames 5 10 19/20) (frames 20 10 1/100) (frames 5 10 19/20))))
      (ok (equalp (concatenate 'vector (make-array 5 :initial-element 10)
                               (make-array 20 :initial-element nil)
                               (make-array 5 :initial-element 10))
                  (decode-pitch obs 20 3))
          "a long stretch of low voicing is unvoiced")
      (ok (every (lambda (b) (eql b 10)) (decode-pitch obs 20 3 :switch-probability (expt 10 -12)))
          "unless switching is made nearly impossible")))
  (testing "certain frames (no unvoiced mass) a jump beyond the band apart still decode"
    (let ((path (decode-pitch (append (frames 3 2 1) (frames 3 17 1)) 20 3)))
      (ok (= 6 (length path)))
      (ok (every #'integerp path) "through the probability floor")))
  (testing "ties go to the smaller bin"
    (ok (equalp #(3) (decode-pitch (frames 1 3 1/2 7 1/2) 20 3)))
    (ok (equalp #(3 3) (decode-pitch (frames 2 3 1/2 7 1/2) 20 3))))
  (testing "the band is triangular: a small jump costs more than staying"
    ;; w(0) = 4/16 and w(3) = 1/16 at B = 3: 9/20 * 4/16 beats 1/2 * 1/16.
    (ok (equalp #(10 10) (decode-pitch '(((10 . 1/2)) ((10 . 9/20) (13 . 1/2))) 20 3)))
    (ok (equalp #(10 13) (decode-pitch '(((10 . 1/2)) ((10 . 1/20) (13 . 1/2))) 20 3))))
  (testing "ties between predecessors go to the smaller bin"
    (ok (equalp #(8 10) (decode-pitch '(((8 . 1/2) (12 . 1/2)) ((10 . 1))) 20 3))))
  (testing "argument checks"
    (ok (signals (decode-pitch (frames 2 20 1/2) 20 3) 'type-error) "bin past n - 1")
    (ok (signals (decode-pitch (frames 2 3 -1/2) 20 3) 'type-error) "negative probability")
    (ok (signals (decode-pitch '(((3 . 1/2)) 7) 20 3) 'type-error) "a frame not a list")
    (ok (signals (decode-pitch (frames 2) 0 3) 'type-error))
    (ok (signals (decode-pitch (frames 2) 20 -1) 'type-error))
    (ok (signals (decode-pitch (frames 2) 20 3 :switch-probability 0) 'type-error))
    (ok (signals (decode-pitch (frames 2) 20 3 :switch-probability 1) 'type-error))))

(defun voiced-pitches (contour start end)
  "The MIDI pitches of CONTOUR's voiced frames START .. END - 1."
  (loop for i from start below end
        for f = (svref (f0-contour-f0 contour) i)
        when f collect (hz-to-midi f)))

(deftest-parametrize pyin-contour-steady-tones
    ((waveform tolerance)
     :ids ("sine" "sawtooth")
     (:sine 3/100)
     (:sawtooth 3/100))
  ;; 220 Hz from 0.25 s to 1.25 s, 1.5 s in all.
  (let* ((audio (synthesize-notes (notes '(1/4 5/4 57)) :waveform waveform :duration 3/2))
         (c (pyin-contour audio)))
    (ok (= 16000 (f0-contour-sample-rate c)))
    (ok (= 128 (f0-contour-hop c)))
    (ok (eql 125 (contour-frame-rate c)))
    (ok (= 188 (contour-length c)) "ceiling(24000 / 128)")
    (ok (every #'null (subseq (f0-contour-f0 c) 0 25)) "silence before the note")
    (ok (every #'null (subseq (f0-contour-f0 c) 165)) "and after it")
    (ok (every (lambda (m) (close-to m 57 tolerance)) (voiced-pitches c 36 150)) "within 3 cents")
    (ok (= (- 150 36) (length (voiced-pitches c 36 150))) "voiced throughout")
    (ok (every (lambda (v) (> v 9/10)) (subseq (f0-contour-voicing c) 36 150)))
    (ok (every (lambda (v) (= v 0)) (subseq (f0-contour-voicing c) 0 20)) "no voicing in silence")))

(deftest pyin-contour-silence-gate-and-range
  (testing "silence"
    (let ((c (pyin-contour (silence-buffer 1))))
      (ok (= 125 (contour-length c)))
      (ok (every #'null (f0-contour-f0 c)))
      (ok (every #'zerop (f0-contour-voicing c)))))
  (testing "no samples"
    (ok (= 0 (contour-length (pyin-contour (audio:make-pcm-buffer :sample-rate 16000))))))
  (testing "the level gate"
    ;; A periodic tone 50 dB under a loud one: gated at -40 dB, kept at -60 dB.
    (let* ((loud (synthesize-notes (notes '(1/5 1 57)) :duration 2))
           (quiet (synthesize-notes (notes '(13/10 9/5 64)) :duration 2
                                                           :amplitude (* 1/2 (expt 10d0 -5/2))))
           (audio (buffer (map 'vector #'+ (audio:pcm-buffer-samples loud)
                               (audio:pcm-buffer-samples quiet)))))
      (ok (null (voiced-pitches (pyin-contour audio) 175 215)))
      (ok (every (lambda (m) (close-to m 64 1/20))
                 (voiced-pitches (pyin-contour audio :gate-db -60) 175 215)))
      (ok (= 40 (length (voiced-pitches (pyin-contour audio :gate-db -60) 175 215))))))
  (testing "the gate reads the integration window, not the whole frame"
    ;; A quiet tone, 50 dB down, right before a loud one at 1 s: the frames up
    ;; to 0.968 s reach the loud tone only past their integration window.
    (let* ((quiet (synthesize-notes (notes '(1/5 1 64)) :duration 2 :amplitude (* 1/2 (expt 10d0 -5/2))))
           (loud (synthesize-notes (notes '(1 9/5 57)) :duration 2))
           (audio (buffer (map 'vector #'+ (audio:pcm-buffer-samples quiet)
                               (audio:pcm-buffer-samples loud)))))
      (ok (null (voiced-pitches (pyin-contour audio) 100 122)))))
  (testing "the pitch range"
    ;; A1 (55 Hz) is below C2: no frame tracks it, unless the range includes it.
    (let ((audio (synthesize-notes (notes '(1/4 5/4 33)) :duration 3/2)))
      (ok (notany (lambda (m) (close-to m 33 1/2)) (voiced-pitches (pyin-contour audio) 0 188)))
      (let ((ms (voiced-pitches (pyin-contour audio :pitch-range '(28 60)) 40 150)))
        (ok (= 110 (length ms)))
        (ok (every (lambda (m) (close-to m 33 1/20)) ms))))
    (let ((c (pyin-contour (synthesize-notes (notes '(1/4 5/4 57)) :duration 3/2)
                           :pitch-range '(60 72))))
      (ok (every (lambda (f) (or (null f) (<= 59 (hz-to-midi f) 73))) (f0-contour-f0 c))
          "f0 stays within a semitone of the range (lags are searched inside it)")))
  (testing "any rate when called directly"
    (let ((c (pyin-contour (synthesize-notes (notes '(1/4 5/4 57)) :sample-rate 22050 :duration 3/2))))
      (ok (= 22050 (f0-contour-sample-rate c)))
      (ok (every (lambda (m) (close-to m 57 1/20)) (voiced-pitches c 50 200)))))
  (testing "argument checks"
    (ok (signals (pyin-contour (silence-buffer 1 :channels 2)) 'type-error) "stereo")
    (ok (signals (pyin-contour #(0 0 0)) 'type-error))
    (ok (signals (pyin-contour (silence-buffer 1) :frame-size 1023) 'type-error))
    (ok (signals (pyin-contour (silence-buffer 1) :hop 0) 'type-error))
    (ok (signals (pyin-contour (silence-buffer 1) :gate-db 3) 'type-error))
    (ok (signals (pyin-contour (silence-buffer 1) :pitch-range '(84 36)) 'type-error))))

(defun true-pitch (truth time &key (offset 0) (drift-rate 0) (vibrato-rate 0) (vibrato-depth 0))
  "The pitch SYNTHESIZE-NOTES sings at TIME, or NIL in a rest or within 40 ms of a note's ends."
  (loop for n across truth
        for on = (note-event-onset n)
        for off = (note-event-offset n)
        when (and (>= time (+ on 1/25)) (< time (- off 1/25)))
          return (+ (note-event-pitch n)
                    (/ (+ offset (* drift-rate time)
                          (* vibrato-depth (sin (* 2 pi vibrato-rate (- time on)))))
                       100))))

(defun in-rest-p (truth time)
  "True when TIME lies at least 50 ms from every note."
  (notany (lambda (n) (< (- (note-event-onset n) 1/20) time (+ (note-event-offset n) 1/20))) truth))

(deftest-parametrize pyin-contour-accuracy
    ((waveform offset drift-rate vibrato-depth noise)
     :ids ("sine" "sawtooth" "sine-detuned" "sawtooth-vibrato" "sine-vibrato-noise")
     (:sine 0 0 0 0)
     (:sawtooth 0 0 0 0)
     (:sine -30 -4 0 0)
     (:sawtooth 0 0 60 0)
     (:sine 10 2 60 3/50))
  ;; Frame by frame against the pitch the synthesis sings: the median error of
  ;; the frames voiced in both, voicing recall inside notes, and false voicing
  ;; in rests (with noise the rests are noise).
  (let* ((truth (melody *ode-to-joy*))
         (audio (synthesize-notes truth :waveform waveform :offset offset :drift-rate drift-rate
                                        :vibrato-rate 11/2 :vibrato-depth vibrato-depth
                                        :noise noise :seed 17))
         (c (pyin-contour audio))
         (errors '()) (inside 0) (voiced-inside 0) (rest 0) (voiced-rest 0))
    (dotimes (i (contour-length c))
      (let* ((time (contour-frame-time c i))
             (f (svref (f0-contour-f0 c) i))
             (p (true-pitch truth time :offset offset :drift-rate drift-rate
                                       :vibrato-rate 11/2 :vibrato-depth vibrato-depth)))
        (when p
          (incf inside)
          (when f
            (incf voiced-inside)
            (push (* 100 (abs (- (hz-to-midi f) p))) errors)))
        (when (in-rest-p truth time)
          (incf rest)
          (when f (incf voiced-rest)))))
    (let ((median (median-of errors)))
      (ok (<= median 20) (format nil "median error ~,2F cents" median)))
    (ok (>= (/ voiced-inside inside) 99/100) (format nil "voicing recall ~,3F" (float (/ voiced-inside inside))))
    (ok (<= (/ voiced-rest (max 1 rest)) 1/5) (format nil "false voicing ~,3F" (float (/ voiced-rest (max 1 rest)))))))
