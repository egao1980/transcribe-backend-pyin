;;;; PYIN-TRANSCRIBER end to end through TRANSCRIBE, on synthesized takes with
;;;; known notes, against the YIN baseline where the baseline fails.

(in-package #:transcribe-backend-pyin/tests)

(defun sung (spec &rest synthesis)
  "(values truth audio) for the melody SPEC sung as SYNTHESIS says."
  (let ((truth (melody spec)))
    (values truth (apply #'synthesize-notes truth synthesis))))

(defun within-cents (truth tr cents)
  "True when every matched note of TR's corrected notes is within CENTS of TRUTH."
  (<= (worst-cents truth (transcription-notes tr)) cents))

(defparameter *clean-takes*
  ;; (melody waveform)
  `((,*ode-to-joy* :sine) (,*ode-to-joy* :sawtooth)
    (,*twinkle* :sine) (,*twinkle* :sawtooth)
    (,*frere-jacques* :sine) (,*frere-jacques* :sawtooth)))

(deftest-parametrize takes-melodies
    ((spec waveform)
     :ids '("ode-sine" "ode-sawtooth" "twinkle-sine" "twinkle-sawtooth" "frere-sine" "frere-sawtooth")
     :rows *clean-takes*)
  (multiple-value-bind (truth audio) (sung spec :waveform waveform)
    (let* ((tr (transcribe (pyin) audio :hints '(:key "C")))
           (f (f-of truth (transcription-notes tr))))
      (ok (>= f 19/20) (format nil "note F ~,3F" f))
      (ok (within-cents truth tr 10) (format nil "worst ~,1F cents" (worst-cents truth (transcription-notes tr)))))))

(deftest-parametrize takes-legato
    ((spec waveform)
     :ids ("scale-sine" "scale-sawtooth" "chromatic-sine" "chromatic-sawtooth")
     (*scale* :sine) (*scale* :sawtooth) (*chromatic* :sine) (*chromatic* :sawtooth))
  ;; Notes that touch without a ramp: no dips, so the pitch alone splits them.
  (multiple-value-bind (truth audio) (sung spec :waveform waveform :attack 0 :release 0)
    (let ((f (f-of truth (transcription-notes (transcribe (pyin) audio)))))
      (ok (>= f 9/10) (format nil "note F ~,3F" f)))))

(deftest-parametrize takes-vibrato
    ((waveform depth rate f-min)
     :ids ("sine-60c" "sawtooth-60c" "sine-100c" "sawtooth-100c")
     (:sine 60 11/2 9/10) (:sawtooth 60 11/2 9/10) (:sine 100 11/2 9/10) (:sawtooth 100 11/2 9/10))
  ;; Deeper than the baseline's 50 cents: YIN splits the notes, pYIN keeps them.
  (multiple-value-bind (truth audio) (sung *ode-to-joy* :waveform waveform
                                                        :vibrato-rate rate :vibrato-depth depth)
    (let ((f (f-of truth (transcription-notes (transcribe (pyin) audio))))
          (baseline (f-of truth (transcription-notes (transcribe (yin) audio)))))
      (ok (>= f f-min) (format nil "note F ~,3F" f))
      (ok (>= f (+ baseline 2/5)) (format nil "pYIN ~,3F, YIN ~,3F" f baseline)))))

(deftest-parametrize takes-one-long-note
    ((depth)
     :ids ("60c" "100c")
     (60) (100))
  ;; Two seconds of E4 with a 5.5 Hz vibrato: one note, where YIN gives many.
  (let* ((truth (notes '(1/5 11/5 64)))
         (audio (synthesize-notes truth :vibrato-rate 11/2 :vibrato-depth depth))
         (ns (transcription-notes (transcribe (pyin) audio)))
         (baseline (transcription-notes (transcribe (yin) audio))))
    (ok (= 1 (length ns)) (format nil "~D notes" (length ns)))
    (ok (and (= 1 (length ns)) (close-to (note-event-pitch (svref ns 0)) 64 1/10)))
    (ok (> (length baseline) 3) (format nil "YIN: ~D notes" (length baseline)))))

(deftest-parametrize takes-repeated-pitches
    ((waveform)
     :ids ("sine" "sawtooth")
     (:sine) (:sawtooth))
  ;; Six repeated pairs: the baseline merges each pair, pYIN splits them at the dip.
  (multiple-value-bind (truth audio) (sung *twinkle* :waveform waveform)
    (let ((f (f-of truth (transcription-notes (transcribe (pyin) audio))))
          (baseline (f-of truth (transcription-notes (transcribe (yin) audio)))))
      (ok (>= f 19/20) (format nil "note F ~,3F" f))
      (ok (<= baseline 3/4) (format nil "YIN ~,3F" baseline))
      (ok (>= f (+ baseline 1/5))))))

(deftest-parametrize takes-drift-and-offset
    ((waveform offset drift-rate)
     :ids ("sine-flat-drifting" "sawtooth-sharp-drifting" "sine-flat")
     (:sine 10 -4) (:sawtooth -10 4) (:sine -35 0))
  ;; FINISH-TRANSCRIPTION removes the tuning and the per-phrase drift: 64
  ;; cents over the take, 32 within each phrase, the take's mean within 25
  ;; cents of A440 (a mean near 50 cents would read a semitone away).
  (multiple-value-bind (truth audio) (sung *ode-to-joy* :waveform waveform :offset offset
                                                        :drift-rate drift-rate)
    (let ((tr (transcribe (pyin) audio :hints '(:key "C"))))
      (ok (>= (f-of truth (transcription-notes tr)) 19/20))
      (ok (within-cents truth tr 20) (format nil "worst ~,1F cents after the correction"
                                             (worst-cents truth (transcription-notes tr))))
      (ok (close-to (transcription-tuning-offset tr)
                    (estimate-tuning (detune-notes truth :offset offset :drift-rate drift-rate))
                    3)
          (format nil "tuning ~D" (transcription-tuning-offset tr)))
      (ok (= 2 (length (transcription-drift tr))) "two phrases"))))

(deftest-parametrize takes-noise
    ((waveform spec)
     :ids ("ode-sine" "ode-sawtooth" "frere-sine")
     (:sine *ode-to-joy*) (:sawtooth *ode-to-joy*) (:sine *frere-jacques*))
  ;; Uniform noise of 3/50: 20 dB under the sine's power (1/8 against 3/2500).
  (multiple-value-bind (truth audio) (sung spec :waveform waveform :noise 3/50 :seed 17)
    (let ((f (f-of truth (transcription-notes (transcribe (pyin) audio)))))
      (ok (>= f 9/10) (format nil "note F ~,3F" f)))))

(deftest takes-silence-and-noise
  (let ((tr (handler-bind ((silent-audio #'continue))
              (transcribe (pyin) (silence-buffer 1)))))
    (ok (equalp #() (transcription-notes tr)))
    (ok (every #'null (f0-contour-f0 (transcription-contour tr))) "every frame unvoiced"))
  (let ((tr (transcribe (pyin) (synthesize-notes #() :noise 3/10 :duration 2 :seed 5))))
    (ok (equalp #() (transcription-notes tr)) "white noise has no notes"))
  (multiple-value-bind (ns contour) (detect-notes (pyin) (audio:make-pcm-buffer :sample-rate 16000) nil)
    (ok (= 0 (length ns)) "no samples")
    (ok (= 0 (contour-length contour)))))

(deftest takes-minimum-length
  ;; Between rests, a note of 60 ms is dropped and one of 160 ms kept; with
  ;; :MIN-DURATION 1/4, the 160 ms note is dropped and one of 300 ms kept.
  (flet ((count-notes (seconds &rest initargs)
           (length (transcription-notes
                    (transcribe (apply #'pyin initargs)
                                (synthesize-notes (notes (list 1/5 (+ 1/5 seconds) 64)) :duration 3/4))))))
    (ok (= 0 (count-notes 3/50)) "60 ms")
    (ok (= 1 (count-notes 4/25)) "160 ms")
    (ok (= 0 (count-notes 4/25 :min-duration 1/4)) "160 ms, minimum 1/4 s")
    (ok (= 1 (count-notes 3/10 :min-duration 1/4)) "300 ms, minimum 1/4 s")))

(deftest takes-pitch-range
  ;; A1 (55 Hz) lies below the default range, C2 .. C6.
  (let ((audio (synthesize-notes (notes '(1/4 5/4 33)) :duration 3/2)))
    (ok (equalp #() (transcription-notes (transcribe (pyin) audio))))
    (let ((ns (transcription-raw-notes (transcribe (pyin) audio :hints '(:pitch-range (28 60))))))
      (ok (= 1 (length ns)))
      (ok (and (= 1 (length ns)) (close-to (note-event-pitch (svref ns 0)) 33 1/20))))))

(deftest takes-hint-distance
  (dolist (case '((0 0) (-100 -1) (100 1) (-30 0)))
    (destructuring-bind (offset distance) case
      (multiple-value-bind (truth audio) (sung *frere-jacques* :waveform :sawtooth :offset offset)
        (declare (ignore truth))
        (let ((tr (transcribe (pyin) audio :hints '(:key "C"))))
          (ok (eql distance (mt:key-fit-hint-distance (transcription-key-fit tr)))
              (format nil "sung ~@D cents: distance ~@D" offset distance)))))))

(deftest takes-protocol
  (multiple-value-bind (truth audio) (sung (subseq *twinkle* 0 7) :waveform :sawtooth :sample-rate 22050)
    (let ((tr (transcribe (pyin) audio)))
      (ok (= 16000 (f0-contour-sample-rate (transcription-contour tr))) "resampled to 16 kHz")
      (ok (eql 125 (contour-frame-rate (transcription-contour tr))))
      (ok (>= (f-of truth (transcription-notes tr)) 19/20))
      (ok (equal (transcriber-info (pyin)) (transcription-backend tr)))))
  (let* ((mono (audio:pcm-buffer-samples (synthesize-notes (notes '(1/5 1 57)))))
         (stereo (audio:pcm-from-channels (list mono mono) :sample-rate 16000))
         (tr (transcribe (pyin) stereo)))
    (ok (= 1 (length (transcription-notes tr))) "a stereo take is downmixed by the protocol")))
