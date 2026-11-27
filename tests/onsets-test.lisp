;;;; ONSET-FRAMES: dips of the level envelope between two louder stretches.

(in-package #:transcribe-backend-pyin/tests)

(defun boundary-frames (truth)
  "The frames (at 8 ms) where a note starts as the one before it ends."
  (loop for (a b) on (coerce truth 'list)
        while b
        when (= (note-event-offset a) (note-event-onset b))
          collect (round (* 125 (note-event-onset b)))))

(defun near-each (expected actual tolerance)
  "True when ACTUAL has one element within TOLERANCE of each of EXPECTED, and no other."
  (and (= (length expected) (length actual))
       (every (lambda (e a) (<= (abs (- e a)) tolerance)) expected actual)))

(deftest onsets-two-notes
  ;; The same pitch twice, touching at 0.6 s (frame 75), each ramped over 10 ms.
  (let ((audio (synthesize-notes (notes '(1/5 3/5 60) '(3/5 1 60)))))
    (ok (equal '(75) (onset-frames audio)))
    (ok (equal '(38) (onset-frames audio :hop 256)) "frames of another hop: 37.5, the later")))

(deftest-parametrize onsets-in-takes
    ((name waveform vibrato-depth noise)
     :ids ("sine" "sawtooth" "vibrato" "noise")
     ("sine" :sine 0 0)
     ("sawtooth" :sawtooth 0 0)
     ("vibrato" :sine 100 0)
     ("noise" :sawtooth 0 3/50))
  (declare (ignore name))
  ;; Every touching pair of notes has its dip, within a frame; nothing else does.
  (let* ((truth (melody *ode-to-joy*))
         (audio (synthesize-notes truth :waveform waveform :vibrato-rate 11/2
                                        :vibrato-depth vibrato-depth :noise noise :seed 17)))
    (ok (near-each (boundary-frames truth) (onset-frames audio) 1)
        (format nil "~D boundaries" (length (boundary-frames truth))))))

(deftest onsets-not-dips
  (testing "legato: notes that touch without a ramp"
    (ok (null (onset-frames (synthesize-notes (melody *scale*) :attack 0 :release 0)))))
  (testing "a note after a rest, and a note's end: a rise or a fall only"
    (ok (null (onset-frames (synthesize-notes (notes '(1/5 3/5 60) '(1 7/5 64)))))))
  (testing "a crescendo"
    (ok (null (onset-frames (sine-buffer 220 1 :envelope (lambda (time) (+ 1/10 (* 9/10 time))))))))
  (testing "a dip of 4 dB is too shallow, one of 12 dB is not"
    (flet ((dipped (gain)
             (sine-buffer 220 1 :envelope (lambda (time) (if (< 49/100 time 51/100) gain 1)))))
      (ok (null (onset-frames (dipped (expt 10 -4/20)))))
      (let ((frames (onset-frames (dipped (expt 10 -12/20)))))
        (ok (and (= 1 (length frames)) (member (first frames) '(62 63))) "the dip at 0.49 .. 0.51 s"))))
  (testing "silence and no samples"
    (ok (null (onset-frames (silence-buffer 1))))
    (ok (null (onset-frames (audio:make-pcm-buffer :sample-rate 16000)))))
  (testing "argument checks"
    (ok (signals (onset-frames (silence-buffer 1 :channels 2)) 'type-error))
    (ok (signals (onset-frames (silence-buffer 1) :hop 0) 'type-error))))
