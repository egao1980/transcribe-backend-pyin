;;;; Helpers and test melodies shared by the suite.

(in-package #:transcribe-backend-pyin/tests)

(defun pyin (&rest initargs)
  (apply #'make-instance 'pyin-transcriber initargs))

(defun yin (&rest initargs)
  (apply #'make-instance 'yin-transcriber initargs))

(defun note (onset offset pitch &optional (confidence 1))
  (make-note-event :onset onset :offset offset :pitch pitch :confidence confidence))

(defun notes (&rest specs)
  "A SIMPLE-VECTOR of NOTE-EVENTs from (onset offset pitch [confidence]) lists."
  (map 'simple-vector (lambda (spec) (apply #'note spec)) specs))

(defun close-to (a b tolerance)
  (<= (abs (- a b)) tolerance))

(defun pitches (notes)
  (map 'list #'note-event-pitch notes))

(defun buffer (samples &key (rate 16000) (channels 1))
  (audio:make-pcm-buffer :sample-rate rate :channels channels :samples samples))

(defun sine-samples (frequency seconds &key (rate 16000) (amplitude 0.5d0) (envelope (constantly 1)))
  "A sine of FREQUENCY Hz for SECONDS, its amplitude times (ENVELOPE time)."
  (let* ((n (round (* seconds rate)))
         (v (make-array n :element-type 'double-float)))
    (dotimes (i n v)
      (setf (aref v i) (* amplitude (funcall envelope (/ i rate))
                          (sin (/ (* 2 (float pi 1d0) frequency i) rate)))))))

(defun sine-buffer (frequency seconds &rest keys &key (rate 16000) &allow-other-keys)
  (buffer (apply #'sine-samples frequency seconds keys) :rate rate))

(defun silence-buffer (seconds &key (rate 16000) (channels 1))
  (audio:make-pcm-buffer :sample-rate rate :channels channels
                         :frames (round (* seconds rate))))

(defun caught (thunk type)
  "The condition of TYPE that THUNK signals with ERROR, unwound; NIL if none."
  (handler-case (progn (funcall thunk) nil)
    (error (c) (and (typep c type) c))))

(defun f-of (reference estimate)
  "Onset-only note F-measure: onsets within 50 ms, pitches within 50 cents."
  (values (note-f-measure reference estimate)))

(defun worst-cents (reference estimate)
  "The largest pitch distance, in cents, over a maximum matching of the two."
  (loop for (i . j) in (match-notes reference estimate)
        maximize (* 100 (abs (- (note-event-pitch (svref reference i))
                                (note-event-pitch (svref estimate j)))))))

(defun median-of (values)
  (let* ((sorted (sort (copy-list values) #'<))
         (n (length sorted)))
    (if (oddp n)
        (nth (floor n 2) sorted)
        (/ (+ (nth (1- (floor n 2)) sorted) (nth (floor n 2) sorted)) 2))))

(defun contour-of (pitch-at frames &key (rate 16000) (hop 128) (voicing 1))
  "An F0-CONTOUR of FRAMES frames: frame i has the pitch (PITCH-AT time), a
MIDI number or NIL for unvoiced, at time i * HOP / RATE; VOICING is a number
or a function of time."
  (let ((f0 (make-array frames :initial-element nil))
        (v (make-array frames :initial-element 0)))
    (dotimes (i frames)
      (let* ((time (/ (* i hop) rate))
             (pitch (funcall pitch-at time)))
        (when pitch
          (setf (svref f0 i) (midi-to-hz pitch)
                (svref v i) (if (functionp voicing) (funcall voicing time) voicing)))))
    (make-f0-contour :sample-rate rate :hop hop :f0 f0 :voicing v)))

(defun steps (&rest spec)
  "A function of time for CONTOUR-OF from (pitch seconds) pieces, pitch NIL for unvoiced."
  (lambda (time)
    (loop for (pitch seconds) in spec
          for end = seconds then (+ end seconds)
          when (< time end) return pitch)))

;;; Melodies --------------------------------------------------------------------
;;; Public-domain tunes as (pitch seconds) entries for MELODY-NOTES, and two
;;; figures written for these tests. Note values: quarter 1/2 s unless stated.

(defparameter *ode-to-joy*
  '((64 1/2) (64 1/2) (65 1/2) (67 1/2) (67 1/2) (65 1/2) (64 1/2) (62 1/2)
    (60 1/2) (60 1/2) (62 1/2) (64 1/2) (64 3/4) (62 1/4) (62 3/4) (nil 3/10)
    (64 1/2) (64 1/2) (65 1/2) (67 1/2) (67 1/2) (65 1/2) (64 1/2) (62 1/2)
    (60 1/2) (60 1/2) (62 1/2) (64 1/2) (62 3/4) (60 1/4) (60 1))
  "Beethoven's Ode to Joy (1824), two phrases in C with a breath: 30 notes,
eleven of them repeating the pitch before them, eighths of 1/4 s.")

(defparameter *twinkle*
  '((60 2/5) (60 2/5) (67 2/5) (67 2/5) (69 2/5) (69 2/5) (67 4/5)
    (65 2/5) (65 2/5) (64 2/5) (64 2/5) (62 2/5) (62 2/5) (60 4/5))
  "The French tune Ah ! vous dirai-je, maman (1761), as Twinkle, Twinkle: 14
notes in six repeated pairs.")

(defparameter *frere-jacques*
  '((60 1/2) (62 1/2) (64 1/2) (60 1/2) (60 1/2) (62 1/2) (64 1/2) (60 1/2)
    (64 1/2) (65 1/2) (67 1) (64 1/2) (65 1/2) (67 1) (nil 3/10)
    (67 1/4) (69 1/4) (67 1/4) (65 1/4) (64 1/2) (60 1/2)
    (67 1/4) (69 1/4) (67 1/4) (65 1/4) (64 1/2) (60 1/2)
    (60 1/2) (55 1/2) (60 1) (60 1/2) (55 1/2) (60 1))
  "The round Frère Jacques (18th century) in C: 32 notes, eighths of 1/4 s.")

(defparameter *scale*
  '((60 1/4) (62 1/4) (64 1/4) (65 1/4) (67 1/4) (69 1/4) (71 1/4) (72 1/4)
    (71 1/4) (69 1/4) (67 1/4) (65 1/4) (64 1/4) (62 1/4) (60 1/4))
  "C major up and down in eighths of 1/4 s, no pitch repeated: for legato takes.")

(defparameter *chromatic*
  '((60 1/4) (61 1/4) (62 1/4) (63 1/4) (64 1/4) (65 1/4) (66 1/4) (67 1/4)
    (66 1/4) (65 1/4) (64 1/4) (63 1/4) (62 1/4) (61 1/4) (60 1/4))
  "Semitone steps up and down, 1/4 s each.")

(defun melody (spec &rest keys)
  (apply #'melody-notes spec :start 1/5 keys))
