;;;; T1 transcriber: PYIN-TRANSCRIBER, a subclass of transcribe-protocol's
;;;; TRANSCRIBER.
;;;;
;;;; Initargs and readers (PYIN-TRANSCRIBER-...), checked at MAKE-INSTANCE with
;;;; CHECK-TYPE (TYPE-ERROR, STORE-VALUE), in this order:
;;;;   :FRAME-SIZE      N, an even integer >= 4, default 1024 (64 ms at 16 kHz).
;;;;                    YIN's integration window is W = N / 2 (32 ms), and lags
;;;;                    run up to W.
;;;;   :HOP             an integer >= 1, default 128 (8 ms: 125 frames a second)
;;;;   :GATE-DB         a REAL <= 0, default -40: frames whose integration window
;;;;                    is quieter than this, in dB below the loudest frame's,
;;;;                    have no pitch candidates
;;;;   :MIN-DURATION    seconds, a non-negative REAL, default 1/10: the shortest
;;;;                    note CONTOUR-NOTES gives (Tony's pruning threshold)
;;;;   :VIBRATO-WINDOW  seconds, a non-negative REAL, default 1/5: the span of
;;;;                    the moving mean that smooths vibrato before the note HMM
;;;; Frames of 1024 every 128 at 16 kHz are pYIN's own time scale (Mauch and
;;;; Dixon used 2048 every 256 at 44.1 kHz: 46 and 5.8 ms) and the YIN
;;;; baseline's, so the two compare frame by frame. librosa's 2048 at
;;;; 22.05 kHz would be 128 ms at 16 kHz and smear onsets by twice as much.
;;;;
;;;; TRANSCRIBER-SAMPLE-RATE is 16000 (the protocol's default): TRANSCRIBE
;;;; resamples to it. DETECT-NOTES itself works at any rate s of the buffer it
;;;; is given; the sample counts above are then samples at s.
;;;;
;;;; (TRANSCRIBER-INFO p) -> (:name "pyin-transcriber" :frame-size N :hop H
;;;;   :gate-db G :min-duration D :vibrato-window V), a fresh list of the
;;;;   initargs' values: the backend's provenance in each TRANSCRIPTION.

(in-package #:transcribe-backend-pyin)

;;; Predicates for the argument checks of every file.

(defun %proper-list-p (x)
  (handler-case (list-length x) (type-error () nil)))

(defun %pitch-range-p (x)
  "True for a list (low high) of reals with 0 <= low < high < 128."
  (and (consp x) (consp (cdr x)) (null (cddr x))
       (realp (first x)) (realp (second x))
       (<= 0 (first x)) (< (first x) (second x)) (< (second x) 128)))

(defun %mono-buffer-p (x)
  "True for a one-channel audio-protocol PCM-BUFFER."
  (and (typep x 'audio:pcm-buffer) (= 1 (audio:pcm-buffer-channels x))))

(defclass pyin-transcriber (tp:transcriber)
  ((frame-size :initarg :frame-size :initform 1024 :reader pyin-transcriber-frame-size)
   (hop :initarg :hop :initform 128 :reader pyin-transcriber-hop)
   (gate-db :initarg :gate-db :initform -40 :reader pyin-transcriber-gate-db)
   (min-duration :initarg :min-duration :initform 1/10 :reader pyin-transcriber-min-duration)
   (vibrato-window :initarg :vibrato-window :initform 1/5
                   :reader pyin-transcriber-vibrato-window))
  (:documentation "pYIN pitch tracking and HMM note segmentation: a
transcribe-protocol TRANSCRIBER. See transcriber.lisp."))

(defmethod initialize-instance :after ((p pyin-transcriber) &key)
  (with-slots (frame-size hop gate-db min-duration vibrato-window) p
    (check-type frame-size (and (integer 4) (satisfies evenp)) "an even integer >= 4")
    (check-type hop (integer 1))
    (check-type gate-db (real * 0))
    (check-type min-duration (real 0))
    (check-type vibrato-window (real 0))))

(defmethod tp:transcriber-info ((p pyin-transcriber))
  (list :name "pyin-transcriber"
        :frame-size (pyin-transcriber-frame-size p)
        :hop (pyin-transcriber-hop p)
        :gate-db (pyin-transcriber-gate-db p)
        :min-duration (pyin-transcriber-min-duration p)
        :vibrato-window (pyin-transcriber-vibrato-window p)))
