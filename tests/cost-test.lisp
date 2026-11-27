;;;; The cost of a long take, as ratios measured on one machine (so they hold
;;;; on a slow CI runner too): against the YIN baseline on the same take, and
;;;; against a take six times shorter.

(in-package #:transcribe-backend-pyin/tests)

(defun seconds-of (thunk)
  "Wall-clock seconds THUNK takes, as a double."
  (let ((start (get-internal-real-time)))
    (funcall thunk)
    (/ (float (- (get-internal-real-time) start) 1d0) internal-time-units-per-second)))

(defun long-take (repeats)
  "Ode to Joy sung REPEATS times in a row, sine at 16 kHz: about 16.6 s a time."
  (let* ((once (melody *ode-to-joy*))
         (span (+ (note-event-offset (svref once (1- (length once)))) 2/5))
         (all (loop for r below repeats
                    append (map 'list (lambda (n)
                                        (note (+ (note-event-onset n) (* r span))
                                              (+ (note-event-offset n) (* r span))
                                              (note-event-pitch n)))
                                once))))
    (synthesize-notes (coerce all 'vector))))

(deftest cost-three-minutes
  ;; 11 times 16.6 s: three minutes (2.9 million samples). pYIN must stay
  ;; within 8 times the YIN baseline on the same take (about 2.5 measured),
  ;; and grow linearly: within 8 times a take of 2 repeats (about 6 expected).
  (let* ((long (long-take 11))
         (short (long-take 2))
         (pyin-long 0d0) (yin-long 0d0) (pyin-short 0d0) (notes-long 0))
    (ok (>= (audio:pcm-buffer-duration long) 180) "three minutes")
    (setf pyin-long (seconds-of (lambda () (setf notes-long (length (transcription-notes
                                                                      (transcribe (pyin) long)))))))
    (setf yin-long (seconds-of (lambda () (transcribe (yin) long))))
    (setf pyin-short (min (seconds-of (lambda () (transcribe (pyin) short)))
                          (seconds-of (lambda () (transcribe (pyin) short)))))
    (ok (= (* 11 30) notes-long) (format nil "~D notes" notes-long))
    (ok (< pyin-long (* 8 yin-long))
        (format nil "pYIN ~,2F s, YIN ~,2F s: ratio ~,2F" pyin-long yin-long (/ pyin-long yin-long)))
    (ok (< pyin-long (* 8 pyin-short))
        (format nil "3 min ~,2F s, 33 s ~,2F s: ratio ~,2F" pyin-long pyin-short (/ pyin-long pyin-short)))))
