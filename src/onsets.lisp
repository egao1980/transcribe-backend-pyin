;;;; T1 onsets: where a note is sung again at the same pitch.
;;;;
;;;; A pitch track can't split a repeated pitch: the f0 stays put. The singer's
;;;; re-articulation can: the level dips between the two notes (a consonant, a
;;;; breath, a fresh attack) and rises again. The cue here is such a dip of a
;;;; short-time level envelope, finer than the pitch frames. A level that only
;;;; rises (the first note after a rest) or only falls (a note's end) is not a
;;;; dip: voicing handles those.
;;;;
;;;; Data:
;;;;   +DIP-DB+    6: a dip is at least 6 dB below the loudest level on each side
;;;;   +DIP-SPAN+  1/20: seconds searched on each side of a dip
;;;;
;;;; (ONSET-FRAMES audio &key (hop 128)) -> a fresh list of frame indices, in
;;;;   increasing order without repeats: the frames where a dip lies. AUDIO is
;;;;   a one-channel PCM-BUFFER and HOP an integer >= 1 (CHECK-TYPE, in that
;;;;   order); frames are the contour's, frame i at i * HOP / s.
;;;;   1. Envelope. q = max(1, floor(HOP / 4)) samples (2 ms at hop 128 and
;;;;      16 kHz) and r = max(1, floor(HOP / 2)). With L samples, the envelope
;;;;      has J = ceiling(L / q) steps (0 for L = 0): e_j is the RMS, in
;;;;      DOUBLE-FLOAT, of the 2r samples x_{jq - r} .. x_{jq + r - 1}, 0 outside
;;;;      the audio (a window of one hop, 8 ms, centred on sample jq).
;;;;   2. Dips. D = max(1, round(+DIP-SPAN+ * s / q)) steps (25 at the
;;;;      defaults). Step j, with L_j = e_{max(0, j-D)} .. e_{j-1} and R_j =
;;;;      e_{j+1} .. e_{min(J-1, j+D)} both non-empty, is a dip when
;;;;        e_j < every value of L_j and e_j <= every value of R_j (the first
;;;;        step of a flat bottom), and
;;;;        max L_j > G e_j and max R_j > G e_j, G = (db-to-gain +DIP-DB+).
;;;;   3. Frames. Dip j lies in frame floor((j q + floor(HOP / 2)) / HOP), the
;;;;      frame whose time is nearest to sample jq (the later one of two).
;;;;
;;;; The synthesized takes of transcribe-protocol/testing ramp each note in
;;;; and out over 10 ms, so two notes that touch dip to near 0 between them:
;;;; more than 10 dB on this envelope, which a steady tone (sine, sawtooth,
;;;; with vibrato, or under noise 20 dB down) never shows.

(in-package #:transcribe-backend-pyin)

(defconstant +dip-db+ 6
  "A re-articulation dip lies at least this many dB below the level on each side.")

(defconstant +dip-span+ 1/20
  "Seconds searched on each side of a dip.")

(defun %dip-p (e j span g)
  "True when step J of the envelope E is a dip: E_j below every value of the
SPAN steps before it and at most every value of the SPAN steps after it (both
sides clipped to E, and not empty), and both sides' loudest values above G
times E_j."
  (declare (type (simple-array double-float (*)) e) (type fixnum j span)
           (type double-float g) (optimize speed))
  (let ((ej (aref e j))
        (lmax 0d0)
        (rmax 0d0))
    (declare (type double-float ej lmax rmax))
    (loop for k of-type fixnum from (max 0 (- j span)) below j
          do (let ((v (aref e k)))
               (unless (< ej v) (return-from %dip-p nil))
               (when (> v lmax) (setf lmax v))))
    (loop for k of-type fixnum from (1+ j) to (min (1- (length e)) (+ j span))
          do (let ((v (aref e k)))
               (unless (<= ej v) (return-from %dip-p nil))
               (when (> v rmax) (setf rmax v))))
    (let ((level (* g ej)))
      (and (> lmax level) (> rmax level)))))

(defun onset-frames (audio &key (hop 128))
  "The frames of AUDIO (a mono PCM-BUFFER) where its level dips between two
louder stretches: re-articulations. See onsets.lisp."
  (check-type audio (satisfies %mono-buffer-p) "a one-channel PCM-BUFFER")
  (check-type hop (integer 1))
  (let* ((s (audio:pcm-buffer-sample-rate audio))
         (q (max 1 (floor hop 4)))
         (r (max 1 (floor hop 2)))
         ;; e_j: the RMS of the 2r samples from jq - r, zeros outside the
         ;; audio; ceiling(L / q) steps.
         (e (dsp:frame-rms audio (* 2 r) q))
         (steps (length e))
         (span (max 1 (round (* +dip-span+ s) q)))
         (g (audio:db-to-gain +dip-db+))
         (half (floor hop 2))
         (frames '()))
    (declare (type (simple-array double-float (*)) e))
    ;; Both sides non-empty: steps 1 .. J - 2. A dip's own neighbours are
    ;; tested first, as most steps fail there.
    (loop for j from 1 below (1- steps)
          do (let ((ej (aref e j)))
               (when (and (< ej (aref e (1- j)))
                          (<= ej (aref e (1+ j)))
                          (%dip-p e j span g))
                 (let ((frame (floor (+ (* j q) half) hop)))
                   ;; J rises, so the frames never fall: drop repeats.
                   (unless (eql frame (first frames))
                     (push frame frames))))))
    (nreverse frames)))
