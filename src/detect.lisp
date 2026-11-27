;;;; T1 detection: the DETECT-NOTES method of PYIN-TRANSCRIBER.
;;;;
;;;; (DETECT-NOTES p audio hints) -> (values notes contour), for AUDIO the mono
;;;;   buffer PREPARE-AUDIO gives (16 kHz through TRANSCRIBE; any rate s and
;;;;   any length, 0 samples included, when called directly) and HINTS a
;;;;   normalized hint list:
;;;;     range   = (tp:hint-pitch-range hints)
;;;;     contour = (pyin-contour audio :frame-size N :hop H :gate-db G
;;;;                                   :pitch-range range)
;;;;     onsets  = (onset-frames audio :hop H)
;;;;     notes   = (contour-notes contour :onsets onsets :pitch-range range
;;;;                              :min-duration D :vibrato-window V)
;;;;   with N, H, G, D, V the transcriber's initargs. NOTES are as sung: the
;;;;   protocol's FINISH-TRANSCRIPTION estimates and removes the tuning and the
;;;;   drift and fits the key; nothing of that happens here. A take of one
;;;;   channel needs no preparation; TRANSCRIBE downmixes and resamples others.

(in-package #:transcribe-backend-pyin)

(defmethod tp:detect-notes ((p pyin-transcriber) audio hints)
  (let* ((range (tp:hint-pitch-range hints))
         (hop (pyin-transcriber-hop p))
         (contour (pyin-contour audio :frame-size (pyin-transcriber-frame-size p)
                                      :hop hop
                                      :gate-db (pyin-transcriber-gate-db p)
                                      :pitch-range range))
         (onsets (onset-frames audio :hop hop))
         (notes (contour-notes contour :onsets onsets
                                       :pitch-range range
                                       :min-duration (pyin-transcriber-min-duration p)
                                       :vibrato-window (pyin-transcriber-vibrato-window p))))
    (values notes contour)))
