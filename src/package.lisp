;;;; transcribe-backend-pyin: pYIN pitch tracking and HMM note segmentation as
;;;; a transcribe-protocol backend. Brief: cl-stack/docs/capabilities/transcribe.md,
;;;; slice T1 (0.1.0).
;;;;
;;;; One sung or hummed line at 16 kHz in; the notes as sung and the f0
;;;; contour out, through transcribe-protocol's DETECT-NOTES. Everything after
;;;; that (the take's tuning, its drift per phrase, the key fit, the checks of
;;;; the audio and the hints) is the protocol's, and this backend never
;;;; computes it: TRANSCRIBE does, in FINISH-TRANSCRIPTION.
;;;;
;;;; The pipeline (each file header is the exact spec of its part):
;;;;   transcriber.lisp  PYIN-TRANSCRIBER: initargs, their checks, TRANSCRIBER-INFO
;;;;   candidates.lisp   pYIN's threshold distribution: THRESHOLD-PRIOR (a beta
;;;;                     prior over 100 thresholds) and PITCH-CANDIDATES (the
;;;;                     periods of one frame's d' troughs, with probabilities)
;;;;   pitch.lisp        DECODE-PITCH (a pitch and voicing HMM over 10-cent
;;;;                     bins, decoded by Viterbi) and PYIN-CONTOUR (frames,
;;;;                     level gate, candidates, decoding -> F0-CONTOUR)
;;;;   onsets.lisp       ONSET-FRAMES: re-articulations, as dips of a short-time
;;;;                     level envelope
;;;;   notes.lisp        SMOOTH-PITCHES (vibrato smoothing) and CONTOUR-NOTES (a
;;;;                     note HMM: attack and stable states per 1/3 semitone, one
;;;;                     silent state, decoded by Viterbi -> NOTE-EVENTs)
;;;;   detect.lisp       the DETECT-NOTES method: PYIN-CONTOUR, ONSET-FRAMES,
;;;;                     CONTOUR-NOTES
;;;;
;;;; Prior art: YIN (de Cheveigné and Kawahara 2002); pYIN (Mauch and Dixon,
;;;; ICASSP 2014) for the thresholds and the pitch HMM; Tony (Mauch et al.,
;;;; TENOR 2015) for the note HMM; librosa 0.10's pyin for conventions (beta
;;;; (2, 18), 100 thresholds, 10-cent bins, a switch probability of 0.01, a
;;;; transition rate of 35.92 octaves a second). The Vamp plugin and Tony are
;;;; GPL and were not used as code.
;;;;
;;;; Conditions: transcribe-protocol's. This backend defines none; argument
;;;; checks are CHECK-TYPE (TYPE-ERROR with STORE-VALUE).
;;;;
;;;; Dependencies: transcribe-protocol 0.1.0 (TRANSCRIBER, DETECT-NOTES,
;;;; NOTE-EVENT, F0-CONTOUR, HINT-PITCH-RANGE, MIDI-TO-HZ, HZ-TO-MIDI) and
;;;; audio-protocol 0.2.0 (PCM-BUFFER, DB-TO-GAIN; the DSP package's
;;;; YIN-CMNDF, FRAME-COUNT, FRAME-START). Tests: rove 0.11.1 and
;;;; transcribe-protocol/testing (every take is synthesized; no fixtures).

(defpackage #:transcribe-backend-pyin
  (:use #:cl)
  (:local-nicknames (#:tp #:transcribe-protocol)
                    (#:audio #:audio-protocol)
                    (#:dsp #:audio-protocol/dsp))
  (:documentation "pYIN pitch tracking and HMM note segmentation as a
transcribe-protocol backend (slice T1).

(make-instance 'pyin-transcriber) is a TRANSCRIBER: pass it to
transcribe-protocol's TRANSCRIBE. Its DETECT-NOTES method runs PYIN-CONTOUR
(YIN's d' per frame, PITCH-CANDIDATES under THRESHOLD-PRIOR's beta prior, a
level gate, DECODE-PITCH's Viterbi over pitch bins and voicing), ONSET-FRAMES
(dips of a short-time level envelope, where a repeated pitch is sung again)
and CONTOUR-NOTES (SMOOTH-PITCHES's vibrato smoothing, then a note HMM with
a minimum note length). Tuning, drift and the key fit are the protocol's.")
  (:export
   ;; The transcriber (transcriber.lisp)
   #:pyin-transcriber
   #:pyin-transcriber-frame-size #:pyin-transcriber-hop #:pyin-transcriber-gate-db
   #:pyin-transcriber-min-duration #:pyin-transcriber-vibrato-window
   ;; Thresholds and candidates (candidates.lisp)
   #:+threshold-count+ #:+beta-prior+ #:+no-trough-probability+
   #:threshold-prior #:pitch-candidates
   ;; The pitch HMM and the contour (pitch.lisp)
   #:+bins-per-semitone+ #:+max-transition-rate+ #:+switch-probability+
   #:+probability-floor+
   #:pitch-bin-count #:pitch-max-jump #:decode-pitch #:pyin-contour
   ;; Onsets (onsets.lisp)
   #:+dip-db+ #:+dip-span+ #:onset-frames
   ;; The note HMM (notes.lisp)
   #:+note-steps-per-semitone+ #:+smoothing-band+ #:+note-sigma+ #:+prior-pitched+
   #:+note-on-probability+ #:+note-off-probability+ #:+note-change-probability+
   #:+min-change-steps+ #:+onset-change-probability+ #:+reattack-probability+
   #:+pitch-band+
   #:smooth-pitches #:contour-notes))

(in-package #:transcribe-backend-pyin)
