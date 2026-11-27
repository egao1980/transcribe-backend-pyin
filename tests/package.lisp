;;;; The T1 suite. The tests are the specification; the file headers in src/
;;;; state the rules they check. Every take is synthesized here with
;;;; transcribe-protocol/testing: no fixture files.

(defpackage #:transcribe-backend-pyin/tests
  (:use #:cl #:rove #:transcribe-backend-pyin
        #:transcribe-protocol #:transcribe-protocol/testing)
  (:local-nicknames (#:audio #:audio-protocol)
                    (#:dsp #:audio-protocol/dsp)
                    (#:mt #:music-theory)))

(in-package #:transcribe-backend-pyin/tests)
