(defsystem "transcribe-backend-pyin"
  :version "0.1.0"
  :description "pYIN pitch tracking and HMM note segmentation as a transcribe-protocol backend: probabilistic YIN thresholds under a beta prior, a pitch and voicing HMM decoded by Viterbi, an energy-dip onset cue, and a note HMM with vibrato smoothing and a minimum note length, for one sung or hummed line at 16 kHz"
  :author "egao1980"
  :license "MIT"
  :depends-on ((:version "transcribe-protocol" "0.1.0")
               (:version "audio-protocol" "0.2.0"))
  :pathname "src"
  :serial t
  :components ((:file "package")
               (:file "transcriber")
               (:file "candidates")
               (:file "pitch")
               (:file "onsets")
               (:file "notes")
               (:file "detect"))
  :in-order-to ((test-op (test-op "transcribe-backend-pyin/tests"))))

(defsystem "transcribe-backend-pyin/tests"
  :depends-on ("transcribe-backend-pyin" "transcribe-protocol/testing"
               (:version "rove" "0.11.1"))
  :pathname "tests"
  :serial t
  :components ((:file "package")
               (:file "support")
               (:file "transcriber-test")
               (:file "candidates-test")
               (:file "pitch-test")
               (:file "onsets-test")
               (:file "notes-test")
               (:file "takes-test")
               (:file "cost-test"))
  :perform (test-op (o c)
             (unless (symbol-call :rove :run c)
               (error "tests failed for ~A" (component-name c)))))
