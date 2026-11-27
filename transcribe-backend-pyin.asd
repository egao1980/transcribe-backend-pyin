(defsystem "transcribe-backend-pyin"
  :version "0.0.1"
  :description "pYIN pitch tracking and HMM note segmentation as a transcribe-protocol backend"
  :author "egao1980"
  :license "MIT"
  :pathname "src"
  :serial t
  :components ((:file "package"))
  :in-order-to ((test-op (test-op "transcribe-backend-pyin/tests"))))

(defsystem "transcribe-backend-pyin/tests"
  :depends-on ("transcribe-backend-pyin" (:version "rove" "0.11.1"))
  :pathname "tests"
  :serial t
  :components ((:file "transcribe-backend-pyin-test"))
  :perform (test-op (o c)
             (unless (symbol-call :rove :run c)
               (error "tests failed for ~A" (component-name c)))))
