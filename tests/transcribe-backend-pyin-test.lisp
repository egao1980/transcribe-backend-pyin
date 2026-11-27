(defpackage #:transcribe-backend-pyin/tests
  (:use #:cl #:rove))

(in-package #:transcribe-backend-pyin/tests)

(deftest skeleton-loads
  (ok (find-package '#:transcribe-backend-pyin)))
