# transcribe-backend-pyin

A [transcribe-protocol](https://github.com/egao1980/transcribe-protocol) backend for [cl-stack](https://github.com/egao1980/cl-stack): one sung or hummed line at 16 kHz in, the notes as sung and the f0 contour out. pYIN tracks the pitch, and a note HMM in the manner of Tony turns the pitch track into notes. The protocol does the rest: the take's tuning, its drift per phrase and the key fit.

- **pYIN** (Mauch and Dixon 2014): YIN's cumulative mean normalized difference d′ per frame (audio-protocol's DSP package), 100 thresholds under a Beta(2, 18) prior, the period candidates with their probabilities, and a hidden Markov model over 10-cent pitch bins and voicing, decoded by Viterbi. Frames of 1024 samples every 128 (64 ms every 8 ms), a level gate 40 dB under the loudest frame.
- **Onsets:** dips of a short-time level envelope (8 ms windows every 2 ms, 6 dB under both sides within 50 ms), where a singer re-articulates a note. A repeated pitch splits there.
- **Notes:** vibrato smoothing (a 200 ms median, then a mean over the frames on the same side of any step, cut at onsets), then a note HMM: attack and stable states every 1/3 semitone and one silent state, a minimum note length of 1/10 s, and a pitch per note that is a robust mean of its frames. Vibrato of a semitone's depth stays one note.

```lisp
(asdf:load-system "transcribe-backend-pyin")

(let ((take (audio-protocol:read-wav #p"take.wav")))
  (transcribe-protocol:transcribe (make-instance 'transcribe-backend-pyin:pyin-transcriber) take
                                  :hints '(:key "F" :pitch-range (48 79))))
;; → a TRANSCRIPTION: NOTES (retuned), RAW-NOTES (as sung), CONTOUR (125
;;   frames a second, voicing per frame), TUNING-OFFSET, DRIFT, KEY-FIT, ...
```

**Status:** slice T1 (0.1.0) is a **contract**: the functions below are stubs, and a team-solve run implements them. The design is in the capability brief [cl-stack/docs/capabilities/transcribe.md](https://github.com/egao1980/cl-stack/blob/main/docs/capabilities/transcribe.md). Each file header in `src/` states its rules, and the Rove suite in `tests/` is the spec.

| File | What | Status |
|------|------|--------|
| `src/transcriber.lisp` | `pyin-transcriber`: initargs (`:frame-size`, `:hop`, `:gate-db`, `:min-duration`, `:vibrato-window`), their checks, `transcriber-info` | T1, given |
| `src/candidates.lisp` | `+threshold-count+`, `+beta-prior+`, `+no-trough-probability+` | T1, given |
| | `threshold-prior` (Beta(2, 18) over 100 thresholds, exact), `pitch-candidates` (one frame's troughs, periods and probabilities) | T1, contract |
| `src/pitch.lisp` | `+bins-per-semitone+`, `+max-transition-rate+`, `+switch-probability+`, `+probability-floor+` | T1, given |
| | `pitch-bin-count`, `pitch-max-jump`, `decode-pitch` (the pitch and voicing HMM, Viterbi), `pyin-contour` | T1, contract |
| `src/onsets.lisp` | `+dip-db+`, `+dip-span+` | T1, given |
| | `onset-frames` (re-articulation dips) | T1, contract |
| `src/notes.lisp` | the note HMM's constants | T1, given |
| | `contour-notes` (vibrato smoothing, the note HMM, note events) | T1, contract |
| `src/detect.lisp` | the `detect-notes` method | T1, contract |

```sh
CL_REPO_DYNAMIC_SPACE_MB=2048 cl-repo test --system transcribe-backend-pyin
```

## License

MIT, see [LICENSE](LICENSE).
