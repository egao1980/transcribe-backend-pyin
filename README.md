# transcribe-backend-pyin

A [transcribe-protocol](https://github.com/egao1980/transcribe-protocol) backend for [cl-stack](https://github.com/egao1980/cl-stack): pYIN pitch tracking (Mauch and Dixon 2014: probabilistic YIN thresholds and a pitch HMM decoded by Viterbi) and HMM note segmentation in the manner of Tony, for one sung or hummed line at 16 kHz. Builds on [audio-protocol](https://github.com/egao1980/audio-protocol)'s DSP package.

**Status:** skeleton (0.0.1). The design is in the capability brief [cl-stack/docs/capabilities/transcribe.md](https://github.com/egao1980/cl-stack/blob/main/docs/capabilities/transcribe.md) (slice T1). The API and tests land through contract branches and team-solve.

```sh
cl-repo test --system transcribe-backend-pyin
```

## License

MIT, see [LICENSE](LICENSE).
