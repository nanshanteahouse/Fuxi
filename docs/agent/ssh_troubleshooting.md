# Pre-push Hooks & SSH Keep-alive

Fuxi ships `pre-push` hooks via pre-commit (see `.pre-commit-config.yaml`):
pyright (type check) + smoke test (`pytest tests/test_smoke.py` ≈ 20–45s).
Total pre-push wall: typically 50–80s on WSL2 (anndata import + RAPIDS
preload + subprocess startup).

GitHub's SSH server closes idle connections after ~60s. Without keep-alive,
the connection drops during the pre-push window and `git push` fails with
`Connection to github.com closed by remote host` — requiring `--no-verify`
to push, which bypasses the smoke gate.

## Fix: One-time setup per clone

```bash
git config --local core.sshCommand "ssh -o ServerAliveInterval=15 -o ServerAliveCountMax=6"
```

This sends a keep-alive packet every 15s (up to 6 retries = 90s tolerance),
keeping the SSH tunnel alive while pre-push hooks run. The smoke test
parallelization in `tests/test_smoke.py` (concurrent.futures for the 3
modality `--list` checks) keeps hook wall under 60s; combined with keep-
alive, push works reliably without `--no-verify`.

## If you still see SSH disconnects

- Verify `git config --local --get core.sshCommand` returns the keep-alive string
- Or run `GIT_SSH_COMMAND="ssh -o ServerAliveInterval=15" git push` for one push
- As a last resort, `git push --no-verify` (skips smoke gate)
