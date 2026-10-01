# CLAUDE.md

## Coding Agent

- In comment replies, avoid `#<numeral>` style references such as "#1" unless
  you mean an issue or pull request, since GitHub turns them into links. Write
  "No. 1" or "number 1" instead.
- Include plots directly in comment replies via
  `![image name](https://github.com/sgbaird-5DOF/interp/blob/<short-hash>/<path>?raw=true)`.
  Use the 7-character commit hash, never a branch name, so the link keeps
  pointing at the exact version you produced.
- When you mention files in a comment reply, hyperlink them the same way
  (`https://github.com/sgbaird-5DOF/interp/blob/<short-hash>/<path>`).
- Never echo, grep, or print environment secrets. This repository is public,
  so anything in a log, commit, or comment is visible to everyone.

## Repository

- MATLAB code for five degree-of-freedom (5DOF) grain boundary interpolation:
  grain boundary octonions (GBOs), the Voronoi fundamental zone (VFZ), GB
  distances, and GPR/barycentric/NN/IDW interpolation. The high-level entry
  point is `code/interp5DOF.m`; `code/Contents.m` describes the other files.
- The default branch is `master`, not `main`.
- Requires MATLAB R2019b or newer (`arguments ... end` validation is used
  throughout). `*_r2018a.m` files are backports for older releases; keep them
  in sync when you change their modern counterparts.
- Tests are the `*_test.m` files under `code/`, run by
  `.github/workflows/test-matlab.yml` (`addpath(genpath('.'))`, then
  `matlab-actions/run-tests`). The quickest end-to-end check is
  `interp5DOF_test`.
- MATLAB is not installed on the Claude runner. If a change needs MATLAB to
  verify, push it and read the `test-matlab.yml` results on the PR (you have
  `actions: read`) rather than claiming it works.
- The codebase is mid-refactor (see the README note dated 2026-02-24),
  including a helper function in one of the GB conversion scripts. Prefer
  small, well-scoped changes and say which functions you touched.

### The GitHub token dies at minute 60 (push early, re-mint to continue)

The GitHub App token a session starts with (`GITHUB_TOKEN`/`GH_TOKEN` in the
session environment, and embedded in the origin remote URL) expires exactly 60
minutes after the Run Claude Code step starts. The session keeps running, but
`git push` fails, `gh` fails, and the MCP tool that updates the progress
comment fails, all with 401, so from the outside the session goes silent while
finished work stops landing.

- Push and update the tracking comment early and often. Treat minute 50 as the
  deadline for anything that must reach GitHub, in case recovery fails.
- At the first 401 from a push or `gh` call (or proactively around minute 55),
  run `python scripts/refresh_github_app_token.py`. It re-runs the action's
  own OIDC exchange, saves a fresh one-hour token to `/tmp/.ghtok` (mode
  0600), and re-points the origin remote at it, so plain `git push` works
  again. Validated live on a 125 minute session that re-minted hourly.
- `gh` keeps reading the dead token from the environment, so prefix each call:
  `GH_TOKEN=$(cat /tmp/.ghtok) gh ...`. Separately, `DEFAULT_WORKFLOW_TOKEN`
  is a distinct token that lasts the whole job and works for reads at any age.
- The MCP comment tool cannot be re-keyed mid-session. After a re-mint, update
  the tracking comment over REST instead: write the body to a file and run
  `GH_TOKEN=$(cat /tmp/.ghtok) gh api -X PATCH
  repos/$GITHUB_REPOSITORY/issues/comments/<comment-id> -F body=@that-file`.
- A re-minted token also lives one hour, so re-run the script each hour it is
  needed. Never echo, log, or commit a token value; the script prints only
  statuses and lengths.
