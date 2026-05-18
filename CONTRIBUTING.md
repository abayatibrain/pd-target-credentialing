# Contributing

This is a portfolio repository maintained by Armin Bayati. External contributions
are welcome but should expect that the maintainer will weigh changes against the
overall portfolio direction (see the [Cowork brief](https://github.com/arminbayati/cowork-brief) — private).

## Ground rules
1. Every non-trivial decision lives in an ADR under `docs/decisions/`.
   Before writing code, write the ADR.
2. Code style is enforced by pre-commit (`pre-commit install`).
3. No live API calls in tests — mock with `responses` or `respx`.
4. Datasets are downloaded by `scripts/download_data.sh`; nothing larger
   than 1 MB lives in the repo.
5. No claim in the README without a backing notebook cell or citation.

## Workflow
- Branches: `feat/<slug>`, `fix/<slug>`, `docs/<slug>`. Squash-merge to `main`.
- Conventional Commits in commit messages.
- Tag `v0.1.0` once the demo notebook runs end-to-end; `v1.0.0` once the
  per-repo acceptance criteria are met (see the relevant section of the brief).
