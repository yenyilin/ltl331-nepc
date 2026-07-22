# Environment capture

**Capture exact versions from the machine that produced the results** — do not
hand-transcribe (Reviewer 2 minor #3 asks for versions).

## File layout

| File | Role | Commit |
|---|---|---|
| `requirements.in` | top-level wishlist (canonical analysis env, server / paper run) | yes |
| `requirements.txt` | locked versions for the canonical env | yes |
| `requirements-local.in` | wishlist for the Intel-Mac local exploration env (CPU-only, slimmer) | yes |
| `requirements-local.txt` | locked versions for the local env | yes |
| `session.txt` | runtime version capture from the machine that produced the results | yes |

## 1. Generate the locked files from the wishlists

```
uv pip compile env/requirements.in       -o env/requirements.txt
uv pip compile env/requirements-local.in -o env/requirements-local.txt
```
Re-run whenever an `.in` changes; commit the resulting `.txt` alongside.

## 2. Install

Canonical (server / paper run):
```
uv venv .venv --python 3.11
source .venv/bin/activate
uv pip install -r env/requirements.txt
```

Local exploration on Intel Mac (CPU only):
```
uv venv .venv-local --python 3.11
source .venv-local/bin/activate
uv pip install -r env/requirements-local.txt
```

## 3. Capture actual versions from the analysis env that produced the results

```
uv pip freeze > env/requirements.txt
python -c "import sys, scanpy; print(sys.version); scanpy.logging.print_versions()" \
  > env/session.txt
```

## Notes on the two Python envs

- **Canonical** (`requirements.in`): full set including moscot[scrnaseq] (uses
  JAX), celltypist, decoupler, scib. This is what the server runs and the
  paper's results are produced from. Pertpy / scvi-tools are intentionally
  excluded — Stack C (CPU-friendly) was chosen during the revision; see
  `review/patient_level_stats_plan.md`.
- **Local-Intel-Mac** (`requirements-local.in`): slimmer; commented-out
  moscot/celltypist/decoupler/scib so they can be added when needed.
  PyTorch-pulling deps (pertpy, scvi-tools) are deliberately omitted because
  PyTorch dropped Intel-Mac x86_64 wheels after torch 2.2.x.

Tool versions for the upstream R-based steps (Cell Ranger / Seurat processing,
inferCNV, GSEA) are recorded in `../config/params.yaml` (`versions:`); that R
code is outside this Python release.
