# Production Deployment

## Runtime

- Host: `wang@192.168.7.10` (SSH alias: `ubuntu-wg`)
- Checkout: `/home/wang/ubuntu_wang_code1/cas-offinder_with_primer3`
- Service: `cas-offinder-primer3.service`
- Application: Streamlit on port `8501`
- Health check: `curl --fail http://127.0.0.1:8501/_stcore/health`

## Release Procedure

1. Run the relevant local tests and push the verified commit to `origin/main`.
2. Connect with `ssh ubuntu-wg` and enter the checkout directory.
3. Inspect `git status --short`; do not reset or overwrite pre-existing server changes.
4. Run `git pull --ff-only origin main`.
5. Run `sudo systemctl restart cas-offinder-primer3.service`.
6. Confirm the service is active, port `8501` is listening, and the health check returns `ok`.

## Preserved Server Changes

As of 2026-07-22, the production checkout has one uncommitted source change:

- `src/otp/annotate.py` enables `pyranges`, which attempts to load the full GTF into memory. Do not deploy it as-is: a production probe consumed about `7 GiB` and caused a host OOM event.

Keep this change during a fast-forward release until a low-memory annotation design has been reviewed and committed. The CMake compatibility flag in `src/otp/cas_offinder.py` was separately verified and is now tracked in the repository.
