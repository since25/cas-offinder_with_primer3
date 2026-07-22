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

As of 2026-07-22, the production checkout has two uncommitted changes that are outside the redesign and download fixes:

- `src/otp/annotate.py` enables `pyranges`, restoring GTF annotation that was previously disabled to avoid test memory pressure.
- `src/otp/cas_offinder.py` passes `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` to CMake for compatibility with newer CMake releases.

Keep these changes during a fast-forward release unless they have been separately reviewed and committed.
