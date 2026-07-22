## 2026-07-22 - Task: Handle blank rows in primer redesign input tables

### What was done

- Made primer redesign ignore completely blank spreadsheet rows.
- Added a clear validation error for partially filled rows, including the data row, Excel row, and missing required fields.
- Documented the input-row behavior for existing off-target tables.

### Testing

- `PYTHONPATH="$PWD/src:$PWD/.venv/lib/python3.14/site-packages" python3 -m pytest -q`: 38 passed.
- Submitted workbook parse check: 30 rows read; 21 rows before the first blank separator were validated.
- `python3 -m otp.redesign --genome hg38 --input human_抽检2次rep反应重合位点.xlsx`: processed 21 rows; primers found for all 21 rows; generated `results.csv` and `results.xlsx` in a system temporary directory.

### Notes

Changed files:
- `src/otp/redesign.py`: skips blank input rows and validates required values before mapping.
- `tests/test_redesign.py`: covers blank-row and partial-row validation behavior.
- `README.md`: documents blank-row handling.
- `docs/redesign-input.md`: records the required input-row rules.
- `progress.md`: records this change and its validation evidence.

Rollback:
- Revert the changes from this task.

## 2026-07-22 - Task: Download complete Excel workbook from the web UI

### What was done

- Removed the CSV download button from the Streamlit results page.
- Kept `results.csv` as an internal flat result file for the page preview, cache, and command-line compatibility.
- Made the sole web download explicitly the Excel workbook containing all result sheets.

### Testing

- `python3 -m py_compile app/streamlit_app.py`: passed.
- Streamlit source assertion: one Excel workbook download button and no CSV download button.
- `PYTHONPATH="$PWD/src:$PWD/.venv/lib/python3.14/site-packages" python3 -m pytest -q`: 38 passed.

### Notes

Changed files:
- `app/streamlit_app.py`: exposes only the complete Excel workbook download.
- `README.md`: documents the web download behavior.
- `progress.md`: records this change and its validation evidence.

Rollback:
- Restore the CSV `st.download_button` in `app/streamlit_app.py`.

## 2026-07-22 - Task: Record production deployment procedure

### What was done

- Recorded the production host, checkout, service unit, health check, and release procedure.
- Recorded the two pre-existing production-only source changes that must be preserved during deployment.
- Added the deployment result to the persistent Codex project memory.

### Testing

- Verified production commit `b550c8d`, service `cas-offinder-primer3.service`, port `8501`, and the Streamlit health response `ok` on `wang@192.168.7.10`.
- Inspected the remote diffs for `src/otp/annotate.py` and `src/otp/cas_offinder.py` without modifying them.

### Notes

Changed files:
- `docs/production-deployment.md`: documents the verified production deployment process.
- `progress.md`: records the deployment evidence and preservation constraint.

Rollback:
- Remove the deployment documentation if the runtime layout is retired.

## 2026-07-22 - Task: Track the CMake compatibility fix

### What was done

- Added the verified CMake policy compatibility flag to the tracked Cas-OFFinder fallback build path.
- Added a regression test for the exact CMake invocation.
- Recorded that the separate pyranges change is unsafe to deploy because full-GTF loading caused a production memory exhaustion event.

### Testing

- `PYTHONPATH="$PWD/src:$PWD/.venv/lib/python3.14/site-packages" python3 -m pytest tests/test_cas_offinder.py -q`: 1 passed.
- `PYTHONPATH="$PWD/src:$PWD/.venv/lib/python3.14/site-packages" python3 -m pytest -q`: 39 passed.
- Temporary CMake build with `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`: completed successfully.

### Notes

Changed files:
- `src/otp/cas_offinder.py`: uses the CMake policy compatibility flag during fallback builds.
- `tests/test_cas_offinder.py`: covers the fallback CMake invocation.
- `docs/production-deployment.md`: records the reviewed status of the two remote changes.
- `progress.md`: records this merge decision and verification evidence.

Rollback:
- Remove `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` from the fallback CMake command.
