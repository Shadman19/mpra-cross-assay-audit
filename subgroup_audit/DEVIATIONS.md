# Protocol deviations

1. 2026-10-13: code/run_audit.py, modality_diagnostics() only. First run (job 7303282) completed all preregistered subgroup contrasts, then crashed in the descriptive per-modality accuracy export because score columns were read as text (TypeError in np.isfinite). Fix: pd.to_numeric(..., errors='coerce') on the score and effect columns. No preregistered test, threshold, family or correction changed. Lock digest updated from 6bb8a0fc6c670e65c9e07c17eed63e2271345d29be1baaa3de2faff81b871667 to f5d8856fcc8dce96a85f10580c17e74bc78124e558d70a0ea00634500426a62e. Failed run logs kept in logs_failed_1/.
Run 1 outputs (job 7303282) kept in results_run1_failed/; fresh results/ restored from the original zip for a clean versioned run.
