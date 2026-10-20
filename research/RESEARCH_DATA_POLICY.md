# Research data integrity

Author instruction: interrupted runs must be marked and excluded from final completed-run statistics and convergence comparisons. Preserve partial checkpoints separately. Never interpret an absent final file as successful completion.

Every new experimental run records its ID, requested parameters and running status before computation. Final status records achieved time and stop reason. A user stop, wall limit, numerical failure and completed negative physical result are distinct outcomes. Numerical failures and negative findings remain visible in the run inventory; they must not be discarded to improve the apparent results.

Only runs with completed status and matching achieved/requested duration enter completed-run comparisons. An interrupted run's continuation or rerun gets a separate ID and explicit provenance. Report the matching rules and run IDs used for each comparison. Main console physical equations remain unchanged; Authorship-dependent dynamics require separately identified experimental work.
