# Source Provenance

Copyright (c) 2026 Ziyi Yang. See `NOTICE.md`.

## Final Meta-Learning Training Path

The final reliable meta-learning training code path is:

- Original working script: `legacy/Main_MAML_Sec_Left3_training_original.m`
- Core classes used by that script:
  - `code/MAML_Nstep_forget.m`
  - `code/MAML_Nstep_forget_S2.m`
  - `code/FxLMS_phys.m`

This version jointly meta-learns:

- `Wc`: control-filter initialization.
- `Sc`: secondary-path initialization.

Earlier variants existed (`Main_MAML_Sec_Left.m`, `Main_MAML_Sec_Left2.m`, `tst4.m`, and single-filter MAML scripts), but the final ICASSP experiment artifacts are from the `Main_MAML_Sec_Left3` code path.

## Training-Set Mapping

The paper uses four three-task training sets:

| Paper set | Training rows | Artifact name |
| --- | --- | --- |
| A | `[32 37 43]` | `Wc_Sc_MAML_RWTH_Max_3.mat` |
| B | `[20 26 31]` | `Wc_Sc_MAML_RWTH_Min_3.mat` |
| C | `[4 25 46]` | `Wc_Sc_MAML_RWTH_C_3.mat` |
| D | `[9 12 21]` | `Wc_Sc_MAML_RWTH_D_3.mat` |

`Fig. 2` uses `Wc_Sc_MAML_RWTH_Min_3.mat`, matching the original `tst9.m` script.

## Figure Mapping

| Paper figure | Clean entry point | Original source |
| --- | --- | --- |
| Fig. 2 | `code/run_figure2_online_switch.m` | `legacy/tst9_three_phase_online_switch_original.m` |
| Fig. 3 | `code/run_figure3_path_responses.m` | `legacy/path_classify_original.m` |
