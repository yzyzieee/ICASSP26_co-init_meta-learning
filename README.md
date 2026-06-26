# ICASSP 2026 Co-Initialization Meta-Learning Code

This folder contains the cleaned MATLAB release package for:

> Co-Initialization of Control Filter and Secondary Path via Meta-Learning for Active Noise Control

The code reproduces the two main result figures from the published ICASSP 2026 paper:

- `Fig. 2`: online secondary-path modeling FxLMS under path switches at 60 s and 120 s.
- `Fig. 3`: measured PANDAR primary and secondary path responses with the two highlighted training sets.

## Copyright And Attribution

Copyright (c) 2026 Ziyi Yang. Contact: ziyi016@e.ntu.edu.sg.

This codebase modifies and extends MATLAB code originally released by Shi Dongyuan at [ShiDongyuan/Meta](https://github.com/ShiDongyuan/Meta). See `NOTICE.md` and `AUTHORS.md` for details.

## Folder Layout

```text
code/
  run_figure2_online_switch.m      Main script for paper Fig. 2
  run_figure3_path_responses.m     Main script for paper Fig. 3
  MAML_Nstep_forget.m              Control-filter MAML initializer
  MAML_Nstep_forget_S2.m           Secondary-path MAML initializer
  OnlineSPM_Zhang03.m              Reusable OSPM-FxLMS routine
  FxLMS.m, FxLMS_phys.m            FxLMS utilities
data/
  ANC_PriSec_pairs.mat             Paired PANDAR primary/secondary paths
  Wc_Sc_MAML_RWTH_Min_3.mat        Training set B co-initialization for Fig. 2
figures/
  Fig2_MSE_switch.jpg              Final figure used in the paper
  Fig3_path_responses.jpg          Final figure used in the paper
legacy/
  Original research scripts used to locate and verify the final figures
results/
  Generated outputs from the reproduction scripts
```

## Requirements

- MATLAB R2019b or newer is recommended.
- Signal Processing Toolbox is required for `fir1`, `freqz`, and `resample`.
- Communications Toolbox is not required by the cleaned scripts; measured AWGN is implemented locally.

## Quick Start

From MATLAB:

```matlab
cd('E:\NTU\AIANC\Meta-main\ICASSP2026_open_source\code')
run_figure3_path_responses
run_figure2_online_switch
```

Generated figures are written to `../results/`.

Verified locally with MATLAB R2024b:

- `run_figure3_path_responses` generated `results/Fig3_path_responses.png`.
- `run_figure2_online_switch` generated `results/Fig2_online_switch.png` in about 5.5 minutes on the test machine.

## Figure-to-Code Map

| Paper figure | Reproduction entry point | Original located file |
| --- | --- | --- |
| Fig. 2, online modeling FxLMS with auxiliary-noise power | `code/run_figure2_online_switch.m` | `legacy/tst9_three_phase_online_switch_original.m` |
| Fig. 3, primary/secondary path magnitude responses | `code/run_figure3_path_responses.m` | `legacy/path_classify_original.m` |

See `SOURCE_PROVENANCE.md` for the final meta-learning training script lineage.

## Data Notes

The release uses a compact paired-path file, `ANC_PriSec_pairs.mat`, rather than the full raw PANDAR database. The very large intermediate training caches from the working directory were intentionally excluded.

See `DATA_MANIFEST.md` for the included data files and exclusions.

The selected training rows are:

- Training set A: rows `[32 37 43]`
- Training set B: rows `[20 26 31]`
- Training set C: rows `[4 25 46]`
- Training set D: rows `[9 12 21]`
- Three-phase evaluation rows: `[10 25 38]`

`Fig. 2` uses the Training set B initializer (`Wc_Sc_MAML_RWTH_Min_3.mat`) because that is the artifact loaded by the final located plotting script.

## Citation

If this code is useful, please cite:

```bibtex
@inproceedings{yang2026coinitialization,
  title     = {Co-Initialization of Control Filter and Secondary Path via Meta-Learning for Active Noise Control},
  author    = {Yang, Ziyi and Rao, Li and Luo, Zhengding and Shi, Dongyuan and Huang, Qirui and Gan, Woon-Seng},
  booktitle = {ICASSP 2026 - 2026 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)},
  pages     = {15297--15301},
  year      = {2026},
  publisher = {IEEE},
  doi       = {10.1109/ICASSP55912.2026.11463219}
}
```

## License

No formal open-source license has been selected yet. Add a `LICENSE` file before public release.
