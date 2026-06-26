# Data Manifest

Copyright (c) 2026 Ziyi Yang. See `NOTICE.md`.

This release includes only compact reproduction data.

| File | Size | Purpose |
| --- | ---: | --- |
| `data/ANC_PriSec_pairs.mat` | 5.4 MB | Paired primary and secondary paths derived from PANDAR measurements. Used by both reproduction scripts. |
| `data/Wc_Sc_MAML_RWTH_Min_3.mat` | 6.1 KB | Pre-trained MAML co-initialization (`Wc`, `Sc`) for Training set B, used by `run_figure2_online_switch.m`. |
| `figures/Fig2_MSE_switch.jpg` | 890 KB | Final paper image for visual comparison. |
| `figures/Fig3_path_responses.jpg` | 265 KB | Final paper image for visual comparison. |

Excluded from the release:

- Large intermediate training caches from the working repository, including multi-GB `MetaTrainData*.mat` files.
- Raw PANDAR database folders. Use `data/ANC_PriSec_pairs.mat` for reproducible paper figures.
- Exploratory figures and temporary MATLAB workspaces.
