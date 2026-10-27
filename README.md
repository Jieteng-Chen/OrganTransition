# OrganTransition

Takes nine raw `[0, 1]` organ FI values for one individual, computes the
overall and per-organ transition scores with the reference networks, and returns
their percentiles. 

## Installation

```r
# From GitHub
remotes::install_github("Jieteng-Chen/OrganTransition")
```

## Quick start

```r
library(OrganTransition)
ref <- load_reference()

fi <- data.frame(
  cardiovascular   = 0.05, hepatic         = 0.02, urinary       = 0.04,
  endocrine        = 0.03, immune          = 0.06, musculoskeletal = 0.06,
  neural           = 0.04, gastrointestinal = 0.03, respiratory   = 0.05
)

# Overall + all nine organ percentiles
transition_percentile(fi, ref, age = 68, sex = "woman")

# Overall percentile only
transition_percentile(fi, ref, age = 68, sex = "woman", organ = NULL)

# Selected organs
transition_percentile(fi, ref, age = 68, sex = "woman",
                      organ = c("hepatic", "cardiovascular"))
```

## Input

- Nine organ FI columns named exactly:
  ```
  cardiovascular, hepatic, urinary, endocrine, immune,
  musculoskeletal, neural, gastrointestinal, respiratory
  ```
- `age` and `sex` are optional. They are used only for age-sex percentile
  matching. 


## Output

A data.frame with:

- `percentile_all`: percentile in the full reference cohort;
- `percentile_age_sex`: percentile in the same sex and nearby age band;
- `n_all`, `n_age_sex`, `matched_by`, `age_band`, `age`, `sex`;
- `percentile_all_<organ>` and `percentile_age_sex_<organ>` for each requested
  organ (`organ = "all"` by default).

