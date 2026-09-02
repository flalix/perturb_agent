"""
Build a survival frame from the cBioPortal clinical/demographic table.

The table merges TCGA and CPTAC3, and the two write survival differently:
TCGA-derived rows carry os_months / os_status ("1:DECEASED"), while CPTAC rows
often carry only follow_up_days / vital_status. Neither column alone covers the
cohort -- reading os_months by itself silently drops CPTAC, and reading
follow_up_days by itself censors everyone who died, since days-to-last-followup
is NaN for deceased patients by construction.
"""

import numpy as np
import pandas as pd


DEAD = {"1:deceased", "deceased", "dead", "1"}
ALIVE = {"0:living", "living", "alive", "0"}


def _status_to_event(v):
    """1 = death observed, 0 = censored, NaN = unknown. Never guesses."""
    if pd.isna(v):
        return np.nan
    s = str(v).strip().lower()
    if s in DEAD:
        return 1.0
    if s in ALIVE:
        return 0.0
    return np.nan


def survival_frame(df_demog, id_col="barcode_case", verbose=True):
    """Return one row per patient with months, event, and covariates.

    Prefers os_months/os_status; falls back to follow_up_days/vital_status
    where the primary pair is missing. Reports which source each row used, so
    a cohort covered entirely by the fallback is visible rather than assumed.
    """
    d = df_demog.copy()
    if id_col in d.columns:
        d = d.set_index(id_col)

    dup = d.index.duplicated()
    if dup.any():
        if verbose:
            print(f"{dup.sum()} duplicate {id_col}; keeping first")
        d = d[~dup]

    m_os = pd.to_numeric(d.get("os_months"), errors="coerce")
    e_os = d.get("os_status", pd.Series(index=d.index, dtype=object)).map(_status_to_event)

    days = pd.to_numeric(d.get("follow_up_days"), errors="coerce")
    m_fu = days / 30.4375
    e_fu = d.get("vital_status", pd.Series(index=d.index, dtype=object)).map(_status_to_event)

    primary = m_os.notna() & e_os.notna()
    months = m_os.where(primary, m_fu)
    event = e_os.where(primary, e_fu)
    source = pd.Series(np.where(primary, "os_months",
                       np.where(months.notna() & event.notna(), "follow_up_days", "none")),
                       index=d.index)

    out = pd.DataFrame({"months": months, "event": event, "source": source})

    # Non-positive follow-up is not a survival time: 0 contributes no risk
    # period and negative values are data entry errors.
    bad = out.months.notna() & (out.months <= 0)
    if bad.any():
        if verbose:
            print(f"dropping {int(bad.sum())} row(s) with months <= 0")
        out.loc[bad, ["months", "event"]] = np.nan

    for col, new in [("age", "age"), ("gender", "sex"),
                     ("ajcc_pathologic_tumor_stage", "stage"),
                     ("cbioportal_study_id", "study"),
                     ("history_neoadjuvant_trtyn", "neoadjuvant"),
                     ("radiation_therapy", "radiation"),
                     ("person_neoplasm_cancer_status", "tumor_status")]:
        if col in d.columns:
            out[new] = d[col]

    if "age" in out:
        out["age"] = pd.to_numeric(out.age, errors="coerce")
    if "sex" in out:
        out["male"] = (out.sex.astype(str).str.lower().str.startswith("m")).astype(float)
        out.loc[out.sex.isna(), "male"] = np.nan
    if "stage" in out:
        # AJCC I/II/III/IV -> ordinal; "Stage IIB" and "STAGE IIB" both map to 2
        rom = {"IV": 4, "III": 3, "II": 2, "I": 1}
        def _stage(v):
            if pd.isna(v):
                return np.nan
            s = str(v).upper().replace("STAGE", "").strip()
            for r, n in rom.items():          # IV before I, longest first
                if s.startswith(r):
                    return float(n)
            return np.nan
        out["stage_num"] = out.stage.map(_stage)

    if verbose:
        n = len(out)
        ok = out.months.notna() & out.event.notna()
        print(f"{n} patients; {ok.sum()} with usable survival "
              f"({int(out.event[ok].sum())} deaths, "
              f"{int((1 - out.event[ok]).sum())} censored)")
        print(out.source.value_counts().to_dict())
        if "study" in out:
            print(out[ok].groupby("study").size().to_dict())
    return out


def normalise_ids(idx):
    """Strip the T-/N- prefix and any TCGA sample suffix, upper-case."""
    return (pd.Index(idx).astype(str)
            .str.replace(r"^[TN]-", "", regex=True)
            .str.replace(r"-\d{2}[A-Z]?$", "", regex=True)
            .str.upper())


def join_groups(surv, disc, cols=("basal_minus_classical", "myCAF_minus_iCAF")):
    """Attach the 2x2 state assignment, reporting what fails to match."""
    g = disc[list(cols)].copy()
    g.index = normalise_ids(g.index)
    s = surv.copy()
    s.index = normalise_ids(s.index)

    inter = s.index.intersection(g.index)
    print(f"{len(inter)} of {len(g)} discretised samples matched to clinical")
    miss = g.index.difference(s.index)
    if len(miss):
        print(f"  unmatched: {list(miss)[:6]}")

    d = s.join(g, how="inner")
    d = d.dropna(subset=["months", "event"] + list(cols))
    print(d.groupby(list(cols))["event"].agg(n="size", deaths="sum"))
    return d
