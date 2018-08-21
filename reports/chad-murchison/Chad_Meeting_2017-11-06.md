# Chad Murchison MBLT Meeting

## Actiware

* Ecosystem?
* Configuring and downloading watches
* Raw data (CSV) Data export

## Brief Overview of IRB #4085

* Morning bright light therapy's (MBLT) effect on sleep-wake disturbances in Veterans with & without chronic mild traumatic brain injury (mTBI). 
* Protocol
    - Consent at sleep clinic before their overnight PSG.
    - Actiwatch for one month.
    - Half of participants get light box, half get placebo.
    - Actigraphy -- one week of baseline recording before therapy starts, three weeks of therapy.
    - Surveys at consent, one month, two months, six months

## R Script

### ETL

* Read CSVs, `lapply()` into one big data.frame
* Date-time information --> sun-rise/-set, clock time, day, noon-noon day, etc.
* Sleep staging
    - Replicated actiware's sleep staging algorithm --> ((0.12 * [i-1]) + (0.5 * [i]) + (0.12 * [i+1])) > threshold
    - Activity > threshold
    - Lotjonen metrics (PSG-derived)
    - Stage smoothing (Cole; recursive)
* Light and activity adherence, smoothing
* Sleep Consoldiation
    - Require 8 uninterrupted epochs to be considered in a sleep/wake bout
* Off-wrist detection
    - If activity doesn't change for 3+ hours, nix the whole noon-noon day. 
    - Patients need at least 3 days of continuous, non-omitted recording to be included in analysis      

### Analysis

* Major sleep metrics -- TST, WASO, SE -- during sleep bouts
* Intradaily variability (IV) and interdaily stability (IS)
    - Non-parametric measures
    - IV = Quantifies how fragmented the rhythm is relative to the overall variance.
    - IS = The extent to which the profiles of individual days resemble each other.
* Relative Amplitude
    - Non-param. -- ((M10 - L5) / (M10 + L5))
* Light Adherence
* Activity Consoldiation

* Older functions/insights
    - Bedtime/waketime variability
    - Deprecated: Daytime activity ratio (7AM-7PM)
    - Deprecated: Nighttime sleep duration (7AM-7PM)
* Near future
    - Sleep latency
* Far Future
    - Use patient PSG sleep staging compared to their actigraphy during PSG to stage sleep for remainder of study
    - Sleep staging HMM

### Stats & Analysis

* Repeated Measures ANOVA
    - Choosing time points
* Time Series?

### Notes and Questions

- Relative Amplitude hours change each day. Since we're doing a repeated measures analysis, I presume that each patient should have the M10 and L5 hours chosed that averages to be highest/lowest across the entire study, not choosing M10 and L5 hours each day... yeah?
- Repeated Measures ANOVA with 2/3/4 points -- is more necessarily better?

## Talking Points

* What are intradaily variability (IV) and interdaily stability (IS)?
    - Non-parametric measures
    - IV
        - Quantifies how fragmented the rhythm is relative to the overall variance.
        - Calculation: Root mean square (RMS) of the difference between each hour's activity and the previous hour's, normalized by population variance.
    - IS
        - The extent to which the profiles of individual days resemble each other.
        - Calculation: Variance of the average daily profile divided by the total variance.

## Meeting Notes

* 