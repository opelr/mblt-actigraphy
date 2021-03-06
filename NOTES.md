## Watch Properties

* Sampling Freq.: 256Hz
    + Does that mean a 2-minute bin has a physical maxima of 30720? (256 $\times$ 60 $\times$ 2)
* Accelerometer Sensitivity: 0.01 g's
* Band Pass Filter: 3 - 11 Hz

## Scoring

* Activity $\Rightarrow A = 0.04E_{−2} + 0.2E_{−1} + 1E_{0} + 0.2E_{+1} + 0.04E_{+2}$ where $E_0$ is the epoch being scored.
* [Lotjonen](http://www.vivago.fr/SAS/PDF/analyses-cliniques/Sleep-wake-analysis.pdf)
    + Lotjonen_Wake is very conservative with sleep prediction --> let's try the more general regression instead of the "middle aged"

## Statistics and Analysis

* [Resampling -- jack-knife vs. bootstrap](https://en.wikipedia.org/wiki/Resampling_(statistics))
* 