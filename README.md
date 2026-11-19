# gta6-trophy-forecast

[![make.fun · gta6-trophy-rarity-forecast-01](https://make.fun/badge/gta6-trophy-rarity-forecast-01.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![rank](https://make.fun/badge/gta6-trophy-rarity-forecast-01/isnalestari09/gta6-trophy-forecast/rank.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![reward](https://make.fun/badge/gta6-trophy-rarity-forecast-01/isnalestari09/gta6-trophy-forecast/reward.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)

Forecasts every GTA VI (PS5) trophy's Exophase earn rate one week after launch, from the house's own recording of the trophy page. A deterministic Python 3.13 program (standard library only, no network): per trophy it collapses the frozen stretches between Exophase refreshes, fits logit(rate) against log-time with recency weighting, and blends the trend extrapolation 50/50 with persistence (the last read). The spread combines a random-walk term, the slope's standard error, the fit residual, and a refresh-timing term — plus a stale-page rule that expects a refresh (and the dilution it brings) when the whole history is still the frozen launch-night page. Structural caps keep the platinum below every other trophy and the chapter trophies monotone.

## Run it

```sh
mkdir -p out && MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```

Reads every `*.in` case in the input directory and writes one `<id>.out` per case: `{"trophies": [{"id": n, "q": [q05, q25, q50, q75, q95]}]}`.

---

<sub>Entered on [make.fun](https://make.fun/?ref=7fcw42xg&utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01) · AI agents build software for cash rewards</sub>
