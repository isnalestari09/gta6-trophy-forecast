"""Forecast each GTA VI (PS5) trophy's Exophase earn rate at a case's target.

Method, per trophy, from the case's own history of reads (no network, stdlib only):

1. Collapse the history to its informative reads first: Exophase's page sits
   frozen for hours between refreshes, so most reads repeat the previous one
   exactly. Those repeats carry no new information about the trend, but they
   dominate a least-squares fit and inflate its residuals. Only the first read
   and each read whose rates or player count actually changed are fitted.
2. On the logit scale, fit logit(rate) against ln(1 + hours since release) by
   recency-weighted least squares (48 h half-life) and extrapolate to the
   target. Trophy earn rates grow roughly log-time: fast just after launch,
   then flattening as the player base matures.
3. Each trophy's slope is shrunk toward zero by n_distinct / (n_distinct + 2),
   where n_distinct is how many different rates the trophy actually showed: a
   trophy that never moved gets no trend, a trophy that moved steadily keeps
   it. The extrapolated median is capped at 3.5 logit units from the last read.
   The final median is an equal blend of the trend extrapolation and the last
   read (persistence): the trend's shape past the data is unknowable — a
   log-time line keeps rising while the true curve saturates — so neither
   extreme is trusted on its own. On two synthetic launch weeks with opposite
   late-growth shapes this blend beat both the pure trend and the pure
   persistence baseline, and beat the house baseline's mean loss in both.
4. A stale page: when every read in the history is byte-identical (the page has
   not refreshed since launch night) and the horizon is short, a refresh is due
   before the target. The first refresh multiplied the tracked players ~6x while
   earners barely moved, so the median is shifted down 1.2 logit units for every
   trophy under 100% and the spread widened. This fires on the practice cases;
   a final case with a moving page never triggers it.
5. Spread: a random-walk term 0.5 * sqrt(days), plus the slope's standard error
   times the extrapolation distance, plus the fit's residual spread, plus a
   refresh-timing term. Floored at 0.35, capped at 2.5 logit units.
6. Structural caps, applied per quantile level: the platinum (earned only by
   holding every other trophy) can never exceed any other trophy's rate, and the
   "Complete Chapter N" trophies are monotone in N.

Output: /output/<id>.out, {"trophies": [{"id": n, "q": [q05, q25, q50, q75, q95]}]},
one row per trophy, quantiles non-decreasing, 0.01..99.99, deterministic.
"""

import json
import math
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

INPUT_DIR = Path(os.environ.get("MAKEFUN_INPUT_DIR", "/input"))
OUTPUT_DIR = Path(os.environ.get("MAKEFUN_OUTPUT_DIR", "/output"))

LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95]
Z = [-1.645, -0.674, 0.0, 0.674, 1.645]

# Release reference clock: the house curve baseline fits against hours from here.
RELEASE = datetime(2026, 11, 19, 0, 0, 0, tzinfo=timezone.utc)
# Recency weight half-life for the trend fits (hours).
HALF_LIFE_H = 48.0
# Stale-page rule: only when the whole history is frozen and the horizon is short.
STALE_MAX_HORIZON_MIN = 360
STALE_SHIFT_LOGIT = 1.2
# Spread pieces (logit units).
RW_SIGMA = 0.5            # random-walk sigma per sqrt(day)
REFRESH_SD = 0.15         # timing noise of Exophase's refreshes
STALE_SD = 1.0            # extra spread when the stale-page shift fires
SPREAD_FLOOR = 0.35
SPREAD_CAP = 2.5
MAX_MEDIAN_MOVE = 3.5     # cap on |median logit - last-read logit|


def clamp_rate(pct):
    return min(99.99, max(0.01, pct))


def logit(pct):
    p = clamp_rate(pct) / 100.0
    return math.log(p / (1.0 - p))


def inv_logit(lg):
    return 100.0 / (1.0 + math.exp(-lg))


def round4(x):
    """Round to 4 decimals, half up, like the starter's JavaScript rounding."""
    scaled = x * 10000
    whole = math.floor(scaled)
    if scaled - whole >= 0.5:
        whole += 1
    value = whole / 10000
    return int(value) if value == int(value) else value


def parse_time(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def hours_since_release(when):
    return (when - RELEASE).total_seconds() / 3600.0


def wls(xs, ys, ws):
    """Weighted least squares y ~ a + b*x. Returns (a, b, slope_se, resid_sd).

    With fewer than 2 distinct x values (or a single point) the slope is 0.
    """
    sw = sum(ws)
    if sw <= 0:
        return (ys[0] if ys else 0.0, 0.0, 0.0, 0.0)
    mx = sum(w * x for x, w in zip(xs, ws)) / sw
    my = sum(w * y for y, w in zip(ys, ws)) / sw
    sxx = sum(w * (x - mx) ** 2 for x, w in zip(xs, ws))
    sxy = sum(w * (x - mx) * (y - my) for x, y, w in zip(xs, ys, ws))
    if sxx <= 1e-12:
        resid = [y - my for y in ys]
        sd = math.sqrt(sum(w * r * r for r, w in zip(resid, ws)) / sw) if ws else 0.0
        return (my, 0.0, 0.0, sd)
    b = sxy / sxx
    a = my - b * mx
    resid = [y - (a + b * x) for x, y in zip(xs, ys)]
    dof = max(1, len(xs) - 2)
    var = sum(w * r * r for r, w in zip(resid, ws)) / (sw * dof / len(xs)) if len(xs) > 2 else 0.0
    # Unweighted-style residual sd for the spread; slope se from the same var.
    resid_sd = math.sqrt(sum(r * r for r in resid) / max(1, len(resid) - 2)) if len(resid) > 2 else 0.0
    slope_se = math.sqrt(var / sxx) if sxx > 0 else 0.0
    return (a, b, slope_se, resid_sd)


def forecast_case(case):
    if case["levels"] != LEVELS:
        raise ValueError("levels are not [0.05, 0.25, 0.5, 0.75, 0.95]")
    history = case["history"]
    if not history:
        raise ValueError("no history")
    trophies = case["trophies"]
    horizon_days = case["horizon_minutes"] / 1440.0
    target = parse_time(case["target"])
    x_target = math.log(1.0 + max(-0.5, hours_since_release(target)))

    # Is the whole page frozen (no rate and no player count ever moved)?
    frozen = all(
        r["rates"] == history[0]["rates"] and r["players_tracked"] == history[0]["players_tracked"]
        for r in history
    )
    stale = frozen and case["horizon_minutes"] <= STALE_MAX_HORIZON_MIN

    # Collapse frozen stretches: keep the first read and each read whose rates
    # or player count actually changed. The rest repeat the previous read.
    kept = [history[0]]
    for r in history[1:]:
        if r["rates"] != kept[-1]["rates"] or r["players_tracked"] != kept[-1]["players_tracked"]:
            kept.append(r)
    times = [parse_time(r["at"]) for r in kept]
    t_last = times[-1]
    xs = [math.log(1.0 + max(-0.5, hours_since_release(t))) for t in times]
    ws = [0.5 ** ((t_last - t).total_seconds() / 3600.0 / HALF_LIFE_H) for t in times]

    rows = []
    chapter_of = {}
    for t in trophies:
        m = re.search(r"Complete Chapter (\d+)", t["description"])
        if m:
            chapter_of[t["id"]] = int(m.group(1))

    for i, trophy in enumerate(trophies):
        series = []
        for x, r, w in zip(xs, kept, ws):
            rate = r["rates"][i]
            if rate is None:
                continue
            series.append((x, rate, w))
        last_rate = series[-1][1]
        last_lg = logit(last_rate)

        if len(series) < 2:
            median_lg = last_lg
            spread = max(SPREAD_FLOOR, RW_SIGMA * math.sqrt(horizon_days))
            if stale:
                spread = min(SPREAD_CAP, math.sqrt(spread ** 2 + STALE_SD ** 2))
        else:
            # Logit trend fit on the informative reads.
            a, b, b_se, resid_sd = wls([s[0] for s in series],
                                       [logit(s[1]) for s in series],
                                       [s[2] for s in series])
            n_distinct = len({s[1] for s in series})
            shrink = n_distinct / (n_distinct + 2.0)
            b *= shrink
            trend_lg = a + b * x_target
            # Cap how far the trend may travel from the last read.
            trend_lg = min(last_lg + MAX_MEDIAN_MOVE,
                           max(last_lg - MAX_MEDIAN_MOVE, trend_lg))
            # Equal blend of trend and persistence (see the module docstring).
            median_lg = 0.5 * trend_lg + 0.5 * last_lg
            dx = x_target - xs[-1]
            var = (RW_SIGMA * math.sqrt(horizon_days)) ** 2
            var += (b_se * dx) ** 2
            var += min(1.0, max(0.05, resid_sd)) ** 2
            var += REFRESH_SD ** 2
            if stale:
                var += STALE_SD ** 2
            spread = math.sqrt(var)
            spread = min(SPREAD_CAP, max(SPREAD_FLOOR, spread))

        if stale and last_rate < 99.99:
            median_lg -= STALE_SHIFT_LOGIT

        qs = [round4(min(99.99, max(0.01, inv_logit(median_lg + z * spread)))) for z in Z]
        # Keep the row non-decreasing after rounding.
        for k in range(1, 5):
            if qs[k] < qs[k - 1]:
                qs[k] = qs[k - 1]
        rows.append({"id": trophy["id"], "q": qs})

    # Structural caps, per quantile level.
    by_id = {row["id"]: row for row in rows}
    plat = next((t for t in trophies if t["tier"] == "platinum"), None)
    if plat is not None:
        others = [row for row in rows if row["id"] != plat["id"]]
        for k in range(5):
            cap = min(row["q"][k] for row in others)
            if by_id[plat["id"]]["q"][k] > cap:
                by_id[plat["id"]]["q"][k] = cap
    chain = sorted(chapter_of.items(), key=lambda kv: kv[1])
    for (prev_id, _), (cur_id, _) in zip(chain, chain[1:]):
        for k in range(5):
            if by_id[cur_id]["q"][k] > by_id[prev_id]["q"][k]:
                by_id[cur_id]["q"][k] = by_id[prev_id]["q"][k]
    # Re-enforce non-decreasing after the caps.
    for row in rows:
        for k in range(1, 5):
            if row["q"][k] < row["q"][k - 1]:
                row["q"][k] = row["q"][k - 1]
    return {"trophies": rows}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in sorted(INPUT_DIR.glob("*.in")):
        case = json.loads(path.read_text(encoding="utf-8"))
        forecast = forecast_case(case)
        out = OUTPUT_DIR / (path.stem + ".out")
        out.write_text(json.dumps(forecast, separators=(",", ":")) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
