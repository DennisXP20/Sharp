"""Sternmaske: Wo darf das Sternmodell wirken?

Idee: Die Differenz zwischen "Both" (Stern- + Nicht-Stern-Modell) und "Non-Stellar Only" zeigt genau die
Stellen, an denen das Sternmodell etwas veraendert hat. Jede zusammenhaengende Stelle wird nur dann
uebernommen, wenn dort laut Katalog wirklich ein Stern sitzt.
"""
import numpy as np
from scipy import ndimage as ndi


def _lookup(labels, x, y, reach=2):
    """Label an den Positionen (mit kleiner Toleranz); 0 = keins."""
    h, w = labels.shape
    xi, yi = np.round(x).astype(int), np.round(y).astype(int)
    out = np.zeros(len(xi), dtype=labels.dtype)
    ok = (xi >= 0) & (xi < w) & (yi >= 0) & (yi < h)
    for dy in range(-reach, reach + 1):
        for dx in range(-reach, reach + 1):
            xx, yy = np.clip(xi + dx, 0, w - 1), np.clip(yi + dy, 0, h - 1)
            lab = labels[yy, xx]
            take = ok & (out == 0) & (lab > 0)
            out[take] = lab[take]
    return out


def _paint_disc(mask, x, y, r, value):
    h, w = mask.shape
    x0, x1 = max(0, int(x - r - 1)), min(w, int(x + r + 2))
    y0, y1 = max(0, int(y - r - 1)), min(h, int(y + r + 2))
    if x0 >= x1 or y0 >= y1:
        return
    yy, xx = np.ogrid[y0:y1, x0:x1]
    mask[y0:y1, x0:x1][(xx - x) ** 2 + (yy - y) ** 2 <= r * r] = value


def _paint_ellipse(mask, x, y, a, b, th, value):
    h, w = mask.shape
    r = max(a, b)
    x0, x1 = max(0, int(x - r - 1)), min(w, int(x + r + 2))
    y0, y1 = max(0, int(y - r - 1)), min(h, int(y + r + 2))
    if x0 >= x1 or y0 >= y1:
        return
    yy, xx = np.ogrid[y0:y1, x0:x1]
    u = (xx - x) * np.cos(th) + (yy - y) * np.sin(th)
    v = -(xx - x) * np.sin(th) + (yy - y) * np.cos(th)
    mask[y0:y1, x0:x1][(u / a) ** 2 + (v / b) ** 2 <= 1] = value


def star_mask(both, nonstellar, cat, fwhm_px, k_sigma=1.0, k_low=0.2, big_galaxy_stars="sicher", feather=1.5,
              verbose=True):
    """Maske 0..1 (H,W): 1 = Ergebnis des Sternmodells uebernehmen.

    both, nonstellar: (3,H,W) oder (H,W), linear.
    k_sigma: ab dieser Aenderung (in Hintergrund-Sigma) zaehlt eine Stelle als Eingriff des Sternmodells
    k_low:   bis zu dieser Aenderung wird die Stelle nach aussen verfolgt (damit Halos komplett mitkommen)
    big_galaxy_stars: 'sicher' = in grossen Galaxien nur Sterne mit sicherer Gaia-Astrometrie,
                      'alle'   = dort alle Punktquellen aus den Katalogen zulassen.
    """
    from scipy.spatial import cKDTree
    d = np.abs(np.asarray(both, dtype=np.float32) - np.asarray(nonstellar, dtype=np.float32))
    if d.ndim == 3:
        d = d.max(axis=0)
    ref = np.asarray(nonstellar[0] if np.ndim(nonstellar) == 3 else nonstellar)[::5, ::5]
    sigma = 1.4826 * np.median(np.abs(ref - np.median(ref)))
    d = ndi.gaussian_filter(d, 1.0)
    low = ndi.binary_opening(d > k_low * sigma, iterations=1)
    labels, n = ndi.label(low)
    peak = ndi.maximum(d, labels, index=np.arange(n + 1))
    del d
    real = peak > k_sigma * sigma          # Stellen mit echtem Eingriff (Hysterese)
    real[0] = False
    areas = ndi.sum_labels(low, labels, index=np.arange(n + 1))
    del low
    if verbose:
        print("  Sternmodell hat an %d Stellen eingegriffen" % real.sum())

    st, gx = cat["stars"], cat["galaxies"]
    allowed = np.ones(len(st["x"]), bool)
    if big_galaxy_stars == "sicher":
        allowed &= ~(st["in_big"] & ~st["sure"])
    s_lab = _lookup(labels, st["x"], st["y"])
    cand = allowed & real[s_lab]

    # vorlaeufige Radius-Helligkeits-Beziehung
    radius = lambda mag: np.clip(fwhm_px * (1.0 + 0.6 * (17.0 - mag)), 0.8 * fwhm_px, 60 * fwhm_px)
    cnt = np.bincount(s_lab[cand], minlength=n + 1)
    single = cand & (cnt[s_lab] == 1) & st["sure"]
    if single.sum() >= 30:
        r_eq = np.sqrt(areas[s_lab[single]] / np.pi)
        m = st["mag"][single]
        good = np.isfinite(m) & (r_eq > 0)
        slope, icpt = np.polyfit(m[good], np.log(r_eq[good]), 1)
        slope = min(slope, 0.0)
        radius = lambda mag: np.clip(np.exp(icpt + slope * mag) * 1.1, 0.8 * fwhm_px, 60 * fwhm_px)

    # Galaxien beanspruchen ihre Ellipse - ausser sie gehen im Halo eines sehr viel helleren Sterns unter
    g_use = (gx["mag"] < 24.0) & ~gx["big"]
    bright = np.flatnonzero(cand & st["sure"] & (st["mag"] < 15.0))
    if len(bright) and g_use.any():
        tree = cKDTree(np.c_[gx["x"], gx["y"]])
        for i in bright:
            for j in tree.query_ball_point([st["x"][i], st["y"][i]], float(radius(st["mag"][i]))):
                dist = np.hypot(gx["x"][j] - st["x"][i], gx["y"][j] - st["y"][i])
                near = dist < 0.6 * float(radius(st["mag"][i]))
                if gx["mag"][j] > st["mag"][i] + (3.0 if near else 7.0):
                    g_use[j] = False
    g_lab = _lookup(labels, gx["x"], gx["y"], reach=1)
    gmag = np.full(labels.shape, np.inf, dtype=np.float32)   # hellste beanspruchende Galaxie je Pixel
    order = np.flatnonzero(g_use & real[g_lab])
    order = order[np.argsort(-gx["mag"][order])]             # hellste zuletzt -> gewinnt
    for j in order:
        a = 2.5 * gx["a"][j] + 1.2 * fwhm_px
        b = 2.5 * gx["a"][j] * max(gx["q"][j], 0.1) + 1.2 * fwhm_px
        _paint_ellipse(gmag, gx["x"][j], gx["y"][j], a, b, gx["theta"][j], gx["mag"][j])

    xi = np.clip(np.round(st["x"]).astype(int), 0, labels.shape[1] - 1)
    yi = np.clip(np.round(st["y"]).astype(int), 0, labels.shape[0] - 1)
    ok_star = cand & (st["sure"] | (st["mag"] < gmag[yi, xi] - 0.5))
    keep = np.zeros(n + 1, bool)
    keep[s_lab[ok_star]] = True
    keep[0] = False

    claim = np.isfinite(gmag)
    del gmag
    for i in np.flatnonzero(ok_star):       # Stern, der auf/neben einer Galaxie steht, behaelt seinen Kern
        _paint_disc(claim, st["x"][i], st["y"][i], 0.5 * float(radius(st["mag"][i])), False)
    mask = keep[labels] & ~claim

    g_in = np.zeros(n + 1, bool)
    g_in[g_lab[g_use]] = True
    stats = dict(stellen=int(real.sum()), uebernommen=int((keep & real).sum()),
                 verworfen=int((real & ~keep).sum()),
                 verworfen_galaxie=int((real & ~keep & g_in).sum()),
                 verworfen_unbekannt=int((real & ~keep & ~g_in).sum()),
                 gemischt=int((keep & real & g_in).sum()))
    if verbose:
        print("  uebernommen: %(uebernommen)d (davon %(gemischt)d mit Galaxie im selben Fleck, Galaxie ausgespart), "
              "verworfen: %(verworfen)d (%(verworfen_galaxie)d Kataloggalaxien, %(verworfen_unbekannt)d ohne Katalogeintrag)" % stats)
    del labels, claim
    mask = ndi.binary_dilation(mask, iterations=2)
    mask = ndi.gaussian_filter(mask.astype(np.float32), feather)
    return np.clip(mask, 0, 1), stats


def blend(both, nonstellar, mask):
    """Sternmodell nur in der Maske: Ergebnis = Non-Stellar + Maske * (Both - Non-Stellar)."""
    return nonstellar + mask * (both - nonstellar)
