"""Beugungsspikes: im Bild finden und vor dem Sternmodell schuetzen.

Cosmic Clarity behandelt Spikes kachelweise unterschiedlich (mal duenn, mal breit), und entlang der Spikes stehen in
den Katalogen Scheinobjekte. Beides fuehrt zu zerhackten Spikes. Deshalb wird auf jedem Spike durchgehend die
Nicht-Stern-Version genommen; nur der Sternkern kommt vom Sternmodell.
"""
import numpy as np
from scipy import ndimage as ndi


def _sampler(img):
    h, w = img.shape

    def sample(xs, ys):
        ok = (xs >= 0) & (xs <= w - 1) & (ys >= 0) & (ys <= h - 1)
        out = np.full(xs.shape, np.nan, dtype=np.float32)
        if ok.any():
            out[ok] = ndi.map_coordinates(img, [ys[ok], xs[ok]], order=1, mode="nearest")
        return out
    return sample


def _smooth(v, n):
    n = max(1, int(n))
    good = np.isfinite(v)
    a = ndi.uniform_filter1d(np.where(good, v, 0).astype(np.float64), n, mode="nearest")
    b = ndi.uniform_filter1d(good.astype(np.float64), n, mode="nearest")
    return np.where(b > 0.3, a / np.maximum(b, 1e-9), np.nan)


def find_angle(lum, stars, fwhm, bg, sigma, verbose=True):
    """Spike-Winkel (Grad, 0..90) aus den hellsten Sternen; None, wenn es keine Spikes gibt.

    Je Stern wird ein Polar-Gitter abgetastet; pro Radius wird der Ring-Median abgezogen und durch die
    Ring-Streuung geteilt. So zaehlen nur schmale radiale Linien, nicht Halo oder Nachbarsterne.
    """
    h, w = lum.shape
    sample = _sampler(lum)
    order = np.argsort(stars["mag"])
    inb = stars.get("in_big", np.zeros(len(stars["mag"]), bool))
    cand = [i for i in order if not inb[i] and 300 < stars["x"][i] < w - 300 and 300 < stars["y"][i] < h - 300][:25]
    if len(cand) < 3:
        return None
    step = 0.25
    ang = np.radians(np.arange(0, 360, step))
    r = np.arange(3 * fwhm, min(40 * fwhm, 300), 1.0)
    total = np.zeros(int(180 / step))
    used = 0
    for i in cand:
        cx, cy = stars["x"][i], stars["y"][i]
        A, R = np.meshgrid(ang, r)                    # (radien, winkel)
        v = sample(cx + R * np.cos(A), cy + R * np.sin(A))
        med = np.nanmedian(v, axis=1, keepdims=True)
        mad = 1.4826 * np.nanmedian(np.abs(v - med), axis=1, keepdims=True) + 1e-12
        z = np.clip((v - med) / mad, -5, 50)
        prof = np.nanmedian(z, axis=0)                # robust gegen einzelne Nachbarsterne
        prof = prof[: len(total)] + prof[len(total):]  # a und a+180 zusammen
        if np.isfinite(prof).all():
            total += prof
            used += 1
    if used < 3:
        return None
    q = total + np.roll(total, -int(90 / step))       # a und a+90 zusammen (Kreuz)
    q = q[: int(90 / step)]
    k = int(np.argmax(q))
    contrast = (q[k] - np.median(q)) / (1.4826 * np.median(np.abs(q - np.median(q))) + 1e-12)
    if verbose:
        print("  Spike-Suche: Winkel %.2f Grad, Kontrast %.1f (%d Sterne)" % (k * step, contrast, used))
    return float(k * step) if contrast > 8 else None


def _arm_profile(sample, cx, cy, phi, r, off):
    c, s = np.cos(phi), np.sin(phi)
    ridge = sample(cx + r * c, cy + r * s)
    s1 = sample(cx + r * c - off * s, cy + r * s + off * c)
    s2 = sample(cx + r * c + off * s, cy + r * s - off * c)
    return ridge - 0.5 * (s1 + s2), ridge


def protect(mask, nonstellar, both, cat, fwhm, verbose=True):
    """Setzt die Maske auf allen Spikes heller Sterne auf 0. Gibt (maske, info) zurueck."""
    lumN = np.asarray(nonstellar, dtype=np.float32)
    lumB = np.asarray(both, dtype=np.float32)
    if lumN.ndim == 3:
        lumN = lumN.mean(axis=0)
        lumB = lumB.mean(axis=0)
    h, w = lumN.shape
    sub = lumN[::5, ::5]
    bg = float(np.median(sub))
    sigma = float(1.4826 * np.median(np.abs(sub - bg)))
    st = cat["stars"]
    bright = (st["mag"] < 13.0) | (st["sure"] & (st["mag"] < 15.5))
    idx = np.flatnonzero(bright & (st["x"] > 0) & (st["x"] < w) & (st["y"] > 0) & (st["y"] < h))
    sel = {k: v[idx] for k, v in st.items()}
    theta = find_angle(lumN, sel, fwhm, bg, sigma, verbose)
    if theta is None:
        if verbose:
            print("  Keine Beugungsspikes gefunden - nichts zu schuetzen")
        return mask, dict(spikes=0)

    sampN, sampB = _sampler(lumN), _sampler(lumB)
    spike = np.zeros((h, w), dtype=np.float32)
    n_star, longest = 0, 0.0
    for j in np.argsort(sel["mag"]):
        cx, cy, mag = sel["x"][j], sel["y"][j], sel["mag"][j]
        lmax = float(np.clip(40 * fwhm * 10 ** (0.2 * (14.0 - mag)), 20 * fwhm, 3000))
        r = np.arange(1.0, lmax, 1.0)
        off = 4.0 * fwhm
        # Rauschen des Spike-Signals an spikefreien Richtungen (45 Grad versetzt)
        noise = []
        for k in range(4):
            e, _ = _arm_profile(sampN, cx, cy, np.radians(theta + 45 + 90 * k), r, off)
            noise.append(_smooth(e, 2 * fwhm))
        noise = np.concatenate(noise)
        noise = noise[np.isfinite(noise)]
        if len(noise) < 50:
            continue
        s_e = 1.4826 * np.median(np.abs(noise - np.median(noise))) + 1e-12
        # Wo stimmen Stern- und Nicht-Stern-Version ausserhalb der Spikes wieder ueberein? Dort beginnt der Schutz.
        dN, dB = [], []
        for k in range(4):
            phi = np.radians(theta + 45 + 90 * k)
            c, s = np.cos(phi), np.sin(phi)
            dN.append(sampN(cx + r * c, cy + r * s))
            dB.append(sampB(cx + r * c, cy + r * s))
        dN, dB = np.nanmean(dN, axis=0), np.nanmean(dB, axis=0)
        same = np.abs(dB - dN) < 0.1 * np.maximum(dN - bg, 0) + 2 * sigma
        run = int(3 * fwhm)
        r0 = 1.5 * fwhm
        for t in range(int(1.5 * fwhm), len(r) - run):
            if same[t:t + run].all():
                r0 = r[t]
                break
        arms = []
        for k in range(4):
            phi = np.radians(theta + 90 * k)
            e, ridge = _arm_profile(sampN, cx, cy, phi, r, off)
            es = _smooth(e, 2 * fwhm)
            above = np.nan_to_num(es, nan=-1) > 3 * s_e
            start = int(max(r0, 3 * fwhm))
            if start >= len(r) or not above[start:start + int(4 * fwhm)].any():
                continue
            gap = int(8 * fwhm)
            last = start
            t = start
            while t < len(r):
                if above[t]:
                    last = t
                elif t - last > gap:
                    break
                t += 1
            L = r[last]
            if L < r0 + 3 * fwhm:
                continue
            # Breite aus dem Querprofil
            rr = r[start:min(last, start + int(20 * fwhm))]
            vs = np.arange(-6 * fwhm, 6 * fwhm + 0.5, 0.5)
            c, s = np.cos(phi), np.sin(phi)
            prof = np.array([np.nanmean(sampN(cx + rr * c - v * s, cy + rr * s + v * c)) for v in vs])
            base = np.nanmedian(np.concatenate([prof[:4], prof[-4:]]))
            pk = np.nanmax(prof) - base
            wide = vs[(prof - base) > 0.2 * pk] if pk > 0 else np.array([fwhm])
            hw = float(np.clip(1.5 * (np.max(np.abs(wide)) if len(wide) else fwhm), 1.5 * fwhm, 8 * fwhm))
            arms.append((phi, L, hw))
        if not arms:
            continue
        n_star += 1
        for phi, L, hw in arms:
            longest = max(longest, L)
            feather = 1.5 * fwhm
            ext = L + 4 * fwhm
            c, s = np.cos(phi), np.sin(phi)
            ra_ = max(0.35 * r0, 3.0 * fwhm)      # Sternkern bleibt unangetastet; Schutz laeuft bis r1 ein
            r1 = max(r0, ra_ + 3.0 * fwhm)
            xs = [cx + ra_ * c, cx + ext * c]
            ys = [cy + ra_ * s, cy + ext * s]
            pad = hw + 3 * feather
            x0, x1 = int(max(0, min(xs) - pad)), int(min(w, max(xs) + pad + 1))
            y0, y1 = int(max(0, min(ys) - pad)), int(min(h, max(ys) + pad + 1))
            if x0 >= x1 or y0 >= y1:
                continue
            yy, xx = np.mgrid[y0:y1, x0:x1].astype(np.float32)
            u = (xx - cx) * c + (yy - cy) * s
            v = np.abs(-(xx - cx) * s + (yy - cy) * c)
            ramp_in = np.clip((u - ra_) / (r1 - ra_), 0, 1)
            ramp_in = ramp_in * ramp_in * (3 - 2 * ramp_in)   # weicher S-Uebergang
            ramp_out = np.clip((ext - u) / (0.15 * L + 4 * fwhm), 0, 1)
            across = np.where(v <= hw, 1.0, np.exp(-0.5 * ((v - hw) / feather) ** 2))
            wgt = (ramp_in * ramp_out * across).astype(np.float32)
            np.maximum(spike[y0:y1, x0:x1], wgt, out=spike[y0:y1, x0:x1])
    if verbose:
        print("  Beugungsspikes: Winkel %.1f / %.1f Grad, an %d Sternen geschuetzt (laengster %d px)"
              % (theta, theta + 90, n_star, longest))
    return mask * (1.0 - spike), dict(spikes=n_star, winkel=theta)
