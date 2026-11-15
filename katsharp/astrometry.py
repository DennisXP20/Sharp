"""WCS aus der PixInsight-Loesung aufbauen und an Gaia-Sternen im Bild nachschaerfen."""
import numpy as np
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord
import astropy.units as u


def wcs_from_pixinsight(im):
    """Lineare TAN-WCS aus den PCL:AstrometricSolution-Eigenschaften.

    Pixelkoordinaten: Array-Konvention (x=Spalte, y=Zeile, 0-basiert, Zeile 0 oben, Pixelmitte ganzzahlig).
    """
    ref_c = im.prop("PCL:AstrometricSolution:ReferenceCelestialCoordinates")
    ref_i = im.prop("PCL:AstrometricSolution:ReferenceImageCoordinates")
    m = im.prop("PCL:AstrometricSolution:LinearTransformationMatrix")
    if ref_c is None or m is None:
        raise ValueError("Keine astrometrische Loesung in der Datei gefunden.")
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [float(ref_c[0]), float(ref_c[1])]
    # PixInsight: Pixelmitte bei +0.5; FITS: erste Pixelmitte = 1  ->  FITS = PI + 0.5
    w.wcs.crpix = [float(ref_i[0]) + 0.5, float(ref_i[1]) + 0.5]
    w.wcs.cd = np.array(m, dtype=float).reshape(2, 2)
    w.array_shape = (im.height, im.width)
    return w


def wcs_from_points(pts, shape, sip_degree=3):
    """WCS aus Stuetzpunkten 'x y ra dec' (x, y in PixInsight-Bildkoordinaten, Pixelmitte bei +0.5)."""
    pts = np.atleast_2d(pts)
    x, y = pts[:, 0] - 0.5, pts[:, 1] - 0.5
    sky = SkyCoord(pts[:, 2] * u.deg, pts[:, 3] * u.deg)
    w = fit_wcs_from_points((x, y), sky, projection="TAN", sip_degree=sip_degree)
    w.array_shape = shape
    return w


def world2pix(wcs, ra, dec, margin=400):
    """Himmel -> Pixel, robust: Objekte weit ausserhalb des Bildes bekommen NaN.

    Die Verzeichnungskorrektur (SIP) konvergiert ausserhalb des Bildes nicht; deshalb wird zuerst mit dem
    linearen Teil grob geprueft, was ueberhaupt in Bildnaehe liegt.
    """
    ra = np.atleast_1d(np.asarray(ra, dtype=float))
    dec = np.atleast_1d(np.asarray(dec, dtype=float))
    x = np.full(ra.shape, np.nan)
    y = np.full(ra.shape, np.nan)
    if not len(ra):
        return x, y
    h, w = wcs.array_shape
    lx, ly = wcs.wcs_world2pix(ra, dec, 0)
    near = (lx > -margin) & (lx < w + margin) & (ly > -margin) & (ly < h + margin)
    if near.any():
        try:
            x[near], y[near] = wcs.all_world2pix(ra[near], dec[near], 0, quiet=True)
        except Exception:
            x[near], y[near] = lx[near], ly[near]
    return x, y


def centroid(lum, x, y, half=7, iters=3):
    """Schwerpunkt in einem Fenster um (x, y); gibt (x, y, fluss, fwhm) zurueck oder None."""
    h, w = lum.shape
    cx, cy = x, y
    for _ in range(iters):
        ix, iy = int(round(cx)), int(round(cy))
        if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
            return None
        win = np.asarray(lum[iy - half:iy + half + 1, ix - half:ix + half + 1], dtype=np.float64)
        edge = np.concatenate([win[0], win[-1], win[1:-1, 0], win[1:-1, -1]])
        win = win - np.median(edge)
        if win[half - 3:half + 4, half - 3:half + 4].max() < 5 * 1.4826 * np.median(np.abs(edge - np.median(edge))):
            return None                       # kein deutlicher Stern an dieser Stelle
        win[win < 0] = 0
        s = win.sum()
        if s <= 0:
            return None
        yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
        dx, dy = (win * xx).sum() / s, (win * yy).sum() / s
        cx, cy = ix + dx, iy + dy
    var = (win * ((xx - dx) ** 2 + (yy - dy) ** 2)).sum() / s / 2
    return cx, cy, s, 2.3548 * np.sqrt(max(var, 1e-6))


def refine_wcs(lum, wcs, ra, dec, sip_degree=4, half=7, clip_px=1.5, verbose=True):
    """Passt die WCS (TAN + SIP) an gemessene Positionen bekannter Sterne an."""
    px, py = world2pix(wcs, ra, dec)
    for rnd, hw in enumerate((half + 5, half, half)):
        mx, my, keep = [], [], []
        for i, (x, y) in enumerate(zip(px, py)):
            if not (np.isfinite(x) and np.isfinite(y)):
                continue
            c = centroid(lum, x, y, half=hw)
            if c is None:
                continue
            if np.hypot(c[0] - x, c[1] - y) > hw * 0.8:
                continue
            mx.append(c[0]); my.append(c[1]); keep.append(i)
        mx, my, keep = np.array(mx), np.array(my), np.array(keep)
        n_cand = int(np.isfinite(px).sum())
        if rnd == 0 and len(keep) < max(15, 0.35 * n_cand):
            raise ValueError("Die astrometrische Loesung passt nicht zum Bild: an nur %d von %d Gaia-Sternpositionen "
                             "sitzt tatsaechlich ein Stern. Bitte das Bild neu loesen (ImageSolver)." % (len(keep), n_cand))
        for _ in range(3):
            sky = SkyCoord(ra[keep] * u.deg, dec[keep] * u.deg)
            new = fit_wcs_from_points((mx, my), sky, projection=wcs, sip_degree=sip_degree)
            fx, fy = new.all_world2pix(ra[keep], dec[keep], 0, quiet=True)
            r = np.hypot(fx - mx, fy - my)
            good = r < max(clip_px, 3 * np.median(r))
            if good.all():
                break
            mx, my, keep = mx[good], my[good], keep[good]
        wcs = new
        wcs.array_shape = lum.shape
        px, py = world2pix(wcs, ra, dec)
        if verbose:
            print("  WCS-Runde %d: %d Sterne, Restfehler median %.2f px, 95%% %.2f px"
                  % (rnd + 1, len(keep), np.median(r), np.percentile(r, 95)))
    return wcs


def measure_fwhm(lum, wcs, ra, dec, half=12, max_stars=400):
    """Mittlere Stern-FWHM in Pixel (Flaeche ueber halber Hoehe)."""
    px, py = world2pix(wcs, ra, dec)
    h, w = lum.shape
    out = []
    for x, y in list(zip(px, py))[:max_stars]:
        if not np.isfinite(x):
            continue
        ix, iy = int(round(x)), int(round(y))
        if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
            continue
        win = np.asarray(lum[iy - half:iy + half + 1, ix - half:ix + half + 1], dtype=np.float64)
        edge = np.concatenate([win[0], win[-1], win[1:-1, 0], win[1:-1, -1]])
        win = win - np.median(edge)
        peak = win[half - 2:half + 3, half - 2:half + 3].max()
        if peak <= 5 * edge.std():
            continue
        out.append(2 * np.sqrt((win > peak / 2).sum() / np.pi))
    return float(np.median(out)) if out else 5.0
