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
    px, py = wcs.all_world2pix(ra, dec, 0)
    for rnd, hw in enumerate((half + 5, half, half)):
        mx, my, keep = [], [], []
        for i, (x, y) in enumerate(zip(px, py)):
            c = centroid(lum, x, y, half=hw)
            if c is None:
                continue
            if np.hypot(c[0] - x, c[1] - y) > hw * 0.8:
                continue
            mx.append(c[0]); my.append(c[1]); keep.append(i)
        mx, my, keep = np.array(mx), np.array(my), np.array(keep)
        for _ in range(3):
            sky = SkyCoord(ra[keep] * u.deg, dec[keep] * u.deg)
            new = fit_wcs_from_points((mx, my), sky, projection=wcs, sip_degree=sip_degree)
            fx, fy = new.all_world2pix(ra[keep], dec[keep], 0)
            r = np.hypot(fx - mx, fy - my)
            good = r < max(clip_px, 3 * np.median(r))
            if good.all():
                break
            mx, my, keep = mx[good], my[good], keep[good]
        wcs = new
        wcs.array_shape = lum.shape
        px, py = wcs.all_world2pix(ra, dec, 0)
        if verbose:
            print("  WCS-Runde %d: %d Sterne, Restfehler median %.2f px, 95%% %.2f px"
                  % (rnd + 1, len(keep), np.median(r), np.percentile(r, 95)))
    return wcs


def measure_fwhm(lum, wcs, ra, dec, half=12, max_stars=400):
    """Mittlere Stern-FWHM in Pixel (Flaeche ueber halber Hoehe)."""
    px, py = wcs.all_world2pix(ra, dec, 0)
    h, w = lum.shape
    out = []
    for x, y in list(zip(px, py))[:max_stars]:
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
