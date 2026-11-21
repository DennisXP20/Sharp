"""Katalogschaerfung: Cosmic Clarity, aber das Sternmodell wirkt nur auf katalogbestaetigte Sterne.

Direkter Aufruf:
    python katalogschaerfung.py "bild.xisf"
  -> bild_katalog.xisf (Metadaten wie das Original) und bild_sternmaske.xisf

Aus PixInsight: Skript Katalogschaerfung.js (ruft dieses Programm auf).

Voraussetzung: lineares Bild mit astrometrischer Loesung (PixInsight ImageSolver), XISF Float32.
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from katsharp import astrometry, catalogs, classify, mask as maskmod, spikes  # noqa: E402
from katsharp.xisf_io import XisfImage, write_xisf  # noqa: E402

SASPRO = r"C:\Program Files\SetiAstroSuitePro\setiastrosuitepro.exe"


def say(*a):
    print(*a, flush=True)


def cc_args(mode, a):
    c = ["--gpu" if not a.kein_gpu else "--no-gpu", "--temp-stretch", "--target-median", "0.25",
         "--chunk-size", "256", "--overlap", "64", "--sharpening-mode", mode,
         "--stellar-amount", "%.2f" % a.stern_staerke, "--nonstellar-amount", "%.2f" % a.nichtstern_staerke]
    c += ["--auto-psf"] if a.auto_psf else ["--no-auto-psf", "--nonstellar-psf", "%.2f" % a.nichtstern_psf]
    if a.korrektur and mode != "Non-Stellar Only":
        c += ["--stellar-correct-mode", "correct_sharpen"]
    return c


def run_cc(work, src, mode, a):
    """Ruft Cosmic Clarity auf; Ergebnis wird je Parametersatz zwischengespeichert."""
    args = cc_args(mode, a)
    key = hashlib.md5(" ".join(args).encode()).hexdigest()[:8]
    dst = os.path.join(work, "cc_%s_%s.fits" % (mode.lower().replace(" ", "").replace("-", ""), key))
    if os.path.exists(dst):
        say("  %s: vorhanden, wird wiederverwendet" % mode)
        return dst
    if not os.path.exists(a.saspro):
        sys.exit("SetiAstroSuite Pro nicht gefunden: %s" % a.saspro)
    tmp = dst + ".tmp.fits"
    t0 = time.time()
    r = subprocess.run([a.saspro, "cc", "sharpen", "-i", src, "-o", tmp] + args,
                       capture_output=True, text=True, errors="replace")
    if r.returncode != 0 or not os.path.exists(tmp):
        sys.exit("Cosmic Clarity (%s) fehlgeschlagen:\n%s" % (mode, (r.stdout + r.stderr)[-1500:]))
    os.replace(tmp, dst)
    say("  %s: fertig nach %.0f s" % (mode, time.time() - t0))
    return dst


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bild")
    p.add_argument("--stern-staerke", type=float, default=0.90, help="Cosmic Clarity Stellar Amount")
    p.add_argument("--nichtstern-staerke", type=float, default=0.50, help="Cosmic Clarity Non-Stellar Amount")
    p.add_argument("--nichtstern-psf", type=float, default=3.0, help="Cosmic Clarity Non-Stellar PSF")
    p.add_argument("--auto-psf", action="store_true")
    p.add_argument("--korrektur", action="store_true", help="Sternkorrektur (Aberration) vor dem Schaerfen")
    p.add_argument("--galaxiensterne", choices=["sicher", "alle"], default="sicher",
                   help="in grossen Galaxien nur Gaia-sichere Sterne schaerfen (sicher) oder alle Punktquellen (alle)")
    p.add_argument("--schwelle", type=float, default=1.0, help="Empfindlichkeit der Stellen-Erkennung in Sigma")
    p.add_argument("--spikes-schaerfen", action="store_true",
                   help="Beugungsspikes dem Sternmodell ueberlassen (Standard: Spikes bleiben wie im Nicht-Stern-Ergebnis)")
    p.add_argument("--kein-gpu", action="store_true")
    p.add_argument("--saspro", default=SASPRO)
    p.add_argument("--neu", action="store_true", help="Zwischenergebnisse verwerfen")
    # Aufruf aus PixInsight:
    p.add_argument("--wcs-punkte", help="Textdatei 'x y ra dec' (PixInsight-Bildkoordinaten) statt Loesung aus der Datei")
    p.add_argument("--datum", help="Aufnahmedatum JJJJ-MM-TT (sonst DATE-OBS aus der Datei)")
    p.add_argument("--ausgabe", help="Ergebnisdatei (einfaches XISF ohne Metadaten)")
    p.add_argument("--maske-ausgabe", help="Datei fuer die Sternmaske")
    p.add_argument("--arbeit", help="Ordner fuer Zwischenergebnisse")
    p.add_argument("--katalog-cache", help="Ordner fuer heruntergeladene Kataloge")
    a = p.parse_args()

    from astropy.io import fits
    src = os.path.abspath(a.bild)
    stem = os.path.splitext(src)[0]
    work = os.path.abspath(a.arbeit) if a.arbeit else stem + "_arbeit"
    if a.neu and os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work, exist_ok=True)

    say("1/5 Bild lesen")
    im = XisfImage(src)
    data = im.data()
    h, w = im.height, im.width
    lum = np.asarray(data, dtype=np.float32).mean(axis=0)

    say("2/5 Astrometrie pruefen und Kataloge abfragen")
    if a.wcs_punkte:
        wcs = astrometry.wcs_from_points(np.loadtxt(a.wcs_punkte), (h, w))
    else:
        wcs = astrometry.wcs_from_pixinsight(im)
    corners = wcs.all_pix2world([[0, 0], [w - 1, 0], [0, h - 1], [w - 1, h - 1], [w / 2, h / 2]], 0)
    ra0, dec0 = corners[4]
    half_dec = np.abs(corners[:4, 1] - dec0).max() + 0.02
    dra = (corners[:4, 0] - ra0 + 180) % 360 - 180
    half_ra = (np.abs(dra) * np.cos(np.radians(corners[:4, 1]))).max() + 0.02
    if a.katalog_cache:   # je Himmelsfeld ein Unterordner, damit wiederholte Laeufe nichts neu laden
        cache = os.path.join(a.katalog_cache, "ra%.2f_dec%+.2f_%.2fx%.2f" % (ra0, dec0, half_ra, half_dec))
    else:
        cache = os.path.join(work, "kataloge")
    field = (cache, float(ra0), float(dec0), float(half_ra), float(half_dec))
    say("  Feld: RA %.3f, Dec %+.3f, %.2f x %.2f Grad" % (ra0, dec0, 2 * half_ra, 2 * half_dec))
    gaia = catalogs.gaia(*field)
    if gaia is None or not len(gaia):
        sys.exit("Gaia DR3 nicht erreichbar - ohne Sternkatalog geht es nicht.")
    date = a.datum or im.fits_keyword("DATE-OBS") or ""
    try:
        obs_year = int(date[:4]) + (int(date[5:7]) - 1) / 12.0
    except ValueError:
        obs_year = time.gmtime().tm_year + 0.0
        say("  ! Kein Aufnahmedatum gefunden - Eigenbewegungen werden auf %d gerechnet" % obs_year)
    ls = catalogs.legacy_survey(*field)
    sdss = catalogs.sdss(*field) if ls is None else None
    pgc = catalogs.hyperleda(*field)

    F = classify._f
    g_mag = F(gaia["phot_g_mean_mag"], 99)
    dt = obs_year - classify.GAIA_EPOCH
    g_dec = F(gaia["dec"])
    g_ra = F(gaia["ra"]) + F(gaia["pmra"], 0) * dt / 3.6e6 / np.cos(np.radians(g_dec))
    g_dec = g_dec + F(gaia["pmdec"], 0) * dt / 3.6e6
    gx, gy = astrometry.world2pix(wcs, g_ra, g_dec)
    sel = (gx > 20) & (gx < w - 20) & (gy > 20) & (gy < h - 20) & (g_mag > 11) & (g_mag < 17.5)
    if sel.sum() < 30:
        sys.exit("Zu wenige Gaia-Sterne im Bild gefunden (%d) - stimmt die astrometrische Loesung?" % sel.sum())
    fwhm = astrometry.measure_fwhm(lum, wcs, g_ra[sel & (g_mag > 14)], g_dec[sel & (g_mag > 14)])
    try:
        wcs = astrometry.refine_wcs(lum, wcs, g_ra[sel], g_dec[sel], half=int(np.clip(round(fwhm), 7, 15)))
    except ValueError as e:
        sys.exit(str(e))
    fwhm = astrometry.measure_fwhm(lum, wcs, g_ra[sel & (g_mag > 14)], g_dec[sel & (g_mag > 14)])
    say("  Stern-FWHM: %.1f px" % fwhm)
    cat = classify.build(wcs, (h, w), obs_year, gaia, ls, sdss, pgc)
    del lum

    say("3/5 Cosmic Clarity (zwei Durchlaeufe ueber das ganze Bild - das dauert einige Minuten)")
    f_in = os.path.join(work, "eingabe.fits")
    if not os.path.exists(f_in) or os.path.getmtime(f_in) < os.path.getmtime(src):
        for old in os.listdir(work):          # neue Eingabe -> alte Cosmic-Clarity-Ergebnisse sind ungueltig
            if old.startswith("cc_") and old.endswith(".fits"):
                os.remove(os.path.join(work, old))
        fits.PrimaryHDU(np.asarray(data, dtype=np.float32)).writeto(f_in, overwrite=True)
    f_ns = run_cc(work, f_in, "Non-Stellar Only", a)
    f_both = run_cc(work, f_in, "Both", a)

    say("4/5 Sternmaske aus Katalog")

    def load(path):
        x = fits.getdata(path).astype(np.float32)
        if x.ndim == 2:
            x = x[None]
        if x.shape[0] not in (1, 3) and x.shape[-1] in (1, 3):
            x = np.moveaxis(x, -1, 0)
        if x.shape[0] != data.shape[0] and data.shape[0] == 1:
            x = x.mean(axis=0, keepdims=True)
        if x.shape != data.shape:
            sys.exit("Cosmic-Clarity-Ausgabe hat andere Bildgroesse als die Eingabe: %s statt %s" % (x.shape, data.shape))
        return x

    ns, both = load(f_ns), load(f_both)
    m, stats = maskmod.star_mask(both, ns, cat, fwhm, k_sigma=a.schwelle, big_galaxy_stars=a.galaxiensterne)
    if not a.spikes_schaerfen:
        m, _ = spikes.protect(m, ns, both, cat, fwhm)

    say("5/5 Zusammensetzen und speichern")
    for c in range(ns.shape[0]):
        ns[c] += m * (both[c] - ns[c])
    del both
    np.clip(ns, 0, 1, out=ns)
    if a.ausgabe:
        out = os.path.abspath(a.ausgabe)
        write_xisf(out, ns, image_id="katalogschaerfung")
    else:
        out = stem + "_katalog.xisf"
        if not im.plain_float32:
            write_xisf(out, ns, image_id="katalogschaerfung")
        else:
            shutil.copyfile(src, out)      # Kopie des Originals behaelt alle Metadaten (Astrometrie, SPCC, ...)
            with open(out, "r+b") as f:
                f.seek(im.offset)
                f.write(ns.astype("<f4").tobytes())
    mout = os.path.abspath(a.maske_ausgabe) if a.maske_ausgabe else stem + "_sternmaske.xisf"
    write_xisf(mout, m, image_id="sternmaske")
    say("ERGEBNIS: %(uebernommen)d Sterne geschaerft, %(verworfen_galaxie)d Galaxien und "
        "%(verworfen_unbekannt)d unbekannte Objekte vor dem Sternmodell geschuetzt" % stats)
    say("Fertig:\n  %s\n  %s" % (out, mout))


if __name__ == "__main__":
    main()
