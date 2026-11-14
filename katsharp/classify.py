"""Katalogobjekte zu Sternen (Punktquellen) und Galaxien (ausgedehnt) zusammenfuehren."""
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

GAIA_EPOCH = 2016.0


def _f(col, fill=np.nan):
    a = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(col, dtype=float)), fill)
    return np.asarray(a, dtype=float)


def _ellipse_to_pixels(wcs, ra, dec, a_arcsec, pa_deg):
    """Grosse Halbachse in Pixel und Winkel der Hauptachse im Bild (rad, gegen +x)."""
    c = SkyCoord(ra * u.deg, dec * u.deg)
    step = np.maximum(a_arcsec, 1.0)
    c2 = c.directional_offset_by(np.nan_to_num(pa_deg) * u.deg, step * u.arcsec)
    x0, y0 = wcs.all_world2pix(ra, dec, 0)
    x1, y1 = wcs.all_world2pix(c2.ra.deg, c2.dec.deg, 0)
    scale = np.hypot(x1 - x0, y1 - y0) / step  # px pro Bogensekunde
    return x0, y0, a_arcsec * scale, np.arctan2(y1 - y0, x1 - x0)


def build(wcs, shape, obs_year, gaia, ls=None, sdss=None, pgc=None, big_arcsec=25.0, verbose=True):
    """Liefert dict mit 'stars' und 'galaxies' (jeweils dict aus Arrays, Pixelkoordinaten im Bild).

    stars:    x, y, mag, sure (Gaia-Astrometrie signifikant), in_big (liegt in grosser Galaxie)
    galaxies: x, y, a (grosse Halbachse px), q, theta, mag, big
    """
    h, w = shape
    margin = 50

    def inside(x, y):
        return (x > -margin) & (x < w + margin) & (y > -margin) & (y < h + margin)

    # ---------- Gaia ----------
    g_ra0, g_dec0 = _f(gaia["ra"]), _f(gaia["dec"])
    pmra, pmde = _f(gaia["pmra"], 0), _f(gaia["pmdec"], 0)
    dt = obs_year - GAIA_EPOCH
    g_ra = g_ra0 + pmra * dt / 3.6e6 / np.cos(np.radians(g_dec0))
    g_dec = g_dec0 + pmde * dt / 3.6e6
    g_mag = _f(gaia["phot_g_mean_mag"], 21.5)
    pm_sig = np.hypot(pmra / _f(gaia["pmra_error"], np.inf), pmde / _f(gaia["pmdec_error"], np.inf))
    plx_sig = _f(gaia["parallax_over_error"], 0)
    g_sure = (pm_sig > 4) | (plx_sig > 4)
    g_galflag = (np.asarray(gaia["in_galaxy_candidates"]).astype(str) == "True") | \
                (_f(gaia["classprob_dsc_combmod_galaxy"], 0) > 0.5)
    g_point = np.ones(len(gaia), bool)   # wird unten durch LS/SDSS ggf. widerlegt
    g_point[g_galflag & ~g_sure] = False

    gal = dict(ra=[], dec=[], a=[], q=[], pa=[], mag=[])
    star_extra = dict(ra=[], dec=[], mag=[])
    gc = SkyCoord(g_ra0 * u.deg, g_dec0 * u.deg)

    def add_extended(ra, dec, a, q, pa, mag):
        gal["ra"].append(ra); gal["dec"].append(dec); gal["a"].append(a)
        gal["q"].append(q); gal["pa"].append(pa); gal["mag"].append(mag)

    def veto_gaia(coords):
        """Gaia-Quellen ohne sichere Astrometrie, die laut tieferem Katalog ausgedehnt sind, sind keine Sterne."""
        if len(coords) == 0:
            return
        idx, d2d, _ = gc.match_to_catalog_sky(coords)
        hit = (d2d.arcsec < 1.0) & ~g_sure
        g_point[hit] = False

    def unmatched(coords):
        if len(coords) == 0:
            return np.zeros(0, bool)
        idx, d2d, _ = coords.match_to_catalog_sky(gc)
        return d2d.arcsec > 1.0

    # ---------- Legacy Survey ----------
    if ls is not None and len(ls):
        typ = np.char.strip(np.asarray(ls["type"]).astype(str))
        ra, dec = _f(ls["ra"]), _f(ls["dec"])
        flux = _f(ls["flux_r"], 0)
        for alt in ("flux_g", "flux_z"):
            bad = ~(flux > 0)
            flux[bad] = _f(ls[alt], 0)[bad]
        mag = np.where(flux > 0, 22.5 - 2.5 * np.log10(np.maximum(flux, 1e-9)), 26.0)
        ext = np.isin(typ, ["REX", "EXP", "DEV", "SER"])
        # Scheingalaxien in den Fluegeln heller Sterne: Licht stammt ueberwiegend vom Nachbarn
        ext &= ~(_f(ls["fracflux_r"], 0) > 5)
        psf = typ == "PSF"
        e1, e2 = _f(ls["shape_e1"], 0), _f(ls["shape_e2"], 0)
        e = np.minimum(np.hypot(e1, e2), 0.95)
        q = (1 - e) / (1 + e)
        pa = np.degrees(0.5 * np.arctan2(e2, e1))
        r = _f(ls["shape_r"], 0)
        add_extended(ra[ext], dec[ext], r[ext], q[ext], pa[ext], mag[ext])
        veto_gaia(SkyCoord(ra[ext] * u.deg, dec[ext] * u.deg))
        pc = SkyCoord(ra[psf] * u.deg, dec[psf] * u.deg)
        um = unmatched(pc)
        star_extra["ra"].append(ra[psf][um]); star_extra["dec"].append(dec[psf][um]); star_extra["mag"].append(mag[psf][um])
        if verbose:
            print("  Legacy Survey: %d ausgedehnt, %d punktfoermig" % (ext.sum(), psf.sum()))
    # ---------- SDSS (nur wenn kein Legacy Survey) ----------
    elif sdss is not None and len(sdss):
        ra, dec = _f(sdss["RA_ICRS"]), _f(sdss["DE_ICRS"])
        cls = _f(sdss["class"], 0)
        mag = _f(sdss["rmag"], 25)
        ext = cls == 3
        psf = cls == 6
        add_extended(ra[ext], dec[ext], _f(sdss["rPrad"], 1.5)[ext] / 2, np.ones(ext.sum()), np.zeros(ext.sum()), mag[ext])
        veto_gaia(SkyCoord(ra[ext] * u.deg, dec[ext] * u.deg))
        pc = SkyCoord(ra[psf] * u.deg, dec[psf] * u.deg)
        um = unmatched(pc)
        star_extra["ra"].append(ra[psf][um]); star_extra["dec"].append(dec[psf][um]); star_extra["mag"].append(mag[psf][um])
        if verbose:
            print("  SDSS: %d Galaxien, %d Sterne" % (ext.sum(), psf.sum()))
    elif verbose:
        print("  ! Kein tiefer Galaxienkatalog fuer dieses Feld - nur Gaia und HyperLEDA.")

    # ---------- HyperLEDA: grosse Galaxien ----------
    if pgc is not None and len(pgc):
        c = SkyCoord(np.asarray(pgc["RAJ2000"]).astype(str), np.asarray(pgc["DEJ2000"]).astype(str),
                     unit=(u.hourangle, u.deg))
        ld, lr, pa = _f(pgc["logD25"]), _f(pgc["logR25"], 0), _f(pgc["PA"])
        ok = np.isfinite(ld)
        a = 3.0 * 10 ** ld[ok]          # D25/2 in Bogensekunden
        q = 10 ** -lr[ok]
        q[~np.isfinite(pa[ok])] = 1.0   # ohne Positionswinkel: Kreis
        add_extended(c.ra.deg[ok], c.dec.deg[ok], a * 0.6, q, np.nan_to_num(pa[ok]), np.full(ok.sum(), 14.0))
        veto_gaia(c[ok])

    # ---------- zusammenfuehren ----------
    for k in gal:
        gal[k] = np.concatenate(gal[k]) if gal[k] else np.zeros(0)
    gx, gy, ga, gth = _ellipse_to_pixels(wcs, gal["ra"], gal["dec"], gal["a"], gal["pa"])
    keep = inside(gx, gy)
    galaxies = dict(x=gx[keep], y=gy[keep], a=ga[keep], q=gal["q"][keep], theta=gth[keep],
                    mag=gal["mag"][keep], a_arcsec=gal["a"][keep])
    galaxies["big"] = galaxies["a_arcsec"] > big_arcsec

    sx, sy = wcs.all_world2pix(g_ra, g_dec, 0)
    s = dict(x=sx[g_point], y=sy[g_point], mag=g_mag[g_point], sure=g_sure[g_point])
    if star_extra["ra"]:
        era, edec, emag = (np.concatenate(star_extra[k]) for k in ("ra", "dec", "mag"))
        ex, ey = wcs.all_world2pix(era, edec, 0)
        s = dict(x=np.concatenate([s["x"], ex]), y=np.concatenate([s["y"], ey]),
                 mag=np.concatenate([s["mag"], emag]), sure=np.concatenate([s["sure"], np.zeros(len(ex), bool)]))
    keep = inside(s["x"], s["y"])
    stars = {k: v[keep] for k, v in s.items()}

    # Punktquellen innerhalb grosser Galaxien (HII-Regionen, Sternhaufen, Kerne) zaehlen nur mit sicherer Gaia-Astrometrie
    in_big = np.zeros(len(stars["x"]), bool)
    for x, y, a, q, th in zip(*(galaxies[k][galaxies["big"]] for k in ("x", "y", "a", "q", "theta"))):
        rad = 3.0 * a   # Halblichtradius -> sichtbare Ausdehnung
        dx, dy = stars["x"] - x, stars["y"] - y
        near = (np.abs(dx) < rad) & (np.abs(dy) < rad)
        u_ = dx[near] * np.cos(th) + dy[near] * np.sin(th)
        v_ = -dx[near] * np.sin(th) + dy[near] * np.cos(th)
        in_big[np.flatnonzero(near)[(u_ / rad) ** 2 + (v_ / (rad * max(q, 0.15))) ** 2 < 1]] = True
    stars["in_big"] = in_big
    if verbose:
        print("  Im Bild: %d Sterne (%d mit sicherer Gaia-Astrometrie, %d in grossen Galaxien), %d Galaxien (%d gross)"
              % (len(stars["x"]), stars["sure"].sum(), in_big.sum(), len(galaxies["x"]), galaxies["big"].sum()))
    return dict(stars=stars, galaxies=galaxies)
