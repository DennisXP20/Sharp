"""Katalogabfragen fuer ein Bildfeld: Gaia DR3, Legacy Survey, SDSS, Pan-STARRS, HyperLEDA.

Alle Abfragen senden nur Himmelskoordinaten. Ergebnisse werden als .ecsv im Cache-Ordner abgelegt.
"""
import os
import warnings

import numpy as np
from astropy.table import Table, vstack

warnings.filterwarnings("ignore")
try:  # Windows-Zertifikatsspeicher nutzen (ESA/NOIRLab liefern keine Zwischenzertifikate mit)
    import truststore
    truststore.inject_into_ssl()
except Exception:
    pass

DATALAB_TAP = "https://datalab.noirlab.edu/tap"
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap"
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"


def _cached(cache_dir, name, fn):
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, name + ".ecsv")
    if os.path.exists(path):
        return Table.read(path)
    try:
        t = fn()
    except Exception as e:  # Netz-/Serverfehler: Katalog fehlt dann einfach
        print("  ! %s nicht abrufbar: %s" % (name, str(e).splitlines()[0][:200]))
        return None
    t.write(path, overwrite=True)
    return t


def _tap(url, query, maxrec=2000000):
    import pyvo
    svc = pyvo.dal.TAPService(url)
    job = svc.submit_job(query, maxrec=maxrec)
    job.run()
    job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=900)
    job.raise_if_error()
    t = job.fetch_result().to_table()
    try:
        job.delete()
    except Exception:
        pass
    for c in t.colnames:  # maskierte Spalten vereinfachen
        if t[c].dtype.kind == "O":
            t[c] = t[c].astype(str)
    return t


def _tap_sync(url, query, timeout=900):
    import pyvo
    t = pyvo.dal.TAPService(url).run_sync(query, timeout=timeout, maxrec=2000000).to_table()
    for c in t.colnames:
        if t[c].dtype.kind == "O":
            t[c] = t[c].astype(str)
    return t


def _box(ra0, dec0, half_ra, half_dec):
    """RA/Dec-Grenzen eines Kastens (half_* in Grad am Himmel)."""
    d1, d2 = dec0 - half_dec, dec0 + half_dec
    cosd = np.cos(np.radians(min(89.0, max(abs(d1), abs(d2)))))
    hr = half_ra / cosd
    return ra0 - hr, ra0 + hr, d1, d2


def gaia(cache, ra0, dec0, half_ra, half_dec):
    r1, r2, d1, d2 = _box(ra0, dec0, half_ra, half_dec)
    q = ("SELECT source_id, ra, dec, pmra, pmdec, pmra_error, pmdec_error, parallax, parallax_over_error, "
         "phot_g_mean_mag, bp_rp, ruwe, astrometric_params_solved, astrometric_excess_noise, "
         "classprob_dsc_combmod_star, classprob_dsc_combmod_galaxy, classprob_dsc_combmod_quasar, "
         "in_galaxy_candidates, in_qso_candidates "
         "FROM gaiadr3.gaia_source WHERE ra BETWEEN %.6f AND %.6f AND dec BETWEEN %.6f AND %.6f" % (r1, r2, d1, d2))
    return _cached(cache, "gaia_dr3", lambda: _tap(GAIA_TAP, q))


def legacy_survey(cache, ra0, dec0, half_ra, half_dec):
    """Legacy Survey Tractor-Katalog (DR10, sonst DR9). Liefert Morphologie-Typ und Ellipsenform."""
    radius = float(np.hypot(half_ra, half_dec))
    cols = ("ra, dec, type, shape_r, shape_e1, shape_e2, sersic, flux_g, flux_r, flux_z, "
            "flux_ivar_r, fracflux_r, maskbits, ref_cat, ref_id")
    for dr in ("ls_dr10", "ls_dr9"):
        q = ("SELECT %s FROM %s.tractor WHERE 't'=Q3C_RADIAL_QUERY(ra, dec, %.6f, %.6f, %.4f) "
             "AND brick_primary = 1" % (cols, dr, ra0, dec0, radius))
        t = _cached(cache, dr, lambda q=q: _tap_sync(DATALAB_TAP, q))  # async haengt dort oft
        if t is not None and len(t):
            t["dr"] = dr
            return t
    return None


VIZIER_SERVERS = ("vizier.cds.unistra.fr", "vizier.cfa.harvard.edu")


def _vizier(catalog, columns, ra0, dec0, half_ra, half_dec, column_filters=None):
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    cosd = np.cos(np.radians(min(89.0, abs(dec0) + half_dec)))
    last = None
    for server in VIZIER_SERVERS:
        try:
            v = Vizier(columns=columns, row_limit=-1, vizier_server=server, timeout=300,
                       column_filters=column_filters or {})
            res = v.query_region(SkyCoord(ra0 * u.deg, dec0 * u.deg), width=2 * half_ra / cosd * u.deg,
                                 height=2 * half_dec * u.deg, catalog=catalog)
            t = res[0] if len(res) else Table()
            for c in t.colnames:
                if t[c].dtype.kind == "O":
                    t[c] = t[c].astype(str)
            return t
        except Exception as e:
            last = e
    raise last


def sdss(cache, ra0, dec0, half_ra, half_dec):
    """SDSS DR16 ueber VizieR (V/154): class 3 = Galaxie, 6 = Stern."""
    return _cached(cache, "sdss_dr16", lambda: _vizier(
        "V/154/sdss16", ["RA_ICRS", "DE_ICRS", "class", "rmag", "e_rmag", "rPrad", "mode"],
        ra0, dec0, half_ra, half_dec, {"mode": "=1"}))


def panstarrs(cache, ra0, dec0, half_ra, half_dec):
    """Pan-STARRS1 DR1 ueber VizieR (II/349): Stern/Galaxie ueber PSF- minus Kron-Helligkeit."""
    return _cached(cache, "panstarrs_dr1", lambda: _vizier(
        "II/349/ps1", ["RAJ2000", "DEJ2000", "rmag", "rKmag", "imag", "iKmag", "Nd"],
        ra0, dec0, half_ra, half_dec, {"Nd": ">2"}))


def hyperleda(cache, ra0, dec0, half_ra, half_dec):
    """PGC/HyperLEDA ueber VizieR (VII/237): grosse Galaxien mit Durchmesser, Achsenverhaeltnis, Positionswinkel."""
    return _cached(cache, "hyperleda_pgc", lambda: _vizier(
        "VII/237/pgc", ["PGC", "RAJ2000", "DEJ2000", "OType", "logD25", "logR25", "PA"],
        ra0, dec0, half_ra + 0.5, half_dec + 0.5))
