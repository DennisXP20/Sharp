"""Hilfsfunktionen: Autostretch fuer Vorschauen, Ausschnitte speichern."""
import numpy as np
from PIL import Image


def mtf(m, x):
    """PixInsight-Mitteltonfunktion."""
    return (m - 1) * x / ((2 * m - 1) * x - m)


def stretch_params(lum, target=0.2, shadows=-2.8):
    sub = np.asarray(lum[::7, ::7], dtype=np.float64)
    med = np.median(sub)
    mad = np.median(np.abs(sub - med)) * 1.4826
    lo = max(0.0, med + shadows * mad)
    m = mtf(target, (med - lo) / (1 - lo))
    return lo, m


def autostretch(img, params):
    lo, m = params
    x = np.clip((np.asarray(img, dtype=np.float64) - lo) / (1 - lo), 0, 1)
    return mtf(m, x)


def save_png(path, img, params, scale=1):
    """img: (H,W) oder (3,H,W) linear."""
    y = autostretch(img, params)
    if y.ndim == 3:
        y = np.moveaxis(y, 0, -1)
    im = Image.fromarray((np.clip(y, 0, 1) * 255).astype("u1"))
    if scale != 1:
        im = im.resize((im.width * scale, im.height * scale), Image.NEAREST)
    im.save(path)
