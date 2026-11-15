"""Minimaler XISF-Leser/-Schreiber (nur unkomprimierte Float32-Bilder, wie PixInsight sie standardmaessig speichert)."""
import base64
import re
import struct

import numpy as np


class XisfImage:
    def __init__(self, path):
        self.path = path
        with open(path, "rb") as f:
            if f.read(8) != b"XISF0100":
                raise ValueError("Keine XISF-Datei: %s" % path)
            n = struct.unpack("<I", f.read(4))[0]
            f.read(4)
            self.xml = f.read(n).decode("utf-8", "replace")
        m = re.search(r"<Image ([^>]*)>", self.xml)
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        self.compression = attrs.get("compression")
        self.compressed = self.compression is not None
        formats = {"Float32": ("<f4", 1.0), "Float64": ("<f8", 1.0), "UInt16": ("<u2", 65535.0),
                   "UInt8": ("u1", 255.0), "UInt32": ("<u4", 4294967295.0)}
        if attrs.get("sampleFormat") not in formats:
            raise ValueError("Pixelformat %s wird nicht unterstuetzt" % attrs.get("sampleFormat"))
        self.dtype, self.scale = formats[attrs["sampleFormat"]]
        w, h, c = (int(v) for v in attrs["geometry"].split(":"))
        self.width, self.height, self.channels = w, h, c
        _, off, size = attrs["location"].split(":")
        self.offset, self.size = int(off), int(size)

    @property
    def plain_float32(self):
        """True, wenn die Pixel unkomprimiert als Float32 in der Datei stehen (dann kann man sie direkt ersetzen)."""
        return not self.compressed and self.dtype == "<f4"

    def data(self):
        """Bilddaten als Float32 0..1, Form (Kanaele, Hoehe, Breite), Zeile 0 = oben."""
        shape = (self.channels, self.height, self.width)
        if self.plain_float32:
            return np.memmap(self.path, dtype="<f4", mode="r", offset=self.offset, shape=shape)
        if not self.compressed:
            raw = np.memmap(self.path, dtype=self.dtype, mode="r", offset=self.offset, shape=shape)
        else:
            parts = self.compression.split(":")
            codec = parts[0].lower()
            if not codec.startswith("zlib"):
                raise ValueError("XISF-Kompression '%s' wird nicht unterstuetzt - bitte ohne Kompression oder mit zlib speichern." % codec)
            import zlib
            with open(self.path, "rb") as f:
                f.seek(self.offset)
                buf = zlib.decompress(f.read(self.size))
            if codec.endswith("+sh"):   # Byte-Shuffling rueckgaengig machen
                item = int(parts[2])
                n = len(buf) // item
                buf = np.frombuffer(buf, dtype=np.uint8)[:n * item].reshape(item, n).T.tobytes()
            raw = np.frombuffer(buf, dtype=self.dtype).reshape(shape)
        out = np.asarray(raw, dtype=np.float32)
        if self.scale != 1.0:
            out /= self.scale
        return out

    def _prop(self, pid):
        m = re.search(r'<Property id="%s"([^>]*?)(?:/>|>(.*?)</Property>)' % re.escape(pid), self.xml, re.S)
        if not m:
            return None, None
        return m.group(1), m.group(2)

    def prop(self, pid):
        a, body = self._prop(pid)
        if a is None:
            return None
        v = re.search(r'value="([^"]*)"', a)
        if v:
            return v.group(1)
        if "inline:base64" in a:
            return np.frombuffer(base64.b64decode(body), dtype="<f8")
        m = re.search(r"attachment:(\d+):(\d+)", a)
        if m:
            with open(self.path, "rb") as f:
                f.seek(int(m.group(1)))
                return np.frombuffer(f.read(int(m.group(2))), dtype="<f8")
        return body

    def fits_keyword(self, name):
        m = re.search(r'<FITSKeyword name="%s" value="([^"]*)"' % re.escape(name), self.xml)
        return m.group(1).strip("' ") if m else None


def write_xisf(path, arr, image_id="image"):
    """Schreibt (C,H,W)- oder (H,W)-Float32 als unkomprimiertes XISF."""
    arr = np.asarray(arr, dtype="<f4")
    if arr.ndim == 2:
        arr = arr[None]
    c, h, w = arr.shape
    cs = "RGB" if c == 3 else "Gray"

    def header(offset):
        return ('<?xml version="1.0" encoding="UTF-8"?>'
                '<xisf version="1.0" xmlns="http://www.pixinsight.com/xisf" '
                'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                'xsi:schemaLocation="http://www.pixinsight.com/xisf http://pixinsight.com/xisf/xisf-1.0.xsd">'
                '<Image id="%s" geometry="%d:%d:%d" sampleFormat="Float32" bounds="0:1" colorSpace="%s" '
                'location="attachment:%d:%d"/></xisf>' % (image_id, w, h, c, cs, offset, arr.nbytes)).encode("utf-8")

    offset = 16 + len(header(10 ** 12)) + 64
    offset = (offset + 4095) // 4096 * 4096
    hdr = header(offset)
    with open(path, "wb") as f:
        f.write(b"XISF0100" + struct.pack("<I", len(hdr)) + b"\0\0\0\0" + hdr)
        f.write(b"\0" * (offset - f.tell()))
        f.write(arr.tobytes())
