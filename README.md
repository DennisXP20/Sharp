# Katalogschaerfung (Sharp)

PixInsight-Skript: schärft lineare, astrometrisch gelöste Bilder mit Cosmic Clarity
(SetiAstroSuite Pro), lässt das **Sternmodell aber nur dort wirken, wo Gaia DR3 oder der
Legacy Survey wirklich einen Stern kennen**. Hintergrundgalaxien bekommen nur die
Nicht-Stern-Schärfung und werden nicht zu Sternen oder Punktketten.

## Voraussetzungen

- PixInsight 1.8.9 oder neuer, Windows
- SetiAstroSuite Pro unter `C:\Program Files\SetiAstroSuitePro` (liefert Cosmic Clarity)
- Python 3.11 oder neuer
- Internetzugang (Katalogabfragen bei ESA/Gaia, NOIRLab, VizieR; es werden nur Himmelskoordinaten gesendet)

## Einrichten

Python-Umgebung anlegen (einmalig, in einer Eingabeaufforderung):

```
python -m venv %USERPROFILE%\Katalogschaerfung\venv
%USERPROFILE%\Katalogschaerfung\venv\Scripts\python -m pip install -r requirements.txt
```

Skript in PixInsight einbinden, entweder

- **über das Repository:** Resources > Updates > Manage Repositories > Add:
  `https://raw.githubusercontent.com/DennisXP20/Sharp/main/repository/`
  dann Resources > Updates > Check for Updates und PixInsight neu starten, oder
- **von Hand:** diesen Ordner herunterladen, Script > Feature Scripts > Add und den Ordner wählen.

Danach: Script > Utilities > Katalogschaerfung.

## Benutzung

Bild muss **linear** sein und eine **astrometrische Lösung** (ImageSolver) haben.
Laufzeit bei 79 Megapixeln und RTX 5060 Ti: etwa 10 Minuten.

Ohne PixInsight: `python katalogschaerfung.py "bild.xisf"`

## Wie es funktioniert

1. Astrometrie des Bildes an Gaia-Sternen nachschärfen (typisch 0,4 px Restfehler).
2. Gaia DR3, Legacy Survey (sonst SDSS) und HyperLEDA für das Bildfeld abfragen.
3. Cosmic Clarity zweimal rechnen: "Non-Stellar Only" und "Both".
4. Die Differenz zeigt, wo das Sternmodell eingegriffen hat. Übernommen wird es nur an
   Stellen mit katalogbestätigtem Stern; Kataloggalaxien werden ausgespart.

## Grenzen

- Die Sternkorrektur (Option im Dialog) ist ungetestet.
- Ohne Legacy-Survey-Abdeckung (Teile des Südhimmels, Milchstraßenebene) ist die
  Galaxienerkennung deutlich schwächer; der SDSS-Ersatzweg ist ungetestet.
- Das Repository-Paket ist nicht signiert.
