// ============================================================================
// Katalogschaerfung.js
//
// Schaerft ein lineares, astrometrisch geloestes Bild mit Cosmic Clarity
// (SetiAstroSuite Pro), laesst das Sternmodell aber nur dort wirken, wo Gaia DR3
// bzw. der Legacy Survey wirklich einen Stern kennen. Hintergrundgalaxien
// bekommen nur die Nicht-Stern-Schaerfung und werden nicht zu Sternen.
//
// Die Rechenarbeit macht katalogschaerfung.py (liegt neben diesem Skript).
// ============================================================================

#feature-id    Katalogschaerfung : Utilities > Katalogschaerfung
#feature-info  Cosmic Clarity mit Katalogpruefung: Sterne werden nur dort geschaerft, wo Gaia oder der Legacy Survey einen Stern kennen. Braucht eine astrometrische Loesung.

#include <pjsr/Sizer.jsh>
#include <pjsr/StdButton.jsh>
#include <pjsr/StdIcon.jsh>
#include <pjsr/NumericControl.jsh>
#include <pjsr/TextAlign.jsh>
#include <pjsr/UndoFlag.jsh>
#include <pjsr/DataType.jsh>

#define KS_TITLE   "Katalogschaerfung"
#define KS_VERSION "1.0"
#define KS_KEY     "Katalogschaerfung/"

/*
 * PixInsight liefert den eigenen Skriptpfad mit zerstoerten Umlauten (z.B. "sch?rfung").
 * Hier wird der echte Ordnername ueber die Verzeichnisliste wiedergefunden.
 */
function ksResolvePath( path )
{
   if ( path.indexOf( "�" ) < 0 )
      return path;
   var parts = path.split( "/" );
   var cur = parts[0];
   for ( var i = 1; i < parts.length; ++i )
   {
      var name = parts[i];
      if ( name.indexOf( "�" ) >= 0 )
      {
         var rx = new RegExp( "^" + name.replace( /[.*+?^${}()|[\]\\]/g, "\\$&" ).replace( /�+/g, ".{1,2}" ) + "$" );
         var found = null;
         var ff = new FileFind;
         if ( ff.begin( cur + "/*" ) )
         {
            do
            {
               if ( ff.isDirectory && rx.test( ff.name ) )
               {
                  found = ff.name;
                  break;
               }
            }
            while ( ff.next() );
            ff.end();
         }
         if ( found === null )
            return path;
         name = found;
      }
      cur += "/" + name;
   }
   return cur;
}

var KS_SCRIPT_DIR = ksResolvePath( File.extractDrive( #__FILE__ ) + File.extractDirectory( #__FILE__ ) );

// ----------------------------------------------------------------------------
// Parameter
// ----------------------------------------------------------------------------

var KSParams = {
   sternStaerke:      0.90,
   nichtsternStaerke: 0.50,
   nichtsternPSF:     3.00,
   autoPSF:           false,
   korrektur:         false,
   galaxiensterne:    0,        // 0 = nur Gaia-sichere Sterne in grossen Galaxien, 1 = alle Punktquellen
   maskeOeffnen:      false,
   neuesBild:         false,
   pythonPfad:        File.homeDirectory + "/Katalogschaerfung/venv/Scripts/python.exe",
   skriptPfad:        KS_SCRIPT_DIR + "/katalogschaerfung.py",

   load: function()
   {
      function rd( key, type, def )
      {
         var v = Settings.read( KS_KEY + key, type );
         return ( Settings.lastReadOK && v !== null && v !== undefined ) ? v : def;
      }
      this.sternStaerke      = rd( "sternStaerke",      DataType_Double,  this.sternStaerke );
      this.nichtsternStaerke = rd( "nichtsternStaerke", DataType_Double,  this.nichtsternStaerke );
      this.nichtsternPSF     = rd( "nichtsternPSF",     DataType_Double,  this.nichtsternPSF );
      this.autoPSF           = rd( "autoPSF",           DataType_Boolean, this.autoPSF );
      this.korrektur         = rd( "korrektur",         DataType_Boolean, this.korrektur );
      this.galaxiensterne    = rd( "galaxiensterne",    DataType_Int32,   this.galaxiensterne );
      this.maskeOeffnen      = rd( "maskeOeffnen",      DataType_Boolean, this.maskeOeffnen );
      this.neuesBild         = rd( "neuesBild",         DataType_Boolean, this.neuesBild );
      this.pythonPfad        = rd( "pythonPfad",        DataType_String,  this.pythonPfad );
      var sp                 = rd( "skriptPfad",        DataType_String,  "" );
      if ( sp.length > 0 && File.exists( sp ) )
         this.skriptPfad = sp;
   },

   save: function()
   {
      Settings.write( KS_KEY + "sternStaerke",      DataType_Double,  this.sternStaerke );
      Settings.write( KS_KEY + "nichtsternStaerke", DataType_Double,  this.nichtsternStaerke );
      Settings.write( KS_KEY + "nichtsternPSF",     DataType_Double,  this.nichtsternPSF );
      Settings.write( KS_KEY + "autoPSF",           DataType_Boolean, this.autoPSF );
      Settings.write( KS_KEY + "korrektur",         DataType_Boolean, this.korrektur );
      Settings.write( KS_KEY + "galaxiensterne",    DataType_Int32,   this.galaxiensterne );
      Settings.write( KS_KEY + "maskeOeffnen",      DataType_Boolean, this.maskeOeffnen );
      Settings.write( KS_KEY + "neuesBild",         DataType_Boolean, this.neuesBild );
      Settings.write( KS_KEY + "pythonPfad",        DataType_String,  this.pythonPfad );
      // Nur merken, wenn von Hand ein anderer Ort gewaehlt wurde - sonst gilt die Datei neben diesem Skript.
      Settings.write( KS_KEY + "skriptPfad",        DataType_String,
                      this.skriptPfad == KS_SCRIPT_DIR + "/katalogschaerfung.py" ? "" : this.skriptPfad );
   }
};

// ----------------------------------------------------------------------------
// Hilfsfunktionen
// ----------------------------------------------------------------------------

function ksTempDir()
{
   var d = File.systemTempDirectory + "/Katalogschaerfung";
   if ( !File.directoryExists( d ) )
      File.createDirectory( d, true );
   return d;
}

/*
 * Laufprotokoll: jede Stufe und die Ausgabe des Python-Programms landen in
 * <Temp>/Katalogschaerfung/letzter_lauf.log - hilfreich, wenn etwas schiefgeht.
 */
var ksTraceLines = [];
function ksTrace( text )
{
   ksTraceLines.push( text );
   try { File.writeTextFile( ksTempDir() + "/letzter_lauf.log", ksTraceLines.join( "\n" ) + "\n" ); } catch ( e ) {}
}

function ksRemove( path )
{
   try { if ( File.exists( path ) ) File.remove( path ); } catch ( e ) {}
}

/*
 * Astrometrie als Stuetzpunkte: ein Gitter aus Bildkoordinaten und den zugehoerigen
 * Himmelskoordinaten. So ist das Python-Programm unabhaengig davon, wie PixInsight
 * die Loesung intern speichert (Splines, Verzeichnung).
 */
function ksAstrometryPoints( window )
{
   var img = window.mainView.image;
   var n = 16, lines = [];
   for ( var j = 0; j <= n; ++j )
      for ( var i = 0; i <= n; ++i )
      {
         var x = 0.5 + i*( img.width  - 1 )/n;
         var y = 0.5 + j*( img.height - 1 )/n;
         var q = null;
         try { q = window.imageToCelestial( new Point( x, y ) ); }
         catch ( e ) { try { q = window.imageToCelestial( x, y ); } catch ( e2 ) { q = null; } }
         if ( q !== null && q !== undefined && isFinite( q.x ) && isFinite( q.y ) )
            lines.push( format( "%.4f %.4f %.9f %.9f", x, y, q.x, q.y ) );
      }
   return lines;
}

function ksObservationDate( window )
{
   var kw = window.keywords;
   for ( var i = 0; i < kw.length; ++i )
      if ( kw[i].name == "DATE-OBS" )
      {
         var v = kw[i].strippedValue;
         if ( v && v.length >= 7 )
            return v.substring( 0, 10 );
      }
   try
   {
      var t = window.mainView.propertyValue( "Observation:Time:Start" );
      if ( t )
         return ( new Date( t ) ).toISOString().substring( 0, 10 );
   }
   catch ( e ) {}
   return "";
}

function ksSaveCopy( window, path )
{
   var src = window.mainView.image;
   var tmp = new ImageWindow( src.width, src.height, src.numberOfChannels, 32, true,
                              src.isColor, "ks_tmp" );
   tmp.mainView.beginProcess( UndoFlag_NoSwapFile );
   tmp.mainView.image.assign( src );
   tmp.mainView.endProcess();
   var ok = tmp.saveAs( path, false, false, false, false );
   tmp.forceClose();
   if ( !ok )
      throw new Error( "Temporaere Datei konnte nicht gespeichert werden: " + path );
}

function ksRunProcess( program, args )
{
   var P = new ExternalProcess;
   var out = "", err = "", buf = "";
   P.onStandardOutputDataAvailable = function()
   {
      var s = String( P.stdout );
      out += s;
      buf += s;
      var lines = buf.split( /\r?\n/ );
      buf = lines.pop();
      for ( var i = 0; i < lines.length; ++i )
         if ( lines[i].length > 0 )
         {
            console.writeln( lines[i] );
            ksTrace( lines[i] );
         }
      console.flush();
   };
   P.onStandardErrorDataAvailable = function()
   {
      err += String( P.stderr );
   };
   var started = false, finished = false, failed = false;
   P.onStarted  = function() { started = true; };
   P.onFinished = function() { finished = true; };
   P.onError    = function( code ) { if ( !started ) failed = true; };
   ksTrace( "Start: " + program + " " + args.join( " " ) );
   try { P.start( program, args ); }
   catch ( e ) { throw new Error( "Python konnte nicht gestartet werden:\n" + program + "\n" + e.message ); }
   while ( P.isStarting ) processEvents();
   processEvents();
   if ( failed || ( !started && !P.isRunning && !finished ) )
      throw new Error( "Python konnte nicht gestartet werden:\n" + program );
   while ( P.isRunning )
   {
      processEvents();
      if ( console.abortRequested )
      {
         P.kill();
         throw new Error( "Abgebrochen." );
      }
      msleep( 50 );
   }
   processEvents();
   if ( buf.length > 0 )
      console.writeln( buf );
   ksTrace( "Python beendet, Code " + P.exitCode + ( err.length > 0 ? "\nstderr: " + err : "" ) );
   return { exitCode: P.exitCode, stdout: out, stderr: err };
}

// ----------------------------------------------------------------------------
// Ausfuehrung
// ----------------------------------------------------------------------------

function ksExecute( window, P )
{
   if ( window === null || window === undefined || window.isNull )
      throw new Error( "Kein Bild ausgewaehlt." );
   if ( !File.exists( P.pythonPfad ) )
      throw new Error( "Python nicht gefunden:\n" + P.pythonPfad );
   if ( !File.exists( P.skriptPfad ) )
      throw new Error( "katalogschaerfung.py nicht gefunden:\n" + P.skriptPfad );

   var points = ksAstrometryPoints( window );
   if ( points.length < 30 )
      throw new Error( "Das Bild hat keine astrometrische Loesung.\nBitte zuerst ImageSolver ausfuehren." );

   var dir  = ksTempDir();
   var base = dir + "/ks_" + ( new Date ).getTime().toString();
   var fIn = base + "_in.xisf", fPts = base + "_wcs.txt", fOut = base + "_out.xisf", fMask = base + "_maske.xisf";
   var work = base + "_arbeit";

   console.show();
   console.abortEnabled = true;
   console.noteln( "<end><cbr><br>" + KS_TITLE + " " + KS_VERSION + ": " + window.mainView.id );
   console.flush();

   try
   {
      ksTraceLines = [];
      ksTrace( "Bild: " + window.mainView.id );
      ksSaveCopy( window, fIn );
      ksTrace( "Kopie gespeichert: " + fIn );
      File.writeTextFile( fPts, points.join( "\n" ) + "\n" );

      var args = [ "-X", "utf8", "-u", P.skriptPfad, fIn,
                   "--wcs-punkte", fPts, "--ausgabe", fOut, "--maske-ausgabe", fMask,
                   "--arbeit", work, "--katalog-cache", dir + "/kataloge",
                   "--stern-staerke", format( "%.2f", P.sternStaerke ),
                   "--nichtstern-staerke", format( "%.2f", P.nichtsternStaerke ),
                   "--nichtstern-psf", format( "%.2f", P.nichtsternPSF ),
                   "--galaxiensterne", P.galaxiensterne == 1 ? "alle" : "sicher" ];
      var date = ksObservationDate( window );
      if ( date.length > 0 ) { args.push( "--datum" ); args.push( date ); }
      if ( P.autoPSF )   args.push( "--auto-psf" );
      if ( P.korrektur ) args.push( "--korrektur" );

      var r = ksRunProcess( P.pythonPfad, args );
      if ( r.exitCode != 0 || !File.exists( fOut ) )
         throw new Error( "Die Katalogschaerfung ist fehlgeschlagen.\n\n" +
                          ( r.stderr.length > 0 ? r.stderr.substring( Math.max( 0, r.stderr.length - 1500 ) )
                                                : "Siehe Ausgabe in der Prozesskonsole." ) );

      ksTrace( "Ergebnis wird geladen" );
      var res = ImageWindow.open( fOut );
      if ( res.length < 1 )
         throw new Error( "Ergebnisdatei konnte nicht geoeffnet werden." );
      var rw = res[0];
      var target = window.mainView;
      if ( rw.mainView.image.width != target.image.width || rw.mainView.image.height != target.image.height )
      {
         rw.forceClose();
         throw new Error( "Ergebnis hat eine andere Groesse als das Bild." );
      }

      if ( P.neuesBild )
      {
         rw.mainView.id = target.id + "_katalog";
         rw.keywords = window.keywords;
         try { rw.copyAstrometricSolution( window ); } catch ( e ) {}
         rw.show();
      }
      else
      {
         target.beginProcess();
         target.image.assign( rw.mainView.image );
         target.endProcess();
         rw.forceClose();
      }

      if ( P.maskeOeffnen && File.exists( fMask ) )
      {
         var mw = ImageWindow.open( fMask );
         if ( mw.length > 0 )
         {
            mw[0].mainView.id = target.id + "_sternmaske";
            mw[0].show();
         }
      }
      console.noteln( KS_TITLE + ": fertig." );
      ksTrace( "fertig" );
   }
   catch ( e )
   {
      ksTrace( "FEHLER: " + e.message );
      throw e;
   }
   finally
   {
      console.abortEnabled = false;
      ksRemove( fIn ); ksRemove( fPts ); ksRemove( fOut ); ksRemove( fMask );
      try
      {
         if ( File.directoryExists( work ) )
         {
            var files = searchDirectory( work + "/*", false );
            for ( var i = 0; i < files.length; ++i )
               ksRemove( files[i] );
            File.removeDirectory( work );
         }
      }
      catch ( e ) {}
   }
}

// ----------------------------------------------------------------------------
// Dialog
// ----------------------------------------------------------------------------

function KSDialog()
{
   this.__base__ = Dialog;
   this.__base__();

   var self = this;
   var labelWidth = this.font.width( "Nicht-Stern-Staerke:" ) + 8;

   this.help = new Label( this );
   this.help.wordWrapping = true;
   this.help.useRichText = true;
   this.help.margin = 6;
   this.help.frameStyle = 1;
   this.help.text = "<b>" + KS_TITLE + " " + KS_VERSION + "</b><br/>" +
      "Cosmic Clarity mit Katalogpr&uuml;fung: Sterne werden nur dort gesch&auml;rft, wo Gaia DR3 oder der " +
      "Legacy Survey einen Stern kennen. Hintergrundgalaxien bekommen nur die Nicht-Stern-Sch&auml;rfung.<br/>" +
      "Das Bild muss <b>linear</b> sein und eine <b>astrometrische L&ouml;sung</b> haben. " +
      "Laufzeit: einige Minuten, beim ersten Mal pro Himmelsfeld zus&auml;tzlich der Katalog-Download.";

   this.viewLabel = new Label( this );
   this.viewLabel.text = "Bild:";
   this.viewLabel.minWidth = labelWidth;
   this.viewLabel.textAlignment = TextAlign_Right | TextAlign_VertCenter;
   this.viewList = new ViewList( this );
   this.viewList.getMainViews();
   if ( !ImageWindow.activeWindow.isNull )
      this.viewList.currentView = ImageWindow.activeWindow.mainView;
   var viewSizer = new HorizontalSizer;
   viewSizer.spacing = 6;
   viewSizer.add( this.viewLabel );
   viewSizer.add( this.viewList, 100 );

   function numeric( text, lo, hi, value, tip, setter )
   {
      var c = new NumericControl( self );
      c.label.text = text;
      c.label.minWidth = labelWidth;
      c.setRange( lo, hi );
      c.slider.setRange( 0, 200 );
      c.slider.minWidth = 200;
      c.setPrecision( 2 );
      c.setValue( value );
      c.toolTip = tip;
      c.onValueUpdated = setter;
      return c;
   }

   this.stern = numeric( "Stern-Staerke:", 0, 1, KSParams.sternStaerke,
      "Cosmic Clarity Stellar Amount", function( v ) { KSParams.sternStaerke = v; } );
   this.nichtstern = numeric( "Nicht-Stern-Staerke:", 0, 1, KSParams.nichtsternStaerke,
      "Cosmic Clarity Non-Stellar Amount", function( v ) { KSParams.nichtsternStaerke = v; } );
   this.psf = numeric( "Nicht-Stern-PSF:", 1, 8, KSParams.nichtsternPSF,
      "Cosmic Clarity Non-Stellar PSF (Radius in Pixel)", function( v ) { KSParams.nichtsternPSF = v; } );
   this.psf.enabled = !KSParams.autoPSF;

   function check( text, value, tip, setter )
   {
      var c = new CheckBox( self );
      c.text = text;
      c.checked = value;
      c.toolTip = tip;
      c.onCheck = setter;
      var s = new HorizontalSizer;
      s.addSpacing( labelWidth + 6 );
      s.add( c );
      s.addStretch();
      return s;
   }

   this.autoPSF = check( "PSF automatisch bestimmen", KSParams.autoPSF, "Cosmic Clarity Auto-PSF",
      function( on ) { KSParams.autoPSF = on; self.psf.enabled = !on; } );
   this.korrektur = check( "Sternkorrektur (verzogene Sterne)", KSParams.korrektur,
      "Schaltet die Aberrationskorrektur von Cosmic Clarity vor dem Schaerfen dazu.",
      function( on ) { KSParams.korrektur = on; } );
   this.maske = check( "Sternmaske als Bild oeffnen", KSParams.maskeOeffnen,
      "Zeigt, wo das Sternmodell uebernommen wurde.", function( on ) { KSParams.maskeOeffnen = on; } );
   this.neu = check( "Ergebnis als neues Bild", KSParams.neuesBild,
      "Aus: das Bild wird ersetzt (mit Rueckgaengig). An: das Ergebnis kommt in ein neues Fenster.",
      function( on ) { KSParams.neuesBild = on; } );

   this.galLabel = new Label( this );
   this.galLabel.text = "In grossen Galaxien:";
   this.galLabel.minWidth = labelWidth;
   this.galLabel.textAlignment = TextAlign_Right | TextAlign_VertCenter;
   this.galCombo = new ComboBox( this );
   this.galCombo.addItem( "nur sichere Vordergrundsterne schaerfen" );
   this.galCombo.addItem( "alle Punktquellen schaerfen (auch Knoten)" );
   this.galCombo.currentItem = KSParams.galaxiensterne;
   this.galCombo.toolTip = "Innerhalb grosser Galaxien wie M51 listen die Kataloge auch HII-Regionen und Sternhaufen als Punktquellen.";
   this.galCombo.onItemSelected = function( i ) { KSParams.galaxiensterne = i; };
   var galSizer = new HorizontalSizer;
   galSizer.spacing = 6;
   galSizer.add( this.galLabel );
   galSizer.add( this.galCombo, 100 );

   function pathRow( text, value, filter, setter )
   {
      var l = new Label( self );
      l.text = text;
      l.minWidth = labelWidth;
      l.textAlignment = TextAlign_Right | TextAlign_VertCenter;
      var e = new Edit( self );
      e.text = value;
      e.minWidth = 380;
      e.onEditCompleted = function() { setter( e.text.trim() ); };
      var b = new ToolButton( self );
      b.icon = self.scaledResource( ":/icons/select-file.png" );
      b.setScaledFixedSize( 22, 22 );
      b.onClick = function()
      {
         var d = new OpenFileDialog;
         d.filters = filter;
         if ( e.text.length > 0 )
            d.initialPath = e.text;
         if ( d.execute() )
         {
            e.text = d.fileName;
            setter( d.fileName );
         }
      };
      var s = new HorizontalSizer;
      s.spacing = 6;
      s.add( l );
      s.add( e, 100 );
      s.add( b );
      return s;
   }

   this.python = pathRow( "Python:", KSParams.pythonPfad, [ [ "Python", "python.exe", "*.exe" ] ],
      function( p ) { KSParams.pythonPfad = p; } );
   this.skript = pathRow( "Programm:", KSParams.skriptPfad, [ [ "Python-Skript", "*.py" ] ],
      function( p ) { KSParams.skriptPfad = p; } );

   this.ok = new PushButton( this );
   this.ok.text = "Ausfuehren";
   this.ok.icon = this.scaledResource( ":/icons/power.png" );
   this.ok.onClick = function() { self.ok(); };
   this.cancel = new PushButton( this );
   this.cancel.text = "Abbrechen";
   this.cancel.icon = this.scaledResource( ":/icons/cancel.png" );
   this.cancel.onClick = function() { self.cancel(); };
   var buttons = new HorizontalSizer;
   buttons.spacing = 6;
   buttons.addStretch();
   buttons.add( this.ok );
   buttons.add( this.cancel );

   this.sizer = new VerticalSizer;
   this.sizer.margin = 8;
   this.sizer.spacing = 6;
   this.sizer.add( this.help );
   this.sizer.addSpacing( 4 );
   this.sizer.add( viewSizer );
   this.sizer.add( this.stern );
   this.sizer.add( this.nichtstern );
   this.sizer.add( this.psf );
   this.sizer.add( this.autoPSF );
   this.sizer.add( this.korrektur );
   this.sizer.add( galSizer );
   this.sizer.add( this.maske );
   this.sizer.add( this.neu );
   this.sizer.addSpacing( 4 );
   this.sizer.add( this.python );
   this.sizer.add( this.skript );
   this.sizer.addSpacing( 4 );
   this.sizer.add( buttons );

   this.windowTitle = KS_TITLE;
   this.adjustToContents();
   this.setFixedHeight();
}

KSDialog.prototype = new Dialog;

function ksMain()
{
   KSParams.load();
   var dialog = new KSDialog;
   if ( !dialog.execute() )
      return;
   KSParams.save();
   var view = dialog.viewList.currentView;
   if ( view === null || view === undefined || view.isNull )
   {
      ( new MessageBox( "Bitte zuerst ein Bild oeffnen und auswaehlen.", KS_TITLE, StdIcon_Error, StdButton_Ok ) ).execute();
      return;
   }
   try
   {
      ksExecute( view.window, KSParams );
   }
   catch ( e )
   {
      console.criticalln( KS_TITLE + ": " + e.message );
      ( new MessageBox( e.message, KS_TITLE, StdIcon_Error, StdButton_Ok ) ).execute();
   }
}

#ifndef KS_NO_MAIN
ksMain();
#endif
