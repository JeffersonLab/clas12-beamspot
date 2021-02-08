import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import org.jlab.groot.math.*;
import org.jlab.groot.ui.TCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

// ==========================================
// Beam spot analysis
// 
// From the previous work by S. Stepanyan 
//  CLAS12 Note 2020-003
// 
// author: fbossu (at jlab.org);
// ==========================================


public class BeamSpot {

  // histograms
  // ----------------------------------------- 
  H1F h1_z;
  H1F h1_phi;
  ArrayList<H2F> a_h2_z_phi;  
  ArrayList<GraphErrors> a_g_results;
  ArrayList<Func1D> a_fits;

  // settings
  // -----------------------------------------
  boolean check_slices;
  double[] theta_bins;

  // 
  // ----------------------------------------- 
  public BeamSpot() {
    //h2_z_phi = new H2F("h2_z_phi","h2_phi_vs_z", 50, 15,40,36,-180,180);

    a_h2_z_phi = new ArrayList<H2F>();
    a_g_results = new ArrayList<GraphErrors>();
    a_fits = new ArrayList<Func1D>();
    
    h1_z   = new H1F( "h1_z",   "z vertex"  , 200, -20, 50 );
    h1_phi = new H1F( "h1_phi", "phi vertex", 200, 0, 360 );

    check_slices = true;
  }

  // initialize histograms
  // -----------------------------------------
  public void init() {
    for( int i = 0; i<theta_bins.length-1; i++ ){
      a_h2_z_phi.add( new H2F("h2_z_phi_"+i, "h2_phi_vs_z "+ theta_bins[i], 50, 15,40,36,0,360) );
      a_g_results.add( new GraphErrors() );
    }
  }  

  // setters
  // -----------------------------------------
  public void setCheckSlices( boolean t ) { check_slices = t; }

  public void setThetaBins( double[] bins ) { theta_bins = bins; }

  // 
  // ----------------------------------------- 
  public boolean processEvent( DataEvent event ){
    
    if( event.hasBank( "REC::Particle" ) == false ||
        event.hasBank( "REC::Track" )    == false    ) return false;

    DataBank bpart = event.getBank( "REC::Particle" );
    DataBank btrk  = event.getBank( "REC::Track" );

    // interfaces to the banks
    Track trk = new Track();
    Particle part = new Particle();

    // loop over tracks
    for( int i=0; i < btrk.rows(); i++){
      
      trk.load( btrk, i );

      // check the quality of the track
      if( checkTrack( trk ) == false ) continue;

      // get the corresponding particle
      part.load( bpart, trk.pindex );

      // check the quality of the particle
      if( checkParticle( part ) == false ) continue;

      // compute phi from the transvers momentum components
      float phi = (float) Math.toDegrees( Math.atan2( part.py, part.px) );
      if( phi < 0 ) phi += 360.0;
      float theta = (float) Math.toDegrees( Math.atan2( Math.sqrt( part.px*part.px + part.py*part.py), part.pz ) );

      // find theta bin
      int bin = Arrays.binarySearch( theta_bins, theta );
      bin = -bin -2;
      if( bin < 0 || bin >= theta_bins.length - 1 ) continue;

      // fill the histograms
      h1_z.fill( part.vz );
      h1_phi.fill( phi );

      a_h2_z_phi.get(bin).fill( part.vz, phi );
    } // end loop over tracks
    return true;
  }

  // track cuts
  // ------------------------------------
  boolean checkTrack ( Track trk ){
    // only use FD tracks
    if( trk.detector !=  org.jlab.detector.base.DetectorType.DC.getDetectorId() ) return false;

    // only negative tracks
    if( trk.q > 0 ) return false;

    // TODO additional cuts
    return true; 
  } 

  // particle cuts
  // ------------------------------------
  boolean checkParticle ( Particle p ){

    if ( p.pid != 11 ) return false;

    // TODO additional cuts
    return true; 
  } 

  // analysis
  // ------------------------------------
  public void analyze() {
    for( int i=0; i<theta_bins.length-1; i++ ){
      analyze( i );
    }
  }

  public void analyze( int i_theta_bin ) {

    GraphErrors g_results = a_g_results.get(i_theta_bin);
    H2F h2_z_phi = a_h2_z_phi.get(i_theta_bin);

    // loop over the phi bins of the 2D histogram phi vs z
    // and fit with a gaussian around the target window position
    double xmin = 18.;
    double xmax = 30.;

    // skip bins that are not populated by checking the maximum 
    // as reference we take the slice at phi=0
    H1F h0 = h2_z_phi.sliceY(h2_z_phi.getYAxis().getBin(0.));
    double max = h0.getBinContent( h0.getMaximumBin() );
    
    TCanvas c = null;
    int ic = 0;
    if( check_slices ){
      c = new TCanvas("cc"+i_theta_bin,900,900);
      c.divide( 6,7 );
    }
    // loop  over the y bins
    for( int i=0;i<h2_z_phi.getYAxis().getNBins(); i++ ){

      H1F h = h2_z_phi.sliceY( i );
      if( h.getBinContent( h.getMaximumBin() ) < 0.7*max ) continue;

      F1D func = new F1D( "func"+i, "[amp]*gaus(x,[mean],[sigma])", xmin, xmax ); 
      func.setParameter(0, h.getBinContent( h.getMaximumBin() ) );
      func.setParameter(1, h.getAxis().getBinCenter( h.getMaximumBin() )  ); 
      func.setParameter(2, 0.5 );
      DataFitter.fit( func, h, "Q" );
      g_results.addPoint( 
            h2_z_phi.getYAxis().getBinCenter( i ),
            func.getParameter(1),
            0,
            func.parameter(1).error() );
      if( check_slices ) {
        func.show();
        c.cd(ic++);
        c.draw(h);
        func.setLineColor( 2 );
        c.draw(func,"same");
        System.out.println( " ----- " + i);
      }
    } // end loop over bins

    if( check_slices ) c.save("fits"+i_theta_bin+".png");
   
    // fit the graph 
    FitFunc func = new FitFunc( "f1", 0., 360. );
    func.setParameter(0,25.0);
    func.setParameter(1,1.0);
    func.setParLimits(1,-0.1,10.);
    func.setParameter(2, Math.toRadians( 90.0 ) );
    func.setParLimits(2, -0.001, 2*Math.PI +0.001 );
    DataFitter.fit( func, g_results,"");
    func.setLineColor(2);
    func.setOptStat(11110);
    func.show();
    a_fits.add( func );
  }

  // plots
  // ------------------------------------
  public void plot() {

    TCanvas c = new TCanvas("c",800,600);
    c.divide(2,1);
    c.cd(0);
    c.draw( h1_z);
    c.cd(1);
    c.draw( h1_phi);

    for( int i=0; i<theta_bins.length-1; i++ ){

      TCanvas ci = new TCanvas( "c"+i, 800, 600 );
      ci.divide(2,1);
      ci.cd(0);
      ci.draw( a_h2_z_phi.get(i) );

      ci.cd(1);
      ci.draw( a_g_results.get(i) );

      ci.save( "bin"+i+".png");      
    }
  
    // plot the parameters
    double[] Z = new double[ theta_bins.length -1 ];
    double[] R = new double[ theta_bins.length -1 ];
    double[] P = new double[ theta_bins.length -1 ];
    double[] X = new double[ theta_bins.length -1 ];
    double[] Y = new double[ theta_bins.length -1 ];
    double[] T = new double[ theta_bins.length -1 ];
    double[] ET = new double[ theta_bins.length -1 ];

    double[] EZ = new double[ theta_bins.length -1 ];
    double[] ER = new double[ theta_bins.length -1 ];
    double[] EP = new double[ theta_bins.length -1 ];
    for( int i=0; i<theta_bins.length-1; i++ ){
      Func1D f = a_fits.get(i);
      
      // average theta
      T[i] = (theta_bins[i] + theta_bins[i+1])/2.;
      ET[i] = 0.0;

      Z[i] = f.getParameter( 0 );
      EZ[i] = f.parameter( 0 ).error();
      R[i] = f.getParameter(1) * Math.tan( Math.toRadians( T[i] ) );
      ER[i] = f.parameter(1).error() * Math.tan( Math.toRadians( T[i] ) );

      P[i] = Math.toDegrees( f.getParameter(2) );
      EP[i] = Math.toDegrees( f.parameter(2).error() );

      X[i] = R[i] * Math.cos( f.getParameter(2) ); 
      Y[i] = R[i] * Math.sin( f.getParameter(2) ); 

      System.out.println( P[i] + " " + f.getParameter(2) + " " + Math.toRadians( P[i] ));
 
    }
    GraphErrors gZ = new GraphErrors("Z", T, Z, ET, EZ );
    gZ.setTitleX( "#theta (degrees)" );
    gZ.setTitleY( "Z (cm)" );
    GraphErrors gR = new GraphErrors("R", T, R );
    gR.setTitleX( "#theta (degrees)" );
    gR.setTitleY( "R (cm)" );
    GraphErrors gP = new GraphErrors("Phi", T, P );
    gP.setTitleX( "#theta (degrees)" );
    gP.setTitleY( "#phi0 (degrees)" );
    GraphErrors gX = new GraphErrors("X", T, X );
    gX.setTitleX( "#theta (degrees)" );
    gX.setTitleY( "X (cm)" );
    GraphErrors gY = new GraphErrors("Y", T, Y );
    gY.setTitleX( "#theta (degrees)" );
    gY.setTitleY( "Y (cm)" );


    TCanvas cp = new TCanvas("cpars", 600,600 );
    cp.divide(2,3);
    cp.cd(0);
    cp.draw( gZ );    
    cp.cd(1);
    cp.draw( gR );    
    cp.cd(2);
    cp.draw( gP );    
    cp.cd(3);
    cp.draw( gX );    
    cp.cd(4);
    cp.draw( gY );    

    cp.save("results.png");
  }

  // fit function
  // -------------
  class FitFunc extends Func1D {
    public FitFunc( String name, double min, double max ) {
      super( name, min, max );
      
      addParameter( "a" );
      addParameter( "b" );
      addParameter( "c" );
    }

    @Override
    public double evaluate( double x ) {
      return this.getParameter(0) - this.getParameter(1) * Math.cos( x * Math.PI/180.0 - this.getParameter(2) );
    }
  }


  // bank interfaces
  // ---------------------------------------
  class Track {
    public short pindex; 
    public byte  detector;
    public byte  sector;
    public short status;
    public byte  q;
    public float chi2;
    public short NDF;

    public Track() { }
    public void load( DataBank bank, int i){
      pindex   = bank.getShort( "pindex",i);     
      detector = bank.getByte(  "detector",i);   
      sector   = bank.getByte(  "sector",i);     
      status   = bank.getShort( "status",i);     
      q        = bank.getByte(  "q",i);          
      chi2     = bank.getFloat( "chi2",i);       
      NDF      = bank.getShort( "NDF",i);        
    }
  }

  class Particle {
    public int pid;    
    public float px;     
    public float py;     
    public float pz;     
    public float vx;     
    public float vy;     
    public float vz;     
    public float vt;     
    public byte charge; 
    public float beta;   
    public float chi2pid;
    public short status; 

    public Particle() { }
    public void load( DataBank bank, int i){
      pid       = bank.getInt( "pid",i);
      px        = bank.getFloat( "px",i);
      py        = bank.getFloat( "py",i);
      pz        = bank.getFloat( "pz",i);
      vx        = bank.getFloat( "vx",i);
      vy        = bank.getFloat( "vy",i);
      vz        = bank.getFloat( "vz",i);
      vt        = bank.getFloat( "vt",i);
      charge    = bank.getByte( "charge",i);
      beta      = bank.getFloat( "beta",i);
      chi2pid   = bank.getFloat( "chi2pid",i);
      status    = bank.getShort( "status",i);
    }
  }


  // ------------------------------------
  // ##### MAIN ######
  // ------------------------------------
  public static void main( String[] args ){

    BeamSpot bs = new BeamSpot();

    // set the theta bin edges
    double[] bins = {5., 10., 20.,  35.};
    bs.setThetaBins( bins );

    // call the init method to properly setup all the parameters
    bs.init();

    // loop over input files
    for( String s : args ) {
      if ( s.endsWith( ".hipo" ) == false ) continue;

      HipoDataSource reader = new HipoDataSource();
      reader.open( s );
      int filecount = 0;
      while(reader.hasEvent() ) {
              DataEvent event = reader.getNextEvent();
        bs.processEvent( event );
        filecount++;
        //if( filecount > 1000 ) break;
      }// end loop on events
      reader.close();

    }// end loop on input files

    bs.analyze();    
    bs.plot();    
  }

}

