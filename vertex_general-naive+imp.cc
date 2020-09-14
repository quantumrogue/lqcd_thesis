#include "qdp.h"   // The core QDP++ library header
#include "fft4d.h"  // fft4d from fftw
#include "chroma.h"
#include <vector>  //adicionado para poder utilizar vetores
#include <list>

using namespace std; // Import from STD namespace (io etc)
using namespace QDP; // Import from QDP namespace (QDP++ things)
using namespace Chroma;
using namespace Landau;



/*
	- Compute GLUON PROPAGATOR using lattice tensor basis with 3 and 5 terms
	- Consider the change of phase to the gluon field

  - Using IMPROVED and NAIVE momenta for the basis
  - RECONSTRUCTION PROCEDURE
*/


//User-defined functions
Double flandau( const multi1d<LatticeColorMatrix>& );
void sd_delta( const multi1d<LatticeColorMatrix>& , LatticeColorMatrix&);

Double sd_theta( const multi1d<LatticeColorMatrix>& ); //given u
Double sd_theta( const LatticeColorMatrix& ); //given Delta

//Declaring vertex function
void gluon_prop(const multi1d<LatticeColorMatrix>&, multi1d<int>&);

void writeGauge_singleprec(XMLBufferWriter& file_xml,
    XMLBufferWriter& record_xml, 
    const multi1d<LatticeColorMatrix>& u, 
    const std::string& file, 
    QDP_volfmt_t volfmt, 
    QDP_serialparallel_t serpar);

// Here is our program
int main(int argc, char *argv[])
{

  Chroma::initialize(&argc, &argv);

  
  QDPIO::cout << "after initialize" << endl;

  // Put this in to enable profiling etc.
  START_CODE();

  QDPIO::cout << "after start code" << endl;

  //Read parameters from a input XML file

    multi1d<int> latinput(Nd);    
    int gf_iter;
    Double gf_prec;
    Double alpha;
    Double dF;

  XMLReader xml_input("gfix-landau.in.xml");
  read(xml_input,"/Params/nrow",latinput);
  read(xml_input,"/Params/gf_iter",gf_iter);
  read(xml_input,"/Params/gf_prec",gf_prec);
  read(xml_input,"/Params/alpha",alpha);
  read(xml_input,"/Params/dF",dF);

    QDPIO::cout << "## Read input parameters from input XML file:" << endl;
    QDPIO::cout << "##   Lattice size: " << latinput[0] << " " << latinput[1] 
       << " " << latinput[2] << " " << latinput[3] << endl;
    QDPIO::cout << "##   Max numb of iters:  " << gf_iter << endl;
    QDPIO::cout << "##   Gaugefix precision: 10^(-" << gf_prec << ")" << endl;
    QDPIO::cout << "##   Alpha: " << alpha  << endl;
    QDPIO::cout << "##   dF: " <<  dF << endl;


 QDPIO::cout << "after reading input file" << endl;

  // Setting up the lattice size: 4x4x4x8 for now
                                                                        

  multi1d<int> latt_size(Nd);
  latt_size[0] = latinput[0];
  latt_size[1] = latinput[1];
  latt_size[2] = latinput[2];
  latt_size[3] = latinput[3];

  Layout::setLattSize(latt_size);
  Layout::create();   // Setup the layout

  //Layout::logicalSize();

 QDPIO::cout << "layout ready" << endl;

  // Reading the configuration 
                                                                        
  XMLReader gauge_file_xml, gauge_xml;

  XMLBufferWriter xml_buf;
  push(xml_buf, "Cfg");  

  //Random input configuration
//  write(xml_buf, "cfg_type", "WEAK_FIELD");
//  write(xml_buf, "cfg_file", "dummy");

  //SZINQIO configuration
  write(xml_buf, "cfg_type", "SZINQIO");
  write(xml_buf, "cfg_file", "testconf.in");


  pop(xml_buf);

  XMLReader xml_in(xml_buf);
  Cfg_t   cfg;
  read(xml_in,"/Cfg", cfg);
                                                                      
  // The gauge field itself
  multi1d<LatticeColorMatrix> u(Nd);

  //readmilc_ps("testconf.in", u);
  //"testconf.in"

  
  QDP::StopWatch time;
  time.start();
  

  // A convenience routine to read various gauge field formats
  // gaugeStartup(gauge_file_xml, gauge_xml, u, cfg);
  readGauge(gauge_file_xml, gauge_xml, u, cfg.cfg_file, QDPIO_PARALLEL);
  


  time.stop();
  QDPIO::cout << "Time spent in read config file:  " << time.getTimeInSeconds()  << " s"<< endl;
  time.reset();



  // Check if the gauge field configuration is unitarized
   unitarityCheck(u);
  // ---------------- Gauge Reading Done -------------------------------


   // Definition of invpsq
                                                                                                    
   Double psqmin=0.00001, psqmax=Nd;
   Double PI=3.14159265358979323846264338;
   double PI_double = 3.14159265358979323846264338;

   //QDPIO::cout << "senode 0 " << sin(0) << endl;
   //QDPIO::cout << "senode PI/2 " << sin(PI/2) << endl;
   //QDPIO::cout << "senode PI " << sin(PI) << endl;
 

   LatticeReal invpsq;                                                                              

   for (int ix=0; ix < latt_size[0]; ix++){
     Double sx=sin(ix*PI/latt_size[0]);
   for (int iy=0; iy < latt_size[1]; iy++){
     Double sy=sin(iy*PI/latt_size[1]);
   for (int iz=0; iz < latt_size[2]; iz++){
     Double sz=sin(iz*PI/latt_size[2]);
   for (int it=0; it < latt_size[3]; it++){
     Double st=sin(it*PI/latt_size[3]);

     Double sinsq=sx*sx+sy*sy+sz*sz+st*st;

     Double prcfact=psqmax/sinsq;

     multi1d<int> posi(Nd);
     posi[0]=ix; posi[1]=iy; posi[2]=iz; posi[3]=it;

     if ( toFloat(sinsq) > toFloat(psqmin) ) {                                                    
          pokeSite(invpsq,prcfact,posi ); }
     else {
        pokeSite(invpsq,Double(0),posi);
     }
     

   }}}}

   //Double alpha = 0.06;

   // Gauge transformation matrix
   LatticeColorMatrix g;
   g=Double(1);

  // F function
   Double F_Landau = flandau(u);
   QDPIO::cout << "F = " << F_Landau << endl;

  // Delta
   LatticeColorMatrix Delta;
    sd_delta(u, Delta);

    // Theta
    Double Theta = sd_theta(Delta);

    QDPIO::cout << "Theta = " << Theta << endl;


    for (int iter=1; iter < gf_iter; iter++) {

     // Delta
     // sd_delta(u, Delta);

     //FFT acceleration
     int fft_for=-1, fft_bck=1;
     LatticeColorMatrix daux;
     daux=Delta;

     fft4d(daux,latt_size,fft_for);
     Delta=invpsq*daux;
     fft4d(Delta,latt_size,fft_bck); 

     Delta*=(Double(1.0)/Layout::vol());

     // g  (define the gauge transformation)
     g = Double(1);
     g += Double(alpha/2) * Delta; //+ Double(0.00125/2)*Delta*Delta;

     Chroma::reunit( g );

     // gauge transformation
     //isto corresponde a U' = g(x)Ug(x+mu) -> transf de gauge do link
     for(int mu = 0; mu < Nd; mu++) {
       LatticeColorMatrix u_tmp = g * u[mu];
       u[mu] = u_tmp * shift(adj(g), FORWARD, mu);                                
     }

     // F function
     Double F_landau_old = F_Landau;
     F_Landau = flandau(u);

     // Delta
     sd_delta(u, Delta);
     // Theta
     Theta = sd_theta(Delta);
     
     //Output the results
     QDPIO::cout << " Iteration nr. " << iter  <<  ";  F, Theta = " << F_Landau << " , " << Theta << endl;

     //Changing alpha if F is decreasing
     if ( (toFloat(F_Landau-F_landau_old) < -toFloat(dF))               //df definido no xml file
                  & (toFloat(alpha)>0.02) ) {
        alpha=0.95*alpha;
        QDPIO::cout << "    Warning: changing alpha down " << alpha << endl;
     } 

     
     //Reading new alpha  and dF if needed
     Double alphanew, dFnew;      
     if (iter % 10 == 0   ) {
      
       XMLReader xml_input("gfix-landau.in.xml");
       read(xml_input,"/Params/alphanew",alphanew); 
       read(xml_input,"/Params/dFnew",dFnew);   

       if ( toFloat(alphanew) > 0.0 ) { 
         QDPIO::cout << "    Warning: reading new alpha " << alphanew << endl;
         alpha=alphanew;
       }
       if ( toFloat(dFnew) > 0.0) {
         QDPIO::cout << "    Warning: reading new dF " << dFnew << endl;
         dF=dFnew;
       }
    
    }

     //Stopping if achieved required precision
     if ( toFloat(log10(Theta)) < toFloat(-gf_prec) ) break;
    }

  //ATÉ AQUI VERIFICOU-SE A GAUGE DE LANDAU

    //NOT writing the config
    //Computing gluon three vertex
    time.start();
    gluon_prop(u, latt_size);
    time.stop();

    QDPIO::cout << "Time spent in computing gluon propagator:  " << time.getTimeInSeconds()  << " s"<< endl;
    time.reset();


  // Clean up QDP++
  QDP_finalize();
  exit(0);  // Normal exit
}

// ###########################################################################################################################

Double flandau( const multi1d<LatticeColorMatrix>& u ){
   // F function
   //Double norm_F = Double(1)/Double(Layout::vol()*Nd*Nc) ; 
   Double norm_F = Double(1)/(Double(Layout::vol())*Double(Nd*Nc)) ; 
   Double F_Landau = Double(0);

   for(int mu=0; mu < Nd; mu++) {
     F_Landau += norm_F*sum(real(trace(u[mu])));
   }
 return F_Landau;
 }


   void sd_delta( const multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& Delta){
     // Delta
     Delta = Double(0);

     for(int nu=0; nu < Nd; nu++) { 
       LatticeColorMatrix u_nu_x_minus3_nu;
       u_nu_x_minus3_nu = shift(u[nu], BACKWARD, nu);
       LatticeColorMatrix tmp = u_nu_x_minus3_nu - u[nu];
       Chroma::taproj(tmp);                                                   
       Delta+=Double(2.0)*tmp;
     }

   }


Double sd_theta( const multi1d<LatticeColorMatrix>& u ){

  // Delta
   LatticeColorMatrix Delta;
    sd_delta(u, Delta);

    // Theta
    Double Norm_th = Double(1)/(Double(Layout::vol())*Double(Nc));
  LatticeColorMatrix tmp3= Delta*adj(Delta);
  Double Theta = Norm_th*sum(real(trace(tmp3)));

  return Theta;
}

Double sd_theta( const LatticeColorMatrix& Delta ){                                   

    // Theta
  Double Norm_th = Double(1)/(Double(Layout::vol())*Double(Nc));
  LatticeColorMatrix tmp3= Delta*adj(Delta);
  Double Theta = Norm_th*sum(real(trace(tmp3)));

  return Theta;
}

//############################################################################################################################ COMPUTE GLUON PROP

void gluon_prop(const multi1d<LatticeColorMatrix>& u, multi1d<int>& latt_size){

  multi1d<LatticeColorMatrix> A_x(Nd);
  Complex imagunit=cmplx(Double(0.0), Double(1.0));
  Double PI=3.14159265358979323846264338;
  double PI_double=3.14159265358979323846264338;
 
    //GLUON FIELD
  for(int nu=0; nu < Nd; nu++) { 
    A_x[nu]=u[nu];
    Chroma::taproj(A_x[nu]);
    //A_x[nu]*=Double(2.0);

    //Division by i
    QDPIO::cout << "Division of A_x[nu] by i: " << nu << endl;
    A_x[nu]*=Double(1)/imagunit;
    }

   multi1d<LatticeInt> indmom(Nd);

   multi1d<LatticeDouble> naivemom(Nd), tlimpmom(Nd);
 
  for(int mu=0; mu < Nd; mu++) {
    
    indmom[mu]=Layout::latticeCoordinate(mu);

    for (int nm=latt_size[mu]/2+1; nm<latt_size[mu]; nm++){
      Set timeslice;
      timeslice.make(TimeSliceFunc(mu));
      Subset tnm = timeslice[nm];
      indmom[mu][tnm]-=latt_size[mu];
    }

    naivemom[mu]=Double(2.0)*PI*LatticeDouble(indmom[mu])/LatticeDouble(latt_size[mu]);
    tlimpmom[mu]=Double(2.0)*sin(PI*LatticeDouble(indmom[mu])/LatticeDouble(latt_size[mu]));
  }

  multi1d<LatticeColorMatrix> A_p(Nd);

  // FFT
  int fft_dir=-1;  //-1 -> Forward; 1 -> Backward

  for(int mu=0; mu < Nd; mu++) {
    
    QDPIO::cout << ">> FFT mu= " << mu << endl;
    fft4d(A_x[mu],latt_size,fft_dir);
    A_p[mu]=A_x[mu];

    // Additional factor exp(-iq\mu/2)
    //LatticeInt ind_mu=Layout::latticeCoordinate(mu);
    //LatticeDouble angle=-PI*LatticeDouble(ind_mu)/LatticeDouble(latt_size[mu]);
    LatticeDouble angle= Double(-0.5)*naivemom[mu];
    LatticeComplex factexp=cmplx(cos(angle), sin(angle));
    A_p[mu]*=factexp;
  }


  //Building the LatticeColorMatrix for the zero momentum point 
  //create latpoint for the zero momentum position
  multi1d<int> zero_mom_latpoint(Nd);
  for (int i = 0; i < Nd; i++){
    zero_mom_latpoint[i] = 0;
  } 

  //RANDOM DEFINITIONS 
  double radian_converter = PI_double/180.0;
  double twopi_d = 2.0*PI_double;
  Double twopi_D = Double(2)*PI;
  Double group_permut = Double(24)*Double(16);
  Double lattice_volume = Double(Layout::vol());
  Double color_factor = Nc*(Nc*Nc-1)/Double(4);
  Double permut_plus_color = group_permut*color_factor;
  multi1d<int> shifted(Nd);

  int a = 0;
  int b = 1;
  int c = 2;
  int d = 3;
  //Z4 AVERAGING PERMUTATIONS
  int z4_ave[24][4]={{a,b,c,d},{a,b,d,c},{a,c,b,d},{a,c,d,b},
    {a,d,b,c},{a,d,c,b},{b,a,c,d},{b,a,d,c},
    {b,c,a,d},{b,c,d,a},{b,d,a,c},{b,d,c,a},
    {c,a,b,d},{c,a,d,b},{c,b,a,d},{c,b,d,a},
    {c,d,a,b},{c,d,b,a},{d,a,b,c},{d,a,c,b},
    {d,b,a,c},{d,b,c,a},{d,c,a,b},{d,c,b,a}};

  //NEGATIVE PERMUTATIONS
  int neg_permuts[16][4] = {{1,1,1,1},
      {-1,1,1,1},{1,-1,1,1},{1,1,-1,1},{1,1,1,-1},
      {-1,-1,1,1},{-1,1,-1,1},{-1,1,1,-1},{1,-1,-1,1},{1,-1,1,-1},{1,1,-1,-1},
      {-1,-1,-1,1},{-1,-1,1,-1},{-1,1,-1,-1},{1,-1,-1,-1},
      {-1,-1,-1,-1}};

  //Identity matrix
  multi2d<int> identity(Nd,Nd);
  for (int mu = 0; mu < Nd; mu++){
    for (int nu = 0; nu < Nd; nu++){
      if (mu == nu){
        identity[mu][nu] = 1;
      }
      else{
        identity[mu][nu] = 0; 
      }
    }
  }
  //--------------------------------------------------------------------------------------------------------------------
  //TENSORS
  //Create gluon propagator
  Double gdi;
  //intermediate tensors - no z4
  Double Ei,Fi,Gi,Hi,Ii;
  Double Ji, Ki, Li;
  Double Ai,Bi;
  //averaged tensors
  Double E,F,G,H,I;
  Double J,K,L;
  Double A,B;
  //imp
  Double im_Ei,im_Fi,im_Gi,im_Hi,im_Ii;
  Double im_Ji, im_Ki, im_Li;
  Double im_Ai,im_Bi;
  //averaged tensors
  Double im_E,im_F,im_G,im_H,im_I;
  Double im_J,im_K,im_L;
  Double im_A,im_B;

  //H4 invariants and momenta
  multi1d<Double> pnaive(Nd), pimp(Nd);
  Double p2, p4, p6, p8, p10;
  Double im_p2, im_p4, im_p6, im_p8, im_p10;

  //MATRICES
  ColorMatrix matrixA1, matrixA2;
  ColorMatrix matrixA1_adj;

  multi1d<int> negpoint(Nd), latpoint(Nd);
  //---------------------------------------------------------------------------------------------------------------------

  
  //GLUON PROPAGATOR with z4 average                                      
  for (int ipx=0; ipx < latt_size[0]/2 + 1; ipx++){                      
  for (int ipy=0; ipy < ipx + 1; ipy++){                                    
  for (int ipz=0; ipz < ipy + 1; ipz++){                                    
  for (int ipt=0; ipt < ipz + 1; ipt++){                                    
  //estes limites são considerados para apenas passar uma unica vez por cada ponto da rede dos momentos
    
    //RECONSTRUCTION
    //original, lattice reconstructed, lattice 3 forms, continuum reconstructed, continuum reconstr - no orthog
    Double D_orig = zero; 
    Double D_recons = zero, D_recons2 = zero, D_recons_cont = zero, D_recons_cont2 = zero;
    Double im_D_recons = zero, im_D_recons2 = zero, im_D_recons_cont = zero, im_D_recons_cont2 = zero;

    Double gd = zero;
    //lattice basis
    E = zero; F = zero; G = zero; H = zero; I = zero;
    J = zero; K = zero; L = zero;
    //continuum basis - no orthog
    A = zero; B = zero;
    //imp
    //lattice basis
    J = zero; K = zero; L = zero;
    im_E = zero; im_F = zero; im_G = zero; im_H = zero; im_I = zero;
    //continuum basis - no orthog
    im_A = zero; im_B = zero;

    //reset invariants
    p2 = zero; p4 = zero; p6 = zero; p8 = zero; p10 = zero;
    im_p2 = zero; im_p4 = zero; im_p6 = zero; im_p8 = zero; im_p10 = zero;
    Double im_p2 = zero;

    multi1d<int> mpoint(Nd);
    mpoint[0] = ipx; mpoint[1] = ipy; mpoint[2] = ipz; mpoint[3] = ipt;
    multi1d<int> minus_mpoint(Nd), minus_negpoint(Nd), minus_latpoint(Nd), minus_shifted(Nd); //momentum -p

    for(int mu = 0; mu < Nd; mu++){
    	minus_mpoint[mu] = -mpoint[mu];
      //momenta
      pnaive[mu] = Double(2)*PI*Double(mpoint[mu])/Double(latt_size[mu]);
      pimp[mu] = Double(2)*sin(mpoint[mu]*PI/latt_size[mu]);
      //invariants 
      Double pmu2 = pnaive[mu]*pnaive[mu]; 
      p2 += pmu2;
      p4 += pmu2*pmu2;
      p6 += pmu2*pmu2*pmu2;
      p8 += pmu2*pmu2*pmu2*pmu2;
      p10 += pmu2*pmu2*pmu2*pmu2*pmu2;
      Double im_pmu2 = pimp[mu]*pimp[mu]; 
      im_p2 += im_pmu2;
      im_p4 += im_pmu2*im_pmu2;
      im_p6 += im_pmu2*im_pmu2*im_pmu2;
      im_p8 += im_pmu2*im_pmu2*im_pmu2*im_pmu2;
      im_p10 += im_pmu2*im_pmu2*im_pmu2*im_pmu2*im_pmu2;
    }

    Double delta1, delta2;
    delta1 = Double(Nd)*(p4*p8 - p6*p6) + p2*(p4*p6 - p2*p8) + p4*(p2*p6 - p4*p4);
    delta2 = Double(2)*(p2*p4 - p6)*(p8 - p4*p4) + Double(2)*(p2*p2 - p4)*(p4*p6 - p10);
    //imp
    Double im_delta1, im_delta2;
    im_delta1 = Double(Nd)*(im_p4*im_p8 - im_p6*im_p6) + im_p2*(im_p4*im_p6 - im_p2*im_p8) + im_p4*(im_p2*im_p6 - im_p4*im_p4);
    im_delta2 = Double(2)*(im_p2*im_p4 - im_p6)*(im_p8 - im_p4*im_p4) + Double(2)*(im_p2*im_p2 - im_p4)*(im_p4*im_p6 - im_p10);

    QDPIO::cerr << p2 << "  " << im_delta1 << " " << im_delta2<< "  " << delta1 << " " << delta2<< endl;

    //ignore if redundant basis
    double zerocond = 0.0000001;
    if(toDouble(delta1) < zerocond || toDouble(delta2) < zerocond){
      break;
    }

    //COMPLETE Z4 AVERAGING
    for(int neg = 0; neg < 16; neg++){

      //perform negative permutations
      for(int mu = 0; mu < Nd; mu++){
        negpoint[mu] = mpoint[mu]*neg_permuts[neg][mu];
        minus_negpoint[mu] = minus_mpoint[mu]*neg_permuts[neg][mu];
      }

      for (int a4=0; a4 < 24; a4++){ 

        for(int mu = 0; mu < Nd; mu++){
          latpoint[mu] = negpoint[z4_ave[a4][mu]];
          pnaive[mu] = Double(2)*latpoint[mu]*PI/Double(latt_size[mu]);
          pimp[mu] = Double(2)*sin(latpoint[mu]*PI/latt_size[mu]); 
          minus_latpoint[mu] = minus_negpoint[z4_ave[a4][mu]];
        }

        //gluon propagator - continuum
        gdi = zero;
        //lattice basis
        Ei = zero; Fi = zero; Gi = zero; Hi = zero; Ii = zero;
        Ji = zero; Ki = zero; Li = zero;
        //continuum basis - no orthog
        Ai = zero, Bi = zero;
        //imp
        im_Ei = zero; im_Fi = zero; im_Gi = zero; im_Hi = zero; im_Ii = zero;
        im_Ji = zero; im_Ki = zero; im_Li = zero;
        //continuum basis - no orthog
        im_Ai = zero, im_Bi = zero;

        //assign k1 to change sign if necessary
        bool k1[Nd], k_1[Nd];
        for(int mu = 0; mu < Nd; mu++){
          k1[mu] = false;
          if(latpoint[mu] == -latt_size[mu]/2){
            k1[mu] = true;
          }
          k_1[mu] = false;
          if(minus_latpoint[mu] == -latt_size[mu]/2){
            k_1[mu] = true;
          }
        }

        //create shifted momenta
        for (int mu=0; mu < Nd; mu++){
            if(latpoint[mu] < 0){
              shifted[mu] = latpoint[mu] + latt_size[mu];
            }
            else{
              shifted[mu] = latpoint[mu];
            }
            if(minus_latpoint[mu] < 0){
              minus_shifted[mu] = minus_latpoint[mu] + latt_size[mu];
            }
            else{
              minus_shifted[mu] = minus_latpoint[mu];
            }
        }

        Double realAAmumu, realAAmunu, imagAAmumu, imagAAmunu;
        Double Dmumu = zero, p2Dmumu = zero, p4Dmumu = zero;
        Double ppDmunu = zero, p3p3Dmunu = zero;
        Double ppDmunu_complete = zero;
        //imp
        Double im_p2Dmumu = zero, im_p4Dmumu = zero;
        Double im_ppDmunu = zero, im_p3p3Dmunu = zero;
        Double im_ppDmunu_complete = zero;


        //EXTRACTING MATRICES
        for(int mu1 = 0; mu1 < Nd; mu1++){

          //matrices
        	matrixA1 = peekSite(A_p[mu1], shifted);
        	matrixA1_adj = adj(matrixA1);
          realAAmumu = real(trace(matrixA1*matrixA1_adj));
          imagAAmumu = imag(trace(matrixA1*matrixA1_adj));
          matrixA1 *= (k1[mu1] ? -1:1);

          //tensors
        	Dmumu +=  realAAmumu;
          p2Dmumu += pnaive[mu1]*pnaive[mu1]*realAAmumu;
          p4Dmumu += pnaive[mu1]*pnaive[mu1]*pnaive[mu1]*pnaive[mu1]*realAAmumu;
          //imp
          im_p2Dmumu += pimp[mu1]*pimp[mu1]*realAAmumu;
          im_p4Dmumu += pimp[mu1]*pimp[mu1]*pimp[mu1]*pimp[mu1]*realAAmumu;

          for(int mu2 = 0; mu2 < Nd; mu2++){

            matrixA2 = peekSite(A_p[mu2], minus_shifted);
            matrixA2 *= (k_1[mu2] ? -1:1);
            realAAmunu = real(trace(matrixA1*matrixA2));
            imagAAmunu = imag(trace(matrixA1*matrixA2));
            D_orig += sqrt(realAAmunu*realAAmunu);
            if(mu1 != mu2){
              ppDmunu += pnaive[mu1]*pnaive[mu2]*realAAmunu;
              p3p3Dmunu += pnaive[mu1]*pnaive[mu1]*pnaive[mu1]*pnaive[mu2]*pnaive[mu2]*pnaive[mu2]*realAAmunu;
              im_ppDmunu += pimp[mu1]*pimp[mu2]*realAAmunu;
              im_p3p3Dmunu += pimp[mu1]*pimp[mu1]*pimp[mu1]*pimp[mu2]*pimp[mu2]*pimp[mu2]*realAAmunu;
            }
            ppDmunu_complete += pnaive[mu1]*pnaive[mu2]*realAAmunu;
            im_ppDmunu_complete += pimp[mu1]*pimp[mu2]*realAAmunu;

          }//mu2
        }//mu1

        //Build form factors
        //assuming continuum form - orthog 
        gd += Dmumu/Double(Nd - 1); 

        //assuming continuum basis with orthog - form factor to reconstruct
        gdi = Dmumu/Double(Nd - 1);

        //LATTICE
        //Diagonal
        Ei = (Dmumu*(p4*p8 - p6*p6) + p2Dmumu*(p4*p6 - p2*p8) + p4Dmumu*(p2*p6 - p4*p4))/delta1;
        Fi = (Dmumu*(p4*p6 - p2*p8) + p2Dmumu*(Double(Nd)*p8 - p4*p4) + p4Dmumu*(p2*p4 - Double(Nd)*p6))/delta1;
        Gi = (Dmumu*(p2*p6 - p4*p4) + p2Dmumu*(p2*p4 - Double(Nd)*p6) + p4Dmumu*(Double(Nd)*p4 - p2*p2))/delta1;
        E += Ei; F += Fi; G += Gi;
        Ji = (Dmumu*p4 - p2Dmumu*p2)/(Double(Nd)*p4 - p2*p2);
        Ki = (-p2*Dmumu + Double(Nd)*p2Dmumu)/(Double(Nd)*p4 - p2*p2);
        J += Ji; K += Ki;
        //Off-Diagonal
        Hi = Double(2)*((ppDmunu)*(p4*p6 - p10) - (p3p3Dmunu)*(p2*p4 - p6))/delta2;
        Ii = ((ppDmunu)*(p8 - p4*p4) + (p3p3Dmunu)*(p2*p2 - p4))/delta2;
        H += Hi; I += Ii;
        Li = (ppDmunu)/(p2*p2 - p4);
        L += Li;
        //continuum - no orthog
        Ai = (Dmumu - ppDmunu_complete/p2)/Double(Nd-1);
        Bi = (-Dmumu/p2 + Double(Nd)*ppDmunu_complete/(p2*p2))/Double(Nd-1);
        A += Ai; B += Bi;
        //imp
        im_Ei = (Dmumu*(im_p4*im_p8 - im_p6*im_p6) + im_p2Dmumu*(im_p4*im_p6 - im_p2*im_p8) + im_p4Dmumu*(im_p2*im_p6 - im_p4*im_p4))/im_delta1;
        im_Fi = (Dmumu*(im_p4*im_p6 - im_p2*im_p8) + im_p2Dmumu*(Double(Nd)*im_p8 - im_p4*im_p4) + im_p4Dmumu*(im_p2*im_p4 - Double(Nd)*im_p6))/im_delta1;
        im_Gi = (Dmumu*(im_p2*im_p6 - im_p4*im_p4) + im_p2Dmumu*(im_p2*im_p4 - Double(Nd)*im_p6) + im_p4Dmumu*(Double(Nd)*im_p4 - im_p2*im_p2))/im_delta1;
        im_E += im_Ei; im_F += im_Fi; im_G += im_Gi;
        im_Ji = (Dmumu*im_p4 - im_p2Dmumu*im_p2)/(Double(Nd)*im_p4 - im_p2*im_p2);
        im_Ki = (-im_p2*Dmumu + Double(Nd)*im_p2Dmumu)/(Double(Nd)*im_p4 - im_p2*im_p2);
        im_J += im_Ji; im_K += im_Ki;
        //Off-Diagonal
        im_Hi = Double(2)*((im_ppDmunu)*(im_p4*im_p6 - im_p10) - (im_p3p3Dmunu)*(im_p2*im_p4 - im_p6))/im_delta2;
        im_Ii = ((im_ppDmunu)*(im_p8 - im_p4*im_p4) + (im_p3p3Dmunu)*(im_p2*im_p2 - im_p4))/im_delta2;
        im_H += im_Hi; im_I += im_Ii;
        im_Li = (im_ppDmunu)/(im_p2*im_p2 - im_p4);
        im_L += im_Li;
        //continuum - no orthog
        im_Ai = (Dmumu - im_ppDmunu_complete/im_p2)/Double(Nd-1);
        im_Bi = (-Dmumu/im_p2 + Double(Nd)*im_ppDmunu_complete/(im_p2*im_p2))/Double(Nd-1);
        im_A += im_Ai; im_B += im_Bi;

        //reconstruction
        //tensor_recons:  using lattice basis - 5 form factors
        //tensor_recons2: using lattice basis - 3 form factors
        //tensor_recons_cont: using continuum basis assuming orthogonaliny
        //tensor_recons_cont2: using continuum basis without orthog
        multi2d<Double> tensor_recons(Nd,Nd), tensor_recons2(Nd,Nd), tensor_recons_cont(Nd,Nd), tensor_recons_cont2(Nd,Nd);
        multi2d<Double> im_tensor_recons(Nd,Nd), im_tensor_recons2(Nd,Nd), im_tensor_recons_cont(Nd,Nd), im_tensor_recons_cont2(Nd,Nd);
        //zero valued tensors
        for(int mu1 = 0; mu1 < Nd; mu1++){
          for(int mu2 = 0; mu2 < Nd; mu2++){
            tensor_recons[mu1][mu2] = zero;
            tensor_recons2[mu1][mu2] = zero;
            tensor_recons_cont[mu1][mu2] = zero;
            tensor_recons_cont2[mu1][mu2] = zero;
            //imp
            im_tensor_recons[mu1][mu2] = zero;
            im_tensor_recons2[mu1][mu2] = zero;
            im_tensor_recons_cont[mu1][mu2] = zero;
            im_tensor_recons_cont2[mu1][mu2] = zero;
          }
        }

        for(int mu1 = 0; mu1 < Nd; mu1++){

          tensor_recons[mu1][mu1] += Ei + pnaive[mu1]*pnaive[mu1]*Fi + pnaive[mu1]*pnaive[mu1]*pnaive[mu1]*pnaive[mu1]*Gi;
          tensor_recons2[mu1][mu1] += Ji + pnaive[mu1]*pnaive[mu1]*Ki;
          tensor_recons_cont[mu1][mu1] += gdi;
          tensor_recons_cont2[mu1][mu1] += Ai;
          //imp
          im_tensor_recons[mu1][mu1] += im_Ei + pimp[mu1]*pimp[mu1]*im_Fi + pimp[mu1]*pimp[mu1]*pimp[mu1]*pimp[mu1]*im_Gi;
          im_tensor_recons2[mu1][mu1] += im_Ji + pimp[mu1]*pimp[mu1]*im_Ki;
          im_tensor_recons_cont[mu1][mu1] += gdi;
          im_tensor_recons_cont2[mu1][mu1] += im_Ai;

          for(int mu2 = 0; mu2 < Nd; mu2++){

            if(mu1 != mu2){
              tensor_recons[mu1][mu2] += pnaive[mu1]*pnaive[mu2]*Hi + pnaive[mu1]*pnaive[mu2]*(pnaive[mu1]*pnaive[mu1] + pnaive[mu2]*pnaive[mu2])*Ii;
              tensor_recons2[mu1][mu2] += pnaive[mu1]*pnaive[mu2]*Li;
              //imp
              im_tensor_recons[mu1][mu2] += pimp[mu1]*pimp[mu2]*im_Hi + pimp[mu1]*pimp[mu2]*(pimp[mu1]*pimp[mu1] + pimp[mu2]*pimp[mu2])*im_Ii;
              im_tensor_recons2[mu1][mu2] += pimp[mu1]*pimp[mu2]*im_Li;
            }
            tensor_recons_cont[mu1][mu2] += -(pnaive[mu1]*pnaive[mu2]/p2)*gdi;
            tensor_recons_cont2[mu1][mu2] += Bi*pnaive[mu1]*pnaive[mu2];
            //imp
            im_tensor_recons_cont[mu1][mu2] += -(pimp[mu1]*pimp[mu2]/im_p2)*gdi;
            im_tensor_recons_cont2[mu1][mu2] += im_Bi*pimp[mu1]*pimp[mu2];
          }//mu2 recons tensor
        }// mu1 recons tensor

        //computed reconstructed averaged tensors
        for(int mu1 = 0; mu1 < Nd; mu1++){
          for(int mu2 = 0; mu2 < Nd; mu2++){
            D_recons += sqrt(tensor_recons[mu1][mu2]*tensor_recons[mu1][mu2]);
            D_recons2 += sqrt(tensor_recons2[mu1][mu2]*tensor_recons2[mu1][mu2]);
            D_recons_cont += sqrt(tensor_recons_cont[mu1][mu2]*tensor_recons_cont[mu1][mu2]);
            D_recons_cont2 += sqrt(tensor_recons_cont2[mu1][mu2]*tensor_recons_cont2[mu1][mu2]);
            //imp
            im_D_recons += sqrt(im_tensor_recons[mu1][mu2]*im_tensor_recons[mu1][mu2]);
            im_D_recons2 += sqrt(im_tensor_recons2[mu1][mu2]*im_tensor_recons2[mu1][mu2]);
            im_D_recons_cont += sqrt(im_tensor_recons_cont[mu1][mu2]*im_tensor_recons_cont[mu1][mu2]);
            im_D_recons_cont2 += sqrt(im_tensor_recons_cont2[mu1][mu2]*im_tensor_recons_cont2[mu1][mu2]);
          }//mu2 recons sum
        }// mu1 recons sum

      }//Z4 AVERAGING
    }//neg permutations

    //Divide by 24 to account for all 24 permutations of momentum and 16 for all negative permutations performed
    //Divide by normalization factors:
    gd *= Double(2)/(Double(group_permut*(Nc*Nc-1)*lattice_volume));    
    E *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    F *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    G *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    H *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    I *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    J *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    K *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    L *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    A *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    B *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    D_orig *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    D_recons *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    D_recons2 *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    D_recons_cont *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    D_recons_cont2 *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    //imp
    im_E *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_F *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_G *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_H *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_I *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_J *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_K *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_L *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_A *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_B *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_D_recons *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_D_recons2 *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_D_recons_cont *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);
    im_D_recons_cont2 *= Double(2)/(Double(group_permut*(Nc*Nc-1))*lattice_volume);

    //lattice momentum
    Double qq=sin(ipx*PI/latt_size[0])*sin(ipx*PI/latt_size[0]);
    qq += sin(ipy*PI/latt_size[1])*sin(ipy*PI/latt_size[1]);
    qq += sin(ipz*PI/latt_size[2])*sin(ipz*PI/latt_size[2]);
    qq += sin(ipt*PI/latt_size[3])*sin(ipt*PI/latt_size[3]);
    qq = sqrt( Double(4) * qq );

    //normal momentum
    Double pp = Double(4)*(ipx*PI/latt_size[0])*(ipx*PI/latt_size[0]);
    pp += Double(4)*(ipy*PI/latt_size[1])*(ipy*PI/latt_size[1]);
    pp += Double(4)*(ipz*PI/latt_size[2])*(ipz*PI/latt_size[2]);
    pp += Double(4)*(ipt*PI/latt_size[3])*(ipt*PI/latt_size[3]);
    pp = sqrt(pp);    

    //int and double vectors to print
    int Doublesize = 18;
    int intsize = Nd;
    Double arraytoprint_double[Doublesize] = {qq, pp, gd, E, F, G, H, I, J, K, L, A, B, D_orig, D_recons, D_recons2, D_recons_cont, D_recons_cont2};
    int arraytoprint_int[intsize] = {ipx, ipy, ipz, ipt};
    //imp
    Double arraytoprint_double_imp[Doublesize] = {qq, pp, gd, im_E, im_F, im_G, im_H, im_I, im_J, im_K, im_L, im_A, im_B, D_orig, im_D_recons, im_D_recons2, im_D_recons_cont, im_D_recons_cont2};

    QDPIO::cout << "GP - naive:  ";
    for (int column = 0; column < intsize; column++){
      QDPIO::cout << arraytoprint_int[column] << "  ";
    }
    for (int column = 0; column < Doublesize; column++){
      QDPIO::cout << arraytoprint_double[column] << "  ";
    }
    QDPIO::cout << "\n";
    //imp
    QDPIO::cout << "GP - improved:  ";
    for (int column = 0; column < intsize; column++){
      QDPIO::cout << arraytoprint_int[column] << "  ";
    }
    for (int column = 0; column < Doublesize; column++){
      QDPIO::cout << arraytoprint_double_imp[column] << "  ";
    }
    QDPIO::cout << "\n";
  }}}} //ipx for
}

void writeGauge_singleprec(XMLBufferWriter& file_xml,
    XMLBufferWriter& record_xml, 
    const multi1d<LatticeColorMatrix>& u, 
    const string& file, 
    QDP_volfmt_t volfmt, 
    QDP_serialparallel_t serpar)
{
  QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
  if (to.bad())
  {
    QDPIO::cerr << __func__ << ": error writing file " << file << endl;
    QDP_abort(1);
  }


  multi1d<LatticeColorMatrixF> u_f(u.size());
  for(int mu=0; mu < u.size(); ++mu)
    u_f[mu] = u[mu];

  write(to,record_xml,u_f);         // Writing in single precision

  close(to);
}
