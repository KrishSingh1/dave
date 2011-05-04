///////////////////////////////////////////////////////////////////////////////////////
/// \file bvoc.cpp
/// \brief The BVOC module
///
/// Calculation of VOC production and emission by vegetation.
///
/// \author Guy Schurgers (using Almut's previous attempts)
/// $Date: $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework
//   header file should define all classes used as arguments to functions in
//   the present module. It may also include declarations of global functions,
//   constants and types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by
//   the present one;
//   (3) type definitions, constants and file scope global variables for use
//   within the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to
//   other modules or to the calling framework should be declared in the module
//   header file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models
// (frameworks). When porting between frameworks, the only change required
// should normally be in the "#include" directive referring to the framework
// header file.

#include "config.h"
#include "bvoc.h"
#include "canexch.h"
#include "q10.h"


// INITBVOC
// called from FRAMEWORK

void initbvoc(Pftlist& pftlist){

  // initialising the VOC calculations: calculating the fraction of electrones
  // available for isoprene production from the isoprene and monoterpene
  // emission capacities (at T=30oC and Q=1000 umol m-2 s-1) as given by e.g.
  // Guenther et al. (1997) for all PFTs

  double tau_s;
  double gammastar_s;
  double pi_co2_s;
  double apar_s;
  double tscal_s;
  double c1_s;
  double c2_s;
  double phi_pi_s;
  double a_Y_s;
  double je_s;
  double J_s;
  double sla_loc;
  double J_s_leaf;
  double Y_eps_iso;
  double Y_eps_mon;
  double rd_g_s;
  double ko_s;
  double kc_s;
  double sigma_c3_s;
  double sigma_c4_s;

  const double Q10TAU=0.57;
  const double TAU25=2600.0; // value of tau at 25 deg C
  const double THETA=0.7; // colimitation (shape) parameter
  const double Q10KO=1.2;
  const double Q10KC=2.1;
  const double KO25=3.0E4; // value of ko at 25 deg C (Pa)
  const double KC25=30.0; // value of kc at 25 deg C (Pa)

  LookupQ10 lookup_tau(Q10TAU,TAU25);
  LookupQ10 lookup_ko(Q10KO,KO25);
  LookupQ10 lookup_kc(Q10KC,KC25);

  double Tstand=30.; // standard temperature, oC
  double Qstand=1.e-3; // radiation, mol m-2 s-1
  double CO2stand=370.; // standard CO2 concentration, ppm

  double fpar_s=1.;
  double LAMBDA_SC4=0.4;
  double Cfrac=0.5; // mass fraction of C in leaves
  double BC3=.015;
  double BC4=.020;
  double frabs_Q=0.42; // fraction of light absorbed in the first canopy layer
                      // (for standard measurements), 25% to 35% (Almut)

  dprintf("\nInitialising VOC calculations\n");

  // loop through PFTs
  pftlist.firstobj();
  while (pftlist.isobj){
    Pft& pft=pftlist.getobj();

    // conversion factor alpha (converting electron flux into isoprene)

    if(pft.pathway==C3){
      // tau (experimentally determined parameter
      tau_s=lookup_tau[Tstand];
      // CO2 compensation point (Pa)
      // GS21082006 in Almut's version, gammastar_s is calculated for C3
      // photosynthesis only
      gammastar_s=PO2/(2.*tau_s);
      // intercellular CO2 partial pressure (Pa)
      // assuming that stomatal opening is always optimal
      pi_co2_s=pft.lambda_max*CO2stand*1.e-6*1.e5;
      // conversion factor to convert electron flux into isoprene
      // equivalents (-), Niinemets et al., 1999, equation 2
      a_Y_s=(pi_co2_s-gammastar_s)/(6.*(4.67*pi_co2_s+9.33*gammastar_s));
    }
    else{
      // assuming a_Y_s for C4 plants to be calculated as C3, but without
      // oxygenation losses
      a_Y_s=1./(6.*4.67);
    }

    // total electron transport for the standard conditions

    // apar for the standard condition, mol m-2 h-1
    apar_s=frabs_Q*Qstand*3600.*fpar_s;

    // calculation of temperature coefficient
    tscal_s=(1.-0.01*exp(4.6/(pft.pstemp_max-pft.pstemp_high)*
			 (Tstand-pft.pstemp_high)))/
      (1.+exp(4.6*((pft.pstemp_min+pft.pstemp_low)/2.-Tstand)/
	      ((pft.pstemp_min+pft.pstemp_low)/2.-pft.pstemp_min)));
    // pathway-dependent calculation of C1
    if(pft.pathway==C3)
      {
	// mitochondrial respiration, g C m-2 h-1
	ko_s=lookup_ko[Tstand];
	kc_s=lookup_kc[Tstand];
	c1_s=tscal_s*ALPHA_C3*(pi_co2_s-gammastar_s)/
	  (pi_co2_s+2.0*gammastar_s);
	c2_s=(pi_co2_s-gammastar_s)/(pi_co2_s+kc_s*(1.0+PO2/ko_s));
	sigma_c3_s=sqrt(max(0.,1.-(c2_s-2.*BC3)/(c2_s-THETA*2.*BC3)));
	rd_g_s=(c1_s/c2_s*(4.*THETA*BC3-2.*BC3-4.*THETA*BC3*sigma_c3_s+
			   c2_s*sigma_c3_s)*apar_s*12.)*.5;

	// photosynthesis rate, g C m-2 h-1
	je_s=ALPHA_C3*(pi_co2_s-gammastar_s)/(pi_co2_s+2.*gammastar_s)
	  *tscal_s*apar_s*12.;
	J_s=(je_s+rd_g_s)*(4.*pi_co2_s+8.*gammastar_s)/(pi_co2_s-
							gammastar_s);

      }
    else
      {
	// mitochondrial respiration, g C m-2 h-1
	sigma_c4_s=sqrt(max(0.0,1.0-(1.0-2.*BC4)/(1.-THETA*2.*BC4)));
	rd_g_s=(tscal_s*ALPHA_C4*(4.*THETA*BC4-2.*BC4-4.*THETA*BC4*
				  sigma_c4_s+sigma_c4_s)*apar_s*12.)*.5;

	phi_pi_s=min(pft.lambda_max/LAMBDA_SC4,1.);
	// photosynthesis rate, g C m-2 h-1
	je_s=ALPHA_C4*phi_pi_s*tscal_s*apar_s*12.;
	J_s=4.*(je_s+rd_g_s);
      }

    // convert SLA from m2 kg-1 to m2 g-1 leaf
    sla_loc=pft.sla*Cfrac*1.e-3;
    // photosynthesis rate per leaf mass, g C g-1 leaf h-1
    J_s_leaf=J_s*sla_loc;

    // electron fraction assigned to isoprene and monoterpenes for the
    // standard case
    Y_eps_iso=pft.eps_iso*1.e-6/(a_Y_s*J_s_leaf);
    Y_eps_mon=pft.eps_mon*1.e-6/(a_Y_s*J_s_leaf);

    pft.Y_eps_iso=Y_eps_iso;
    pft.Y_eps_mon=Y_eps_mon;

    pftlist.nextobj();
  }
}


// DAYT
// called from BVOC

double dayT(double temp, double daylength, double tempamp){

  // convert daily temperature to daytime temperature using the actual
  // amplitude the conversion assumes a perfect sine function for
  // temperature, with the maximum temperature reached at noon (derivation by
  // Colin Prentice, checked)

  double hdl; // daylength expressed as fraction of PI (rad)
  double dtT; // daytime temperature (oC)

  hdl=daylength*3.14/24.;
  dtT=temp+(tempamp/2)*sin(hdl)/hdl;
  return dtT;
}


// ISOPRMONOT1
// called from VOCCALC

void isoprmonot1(double co2, double agdd5, double temp, double daylength,
                 const PhotosynthesisResult& photosynthesis,
                 const Pft& pft,double fvocseas,double temprel,
                 double& dmonstor, double& iso, double& mon){

  // Calculation of isoprene and monoterpene emissions coupled to
  // photosynthesis as described in Arneth et al. (2007) for isoprene and
  // Schurgers et al. (2009) for monoterpenes.

  double pstemp_max=pft.pstemp_max;
  double pstemp_min=pft.pstemp_min;
  double pstemp_high=pft.pstemp_high;
  double pstemp_low=pft.pstemp_low;

  double rd_g       = photosynthesis.rd_g;
  double pi_co2_opt = photosynthesis.pi_co2_opt;
  double gammastar  = photosynthesis.gammastar;
  double apar       = photosynthesis.apar;
  double phi_pi     = photosynthesis.phi_pi;

  double f_co2;         // CO2 scaling factor for BVOC production (-)
  double f_temp;        // temperature scaling factor for BVOC production (-)
  double f_season;      // seasonality scaling factor for isoprene production
                        // (-)
  double tscal;         // temperature scaling coefficient
  double je_opt;        // PAR-limited photosynthesis rate under optimum
                        // conditions (molC/m2/h)
  double J_opt;         // photosynthetic electron transport (mol/m2/h)
  double a_Y_opt;       // terpenoid electron transport (mol/m2/h)
  double Y_eps_eff_iso; // effective electron fraction for isoprene production
                        // (-)
  double Y_eps_eff_mon; // effective electron fraction for monoterpene
                        // production (-)
  double tcstor;        // time coefficient of monoterpene storage (d)
  double f_temp_mstor;  // temperature scaling factor for monoterpene storage
                        // release (-)

  double co2stan=370.;  // standard atmospheric CO2 concentration (ppm)
  double Tstand=30.;    // standard temperature (oC)

  double tcstor_s=80.;  // time constant for monoterpene storage under standard
                        // temperature (d)
  double tcstor_max=365.;// maximum time constant for monoterpene storage (d)
  double tcstor_min=2.; // minimum time constant for monoterpene storage (d)
  double q10_mstor=1.9; // Q10 value for monoterpene storage, -
  double f_tempmax=2.3; // maximum temperature scaling factor for BVOC
                        // production (-)
  double epsT=0.1;      // temperature sensitivity of BVOC production (-)

  // CO2 scaling factor, -
  f_co2=co2stan/co2;

  // growing season scaling factor (isoprene only), -
  f_season=fvocseas;

  // temperature scaling factor, -
  f_temp=min(f_tempmax,exp(epsT*(temp-Tstand)));

  // effective electron fraction for isoprene and monoterpene production
  // (corrected for CO2 concentration and temperature, and in case of isoprene
  // including a seasonality correction)
  Y_eps_eff_iso=pft.Y_eps_iso*f_co2*f_temp*f_season;
  Y_eps_eff_mon=pft.Y_eps_mon*f_co2*f_temp;

  // optimum photosynthetic electron transport rate (Niinemets et al., 1999,
  // eq. 1)
  if(temp<pstemp_max)
    tscal=(1.-0.01*exp(4.6/(pstemp_max-pstemp_high)*(temp-pstemp_high)))/
      (1.+exp(4.6*((pstemp_min+pstemp_low)/2.-temp)/
	      ((pstemp_min+pstemp_low)/2.-pstemp_min)));
  else
    tscal=0.;

  if(pft.pathway==C3){
    je_opt=ALPHA_C3*(pi_co2_opt-gammastar)/(pi_co2_opt+2.*gammastar)*
      tscal*apar*CQ*12.;
    J_opt=(je_opt+rd_g)*(4.*pi_co2_opt+8*gammastar)/(pi_co2_opt-gammastar);
    a_Y_opt=(pi_co2_opt-gammastar)/(6.*(4.67*pi_co2_opt+9.33*gammastar));
  }
  else{
    je_opt=ALPHA_C4*phi_pi*tscal*apar*CQ*12.;
    J_opt=4.*(je_opt+rd_g);
    a_Y_opt=1./(6.*4.67);
  }

  // isoprene production, g C m-2 d-1
  iso=Y_eps_eff_iso*J_opt*a_Y_opt;
  // monoterpene production, g C m-2 d-1
  // (only the production part is given here)
  mon=Y_eps_eff_mon*J_opt*a_Y_opt;

  // release from monoterpene storage, g C m-2 d-1
  f_temp_mstor=1./pow(q10_mstor,((temprel-Tstand)/10.));
  tcstor=tcstor_s*f_temp_mstor;
  tcstor=min(tcstor,tcstor_max);
  tcstor=max(tcstor,tcstor_min);
  dmonstor=1./tcstor;
}


// LEAFT
// called from BVOC

double leafT(double temp, double daylength, double adtmm, double co2,
	     double lambda, double eet, double ga, double rs_day, double gmin,
	     double lai){

  // Canopy temperature is calculated from the air temperature and leaf
  // latent heat loss, using a weighted average temperature within the canopy.
  // Revised version compared to Arneth et al. (2007).

  double gc;         // canopy conductance for water vapour, mm s-1
  double trans;      // transpiration, mm s-1
  double delT;       // temperature increase (air --> canopy), K
  double rs_net;     // net shortwave radiation, W m-2
  double canopytemp; // canopy temperature, oC
  double lhloss;     // latent heat loss, W m-2

  double sigma=5.67e-8;  // Stefan-Boltzmann constant, W m-2 K-4
  double emiss_leaf=.97; // average emissivity for leaves, Campbell and Norman,
                         // 1998
  double lam=44100.;     // latent heat loss of vapourization, J mol-1 (for
                         // 20oC)
  double ALPHAM=1.4;
  double GM=5.0;
  double rhoair=1.204;   // air density, kg m-3
  double cp=1010.;       // specific heat capacity of air, J kg-1 K-1

  // canopy conductance for water vapour (mm s-1)
  gc=gmin+adtmm/(daylength*3600.)*1.6/(co2*1.e-6*(1.-lambda));

  // transpiration, corrected for the fraction of the ground covered by
  // vegetation (mm s-1)
  trans=eet/(daylength*3600.)*ALPHAM*(1.-exp(-gc/GM))*
    (1.-exp(-LAMBERTBEER_K*lai));

  // net shortwave radiation (W m-2)
  rs_net=rs_day/(3600.*daylength);

  // latent heat loss of the whole canopy (W m-2)
  lhloss=trans*lam*1.e3/18.;

  // temperature increase (oC)
  if(lai>1.e-2){
    delT=(rs_net-(lhloss/(1.-exp(-LAMBERTBEER_K*lai))))/
      (4*emiss_leaf*sigma*pow((temp+273.15),3.)+rhoair*cp*ga)*
      (1.+exp(-LAMBERTBEER_K*lai))/2.;
  }
  else{
    delT=0.;
  }

  // canopy temperature, (oC)
  canopytemp=temp+delT;

  return canopytemp;
}



// VOCSEAS
// called from VOCCALC

double vocseas(double& f_season, double temp, double daylength, int nday,
	       double agdd5, const Pft& pft){

  // calculating the seasonality for VOCs (isoprene and monoterpene) for PFTs
  // Revised version compared to Arneth et al. (2007).
  // Seasonality switch pft.seas_iso read in from .ins file

  double vocgdd5ramp; // GDD sum required for full expression of isoprene
                      // synthase/isoprene production
  double rdr=0.05;    // relative decay rate (d-1)
  double tmin=5.;     // minimum temperature for end of growing season (oC)
  double dmin=11.;    // minimum daylength for end of growing season (h)

  // required GDD sum for VOCs is assumed to be twice as large as for
  // phenology
  double mulgdd=2.;
  vocgdd5ramp=mulgdd*pft.phengdd5ramp;

  if(pft.seas_iso){
    f_season=1.;
  }
  else{
    if(agdd5<=vocgdd5ramp){
      if(vocgdd5ramp!=0.){
	f_season=agdd5/vocgdd5ramp;
      }
      else{
	f_season=1.;
      }
    }
    else if(temp<tmin || daylength<dmin){
      f_season=f_season*pow((1.-rdr),((double)nday));
    }
    else{
      f_season=1.;
    }
  }
  return f_season;
}


// BVOC
// called from NPP

void bvoc(double daylength, double temp, double tempamp,
     const PhotosynthesisResult& photosynthesis,
     double co2, double lambda, double eet, double agdd5, int nday,
     double rs_day, double lai, const Pft& pft,
     double& iso, double& mon, double& dmonstor, double& dleaftemp,
     double& fvocseas){

  // Calculation of isoprene and monoterpene production in leaves as a function
  // of photosynthetis. Isoprene and monoterpenes are calculated from a
  // standardized fraction of the total photosynthesis, which is adjusted as
  // a function of temperature, CO2 concentration and (for isoprene)
  // seasonality.
  //
  // Isoprene calculations following Arneth et al. (2007), monoterpene
  // calculations following Schurgers et al. (2009). Compared to the original
  // publications, adjustments were made in the calculation of leaf
  // temperatures (which account now for longwave radiation as well, and
  // perform a weighted averaging over the whole canopy), and in the calculation
  // of isoprene seasonality (which requires a GDD sum twice as large as
  // required for phenology, and decreases with a relative rate at the end of
  // the growing season).

  double temp_leaf_daytime; // leaf temperature during daytime (oC)
  double temp_leaf;         // leaf temperature (diurnal average) (oC)

  // perform daily to daytime correction
  temp_leaf_daytime=dayT(temp,daylength,tempamp);

    // perform air temperature to leaf temperature correction
  temp_leaf=leafT(temp,daylength,photosynthesis.adtmm,co2,lambda,eet,pft.ga,
		  rs_day,pft.gmin,lai);
  temp_leaf_daytime=leafT(temp_leaf_daytime,daylength,photosynthesis.adtmm,co2,
			  lambda,eet,pft.ga,rs_day,pft.gmin,lai);

  // calculate seasonality for VOC emissions, -
  fvocseas=vocseas(fvocseas,temp,daylength,nday,agdd5,pft);

  // determine temperature difference between leaf temperature and daytime
  // temperature for output
  dleaftemp=temp_leaf-temp;

  // calculate isoprene and monoterpene emissions, g C m-2 d-1
  isoprmonot1(co2,agdd5,temp_leaf_daytime,daylength,photosynthesis,
              pft,fvocseas,temp_leaf,
              dmonstor,iso,mon);

  // convert from g C m-2 d-1 to mg C m-2 d-1
  iso=iso*1.e3;
  mon=mon*1.e3;


}





//
// REFERENCES
//
// Arneth, A., Niinemets, Ü, Pressley, S., Bäck, J., Hari, P., Karl, T., Noe,
//   S., Prentice, I.C., Serca, D., Hickler, T., Wolf, A., Smith, B., 2007.
//   Process-based estimates of terrestrial ecosystem isoprene emissions:
//   incorporating the effects of a direct CO2-isoprene interaction.
//   Atmospheric Chemistry and Physics, 7, 31-53.
// Campbell, G.S., Norman, J.M., 1998. An introduction to environmental
//   biophysics, 2nd edition. Springer, New York.
// Guenther, A., 1997. Seasonal and spatial variations in natural volatile
//   organic compound emissions. Ecological Applications, 7, 34-45.
// Niinemets, Ü, Tenhunen, J.D., Harley, P.C., Steinbrecher, R., 1999. A model
//   of isoprene emission based on energetic requirements for isoprene
//   synthesisand leaf photosynthetic properties for Liquidambar and Quercus.
//   Plant, Cell and Environment, 22, 1319-1335.
// Schurgers, G., Arneth, A., Holzinger, R., Goldstein, A., 2009. Process-
//   based modelling of biogenic monoterpene emissions combining production
//   and release from storage. Atmospheric Chemistry and Physics, 9, 3409-3423.

