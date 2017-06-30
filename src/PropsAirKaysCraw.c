
#include "PropsAirKaysCraw.h"
#include "Misc.h"

static const double Rgas = 286.9; //air gas constant J/kg-K

const static double T1 = 4400 / 1.8;  //temp used for extrapolation only
const static double T2 = 4500 / 1.8;  //maximum temp (kelvin) used in curve fit

double CpAir(const double T) {
        
        //   Correlation of Crawford Data for Specific Heat of Air
        //   Data Source: "Convective Heat and Mass Transfer", Text by Kays and Crawford
        //   Note: Valid in range 180 to 4500 degrees R

		double cp, cp_out;
        double y1, y2;
        
        //convert temp from Kelvin to Rankine for now since curve fit was generated using english units
        double temp = T * 1.8;
        
        //Polynomial curve-fit constants
        double A[6];
        A[0] = 0.2567471;
        A[1] = -0.00008577731;
        A[2] = 0.0000001373991;
        A[3] = -0.00000000008044221;
        A[4] = 0.00000000000002469929;
        A[5] = -3.925129E-18;
        A[6] = 2.605416E-22;

        if ( temp <= T2 ) {
            cp = PolyFit(temp, A);
         } else {
         	//temp is greater than T2
            y1 = PolyFit(T1, A);
            y2 = PolyFit(T2, A);
            cp = (y2 - y1) * (temp - T1) / (T2 - T1) + y1;
        }       
        
        //cp has been calculated in btu/lbm-R, convert to J/kg-K
        cp_out = cp * 4184; //to J/kg-K
        
        return cp_out;
}

double EnthalpyAir(const double T) {
       
        //   Correlation of Crawford Data for Enthalpy of Air
        //   used to get Enthalpy of air (by integrating polynomial from Cp [Cp= dh/dt])

        //   Note: Valid in range 180 to 4500 degrees R

		double h;
        double y1, y2;
        
        //convert temp from Kelvin to Rankine for now since curve fit was generated using english units
        double temp = T * 1.8;
        
        //Polynomial curve-fit constants
        double Aprime[7];
        Aprime[0] = 0.0f;
        Aprime[1] = 0.2567471 / 1;
        Aprime[2] = -0.00008577731 / 2;
        Aprime[3] = 0.0000001373991 / 3;
        Aprime[4] = -0.00000000008044221 / 4;
        Aprime[5] = 0.00000000000002469929 / 5;
        Aprime[6] = -3.925129E-18 / 6;
        Aprime[7] = 2.605416E-22 / 7;

        if (temp <= T2) {
            h = PolyFit(temp, Aprime);
        } else {
            y1 = PolyFit(T1, Aprime);
            y2 = PolyFit(T2, Aprime);
            h = (y2 - y1) * (temp - T1) / (T2 - T1) + y1;
        }
        
        //h calculated in btu/lbm, convert to J/kg
        double h_out = h * 2326;
        
		return h_out;
 }

double GammaAir(const double T) {
        
        //   Correlation of Crawford Data for Gamma of Air
        //   Data Source: "Convective Heat and Mass Transfer", Text by Kays and Crawford
        //   Note: Valid in range 180 to 4500 degrees R

		double gam;
        double y1, y2;
        
        //convert temp from Kelvin to Rankine for now since curve fit was generated using english units
        double temp = T * 1.8;
        
        //Polynomial curve-fit constants
        double A[6];
        A[0] = 1.362246;
        A[1] = 0.0001949838;
        A[2] = -0.0000003178944;
        A[3] = 0.0000000001905998;
        A[4] = -0.00000000000005777532;
        A[5] = 8.785717E-18;
        A[6] = -5.377543E-22;

        if (temp <= T2 ) {
           gam = PolyFit(temp, A);
        } else {
            y1 = PolyFit(T1, A);
            y2 = PolyFit(T2, A);
            gam = (y2 - y1) * (temp - T1) / (T2 - T1) + y1;
        }
        
        //gamma is unitless so no conversion necessary, just check that temp is in Rankine
		return gam;
}

double kAir(const double T) {
        
        //   Correlation of Crawford Data for Thermal Conductivity of Air
        //   Data Source: "Convective Heat and Mass Transfer", Text by Kays and Crawford
        //   Note: Valid in range 180 to 4500 degrees R
		// english units BTU/(s*ft*deg. R)

        double k;
        double y1, y2;

        //convert temp from Kelvin to Rankine for now since curve fit was generated using english units
        double temp = T * 1.8; 
        
        //Polynomial curve-fit constants
        double A[6];
        A[0] = -0.0000002575394;
        A[1] = 0.00000001038856;
        A[2] = -0.000000000005120556;
        A[3] = 0.000000000000002387962;
        A[4] = -5.654767E-19;
        A[5] = 5.395641E-23;
        A[6] = 0;

        if (temp <= T2) {
            k = PolyFit(temp, A);
        } else {
            y1 = PolyFit(T1, A);
            y2 = PolyFit(T2, A);
            k = (y2 - y1) * (temp - T1) / (T2 - T1) + y1;
        }
        
        //thermal conductivity is calculated in btu/s-ft-R, convert to J/s-m-K
        double k_out = k * 44.371;
        
       return k_out;
}

double PrAir(const double T) {
        
        //   Correlation of Crawford Data for Prandtl Number for Air
        //   Data Source: "Convective Heat and Mass Transfer", Text by Kays and Crawford
        //   Note: Valid in range 180 to 4500 degrees R

        //convert temp from Kelvin to Rankine for now since curve fit was generated using english units
        double temp = T * 1.8;

		double pr;
        double y1, y2;
        
        //Polynomial curve-fit constants
        double A[6];
        A[0] = 0.869641;
        A[1] = -0.0005661112;
        A[2] = 0.0000007014822;
        A[3] = -0.0000000004185312;
        A[4] = 0.0000000000001298091;
        A[5] = -2.014177E-17;
        A[6] = 1.231103E-21;

        if (temp <= T2) {
            pr = PolyFit(temp, A);
        } else {
            y1 = PolyFit(T1, A);
            y2 = PolyFit(T2, A);
            pr = (y2 - y1) * (temp - T1) / (T2 - T1) + y1;
        }
        
        //pradtl number is non-dimensional, just check that temp has been converted
        return pr;
}

double ViscAir(const double T) {
        
        //   Correlation of Crawford Data for Dynamic Viscosity of Air
        //   Data Source: "Convective Heat and Mass Transfer", Text by Kays and Crawford
        //   Note: Valid in range 180 to 4500 degrees R

        double mu;
        double y1, y2;

        //convert temp from Kelvin to Rankine for now since curve fit was generated using english units
        double temp = T * 1.8;

        //Polynomial curve-fit constants
        double A[6];
        A[0] = 0.00000003028189;
        A[1] = 0.00000002911495;
        A[2] = -0.00000000001383959;
        A[3] = 0.000000000000004897591;
        A[4] = -8.809283E-19;
        A[5] = 6.235209E-23;
        A[6] = 0;

        if (temp <= T2) {
            mu = PolyFit(T, A);
        } else {
            y1 = PolyFit(T1, A);
            y2 = PolyFit(T2, A);
            mu = (y2 - y1) * (temp - T1) / (T2 - T1) + y1;
        }
        
        //mu calculated in terms of lbm/ft-s, convert to kg/m-s or Pa-s or N-s/m^2
        double mu_out = mu * 1.488163944;
        
        return mu_out;
}

double DensityAir(const double T, const double P) {
		double rho = P / (T * Rgas);
        return rho;
}
