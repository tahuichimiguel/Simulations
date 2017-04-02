import java.lang.Math;

class Radiation{

	
	public Radiation(){}
	
	////Returns Black Body Emissive Power at a Given Wavelength and Temperature
	////UNIT REQUIREMENTS: T is in Kelvin
	////UNIT REQUIREMENTS: lambda is in microns
	public double Planck(double T, double lambda, double RI){
		
		double lambda_meters = lambda*Math.pow(10,-6);
		
		final double h = 6.626*Math.pow(10,-34); //[units] = J*s
		final double c0 = 2.998*Math.pow(10,8); //[units] = m/s
		final double k  = 1.3806*Math.pow(10,-23); //[units] = J/K
		
		double C1 = 2*Math.PI*h*c0*c0;
		double C2 = h*c0/(k);
		
		double out = (C1/(RI*RI*Math.pow(lambda_meters,5)))*(1/(Math.exp(C2/(RI*lambda_meters*T))-1));

		//double out = (C1/(RI*RI*Math.pow(lambda,5)))*(1/(Math.exp(C2/(RI*lambda*T))-1));

		
		return out;
		
	}


}