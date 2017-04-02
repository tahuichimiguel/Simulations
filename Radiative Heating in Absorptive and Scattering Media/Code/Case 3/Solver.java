import java.lang.Math;
import java.io.*;

class Solver {

	final double Thot = 1000;//degrees Kelvin
	final double Tcool = 500;//degrees Kelvin
	final double sigma = 5.67*Math.pow(10,-8);//Stefan-Boltzman Constant Watts/(m^2*Kelvin^4)
	
	//I Matrices
	private double[][] A;
	private double[] b ;
	
	private double[] IPlus;
	private double[] IMinus;
	
	private double[] solution;
	
	private int nodecount;
	private double convergence = 0.001; // 0.1%
	
	private double kappa, beta, scattC;
	private double dels;

	
	private double Aconst1;
	private double Aconst2;
	private double Bconst1;
	private double Bconst2;
	
	private double lambda;
	private double T;
		
	//////////////////////////////////////////
	//////////////Constructor/////////////////
	//////////////////////////////////////////
	
	public Solver(int nodes){
		nodecount = nodes;
		IPlus = new double[nodecount];
		IMinus = new double[nodecount];
	}
	
	//////////////////////////////////////////////
	////////////////SOLVER METHODS////////////////
	//////////////////////////////////////////////
	
	//Constants relevant for the Differencing Equations of the I+ & I- Matrices
	public void DiffEqConstants(double SIGMA, double KAPPA, double thickness){
		dels = thickness/(nodecount+1);//divide entire thickness over discrete node spacings
		beta = SIGMA+KAPPA;
		scattC = SIGMA;
		kappa = KAPPA;
		
		Aconst1 = 1/(2*dels);
		Aconst2 = beta;
		
		Bconst1 = 5/8*scattC;
		Bconst2 = 3/8*scattC;
	}
	
	//Iteratively Solves the I+ and I- Equations using an Explicit Euler Scheme
	public void InitializeIPlusIMinus(int maxIter){
		
		double error = 1; //Initializes Error Variable
		double temp=0;
		
		//Zeroth Indexed I+ Element (To right of the cool wall)
		IPlus[0]=(-2*beta*sigma*Math.pow(Tcool,4)/Math.PI + 2*Bconst1*sigma*Math.pow(Tcool,4)/Math.PI)*dels+sigma*Math.pow(Tcool,4)/Math.PI;
		
		//Last Indexed I- Element (To left of the hot wall)
		IMinus[nodecount-1]=(2*beta*sigma*Math.pow(Thot,4)/Math.PI - 2*Bconst2*sigma*Math.pow(Thot,4)/Math.PI)*(-dels)+(sigma*Math.pow(Thot,4)/Math.PI);
		
		//Iterations continue until the maximum error is "converged" OR a maxIter # of iterations has passed
		int count=0;
		//System.out.println(count);
		while(count<maxIter){
			//System.out.println(count);
			//System.out.println("Error = " + error);
			
			if(error<=convergence){
				System.out.println(error);
				break;
			}
			 
			error=0;//reset error to evaluate next iteration
			
			//System.out.println("\n IPLUS");

			//I+ Euler Scheme
			for (int i = 1;i<nodecount; i++){
				//Temperature at Node to the LEFT!!!
				//SINCE	the derivative is dependent on the previous node...the temperature is evaluated at the node to the left!!
				T = Tcool + (i)*(Thot-Tcool)/(nodecount+2);

				temp = (-2*beta*IPlus[i-1] + 2*Bconst1*IPlus[i-1] + 2*Bconst2*IMinus[i-1])*dels + IPlus[i-1];
				
				if(Math.abs((temp-IPlus[i])/temp) > Math.abs(error)){
					error = Math.abs((temp-IPlus[i])/temp);
				}
			
				//System.out.println(temp);
				//System.out.println(IPlus[i] + "\n");
				
				IPlus[i] = temp;
				//System.out.println(IPlus[i]);
			}
		
			//System.out.println("\n IMINUS");
			
			//I-Euler Scheme
			for (int i = nodecount-2;i>-1; i--){
				T = Tcool + (i+2)*(Thot-Tcool)/(nodecount+2);
		
				// "-dels" indicates integration in the negative direction)
				temp = (2*beta*IMinus[i+1] - 2*Bconst1*IPlus[i+1] - 2*Bconst2*IMinus[i+1])*(-dels)+IMinus[i+1];
			
				if(Math.abs((temp-IMinus[i])/temp) > Math.abs(error)){
					error = Math.abs((temp-IMinus[i])/temp);
				}
			
				//System.out.println(temp);
				//System.out.println(IMinus[i]);
				
				
				IMinus[i]=temp;
				//System.out.println(IMinus[i]);

			}
		
			count++;
		}
		
	}
	
	public double[] getIPlus(){
		double[] out = new double[nodecount];
				
		for(int i = 0;i<nodecount;i++){
			out[i]=IPlus[i];
		}
		
		
		return out;
	}
	
	/////SOLVES THE CONSTRUCTED I- SYSTEM
	public double[] getIMinus(){
		double[] out = new double[nodecount];
				
		for(int i = 0;i<nodecount;i++){
			out[i]=IMinus[i];
		}
				
		return out;
	}
	
	/*
	/////SOLVES THE CONSTRUCTED I SYSTEM
	public double[] SolveI(){
		double[] temp = new double[nodecount+2];
		double[] out = new double[nodecount+2];
		temp = 	Thomas(A,b);
						
		out[0]=sigma*Math.pow(Tcool,4)/Math.PI;
		
		for(int i = 1;i<nodecount+1;i++){
			out[i]=temp[i-1];
		}
		
		out[nodecount+1]=sigma*Math.pow(Thot,4)/Math.PI;
		
		
		return out;
		
	}
	*/
	
	/////////////////////////////////////////////////////
	////////////////////Internal Methods/////////////////
	/////////////////////////////////////////////////////
	
	///////////////////---Thomas Algorithm---////////////////////
	public double[] Thomas(double[][] AA, double[] bb){
		int n = AA.length;
		
		double gamma, rho;
		double[] xx = new double[n];		

		
		//////////////////Forward Elimination Phase/////////////////
		
		for (int i = 1; i < n; i++)
        {
			gamma = AA[i][i-1]/(AA[i-1][i-1]);
			
			AA[i][i] = AA[i][i]-gamma*AA[i-1][i];   //Diagonal Term
			bb[i] = bb[i]-gamma*bb[i-1];
			
		
			
		}
		
		//////////////////Back Substitution Phase/////////////////

		xx[n-1] = bb[n-1]/AA[n-1][n-1];
		
		for (int i = n-2; i>-1; i--)
		{
			xx[i]=(bb[i]-AA[i][i+1]*xx[i+1])/AA[i][i];		
			
		}
			
		return xx;
	}
	
	
	
}


