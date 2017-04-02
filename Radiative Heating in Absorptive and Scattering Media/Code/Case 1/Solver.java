import java.lang.Math;
import java.io.*;

class Solver {

	final double Thot = 1000;//degrees Kelvin
	final double Tcool = 500;//degrees Kelvin
	
	//I+ Matrices
	public double[][] A1;
	public double[] b1;
	
	//I- Matrices
	public double[][] A2;
	public double[] b2;
	
	//I Matrices
	public double[][] A3;
	public double[] b3;
	
	public double[] solution;
	public int nodecount;
	
	
	private double sigma,margin;
	private double[] temp;

	private double Aconst1;
	private double Aconst2;
	private double Bconst;
	private double lambda;
	private double T;
	private double dels;
	
	private Radiation rad;
	
	//////////////////////////////////////////
	//////////////Constructor/////////////////
	//////////////////////////////////////////
	
	public Solver(int nodes){
		margin = (double)0;
		sigma = (double)0;
		temp = new double[1];
		rad = new Radiation();
		nodecount = nodes;
	}
	
	//////////////////////////////////////////////
	////////////////SOLVER METHODS////////////////
	//////////////////////////////////////////////
	
	//Constants relevant for the Differencing Equations of the I+ & I- Matrices
	public void DiffEqConstants(double beta, double thickness){
		dels = thickness/(nodecount+1);//divide entire thickness over discrete node spacings
		Aconst1 = 1/(4*dels);
		Aconst2 = beta;
	}
	
	//Constants relevant for the BlackBody terms in the B vector
	public void BlackBodyConstants(double wavelength, double kappa){
		Bconst = kappa;
		lambda = wavelength;//In Microns 
	}
	
	
	public void AssembleMatrices(){
		
		///////////////////////////////////////////////////////////////////////////
		////----------Initialize and Assemble the Matrices for I+ -----------//////
		///////////////////////////////////////////////////////////////////////////
		
		//-----------------A-Matrix---------------//
		A1 = new double[nodecount][nodecount];
		
		// Row[0]
		A1[0][0] = Aconst2;//  beta
		A1[0][1] = Aconst1;//  1/4*dels
		
		// Row[1] to Row[nodecount-2]
		for (int i = 1; i<nodecount-1; i++) {
			A1[i][i-1] = -Aconst1;// [-1/(4*dels)]
			A1[i][i] = Aconst2;//    [beta]
			A1[i][i+1] = Aconst1;//  [1/(4*dels)]
		}
		
		// Row[nodecount-1] -- First Order Backward Difference Derivative
		A1[nodecount-1][nodecount-2] = -2*Aconst1;//  [-1/(2*dels)]
		A1[nodecount-1][nodecount-1] = Aconst2 + 2*Aconst1;//   [beta + 1/(2*dels)]
		
		////----------Initialize and Assemble the B-Matrix-----------//////
		
		b1 = new double[nodecount];
		
		T = Tcool + (Thot-Tcool)/(nodecount+2);
		//System.out.print(T);
		
		//First Element
		b1[0] = rad.Planck(T,lambda,1)*Bconst/Math.PI + Aconst1*rad.Planck(Tcool,lambda,1)/Math.PI;
		
		//Second Element to Row[nodecount-2]
		for (int i = 1; i<nodecount-1; i++) {
			T = Tcool + (i+1)*(Thot-Tcool)/(nodecount+2);
															
			b1[i] = rad.Planck(T,lambda,1)*Bconst/Math.PI;
		}
		
		T = Tcool + (nodecount)*(Thot-Tcool)/(nodecount+2);
					
		//Last Element
		b1[nodecount-1]= rad.Planck(T,lambda,1)*Bconst/Math.PI;
		
		///////////////////////////////////////////////////////////////////////////
		////----------Initialize and Assemble the Matrices for I- -----------//////
		///////////////////////////////////////////////////////////////////////////
		
		//-----------------A-Matrix---------------//
		A2 = new double[nodecount][nodecount];
		
		// Row[0]
		A2[0][0] = Aconst2+2*Aconst1;//  beta + 1/(2*dels)
		A2[0][1] = -2*Aconst1;//  [-1/(2*dels)]
		
		// Row[1] to Row[nodecount-2]
		for (int i = 1; i<nodecount-1; i++) {
			A2[i][i-1] = Aconst1;// [1/(4*dels)]
			A2[i][i] = Aconst2;//    [beta]
			A2[i][i+1] = -Aconst1;//  [-1/(4*dels)]
		}
		
		// Row[nodecount-1]
		A2[nodecount-1][nodecount-2] = Aconst1;//  [1/(4*dels)]
		A2[nodecount-1][nodecount-1] = Aconst2;//   [beta]
		
		////----------Initialize and Assemble the B-Matrix-----------//////
		
		b2 = new double[nodecount];
		
		T = Tcool + (Thot-Tcool)/(nodecount+2);
										
		//First Element
		b2[0] = rad.Planck(T,lambda,1)*Bconst/Math.PI;
		
		//First Element to Row[nodecount-2]
		for (int i = 1; i<nodecount-1; i++) {
			T = Tcool + (i+1)*(Thot-Tcool)/(nodecount+2);
			b2[i] = rad.Planck(T,lambda,1)*Bconst/Math.PI;
		}
		
		T = Tcool + (nodecount)*(Thot-Tcool)/(nodecount+2);
		
		//Last Element
		b2[nodecount-1] = rad.Planck(T,lambda,1)*Bconst/Math.PI + Aconst1*rad.Planck(Thot,lambda,1)/Math.PI;
		
	}//End of Assemble Matrices
	
	
	/////SOLVES THE CONSTRUCTED I+ SYSTEM
	public double[] SolveIPlus(){
		double[] out = new double[nodecount];
				
		temp = Thomas(A1,b1);//length = nodecount
		
		for(int i = 0;i<nodecount;i++){
			out[i]=temp[i];
		}
		
		return out;
	}
	
	/////SOLVES THE CONSTRUCTED I- SYSTEM
	public double[] SolveIMinus(){
		double[] out = new double[nodecount];
		
		temp = 	Thomas(A2,b2);				
		
		for(int i = 0;i<nodecount;i++){
			out[i]=temp[i];
		}
				
		return out;
	}
	
	/*
	
	/////SOLVES THE CONSTRUCTED I SYSTEM
	public double[] SolveI(){
		double[] out = new double[nodecount+2];

		temp = 	Thomas(A3,b3);				
		
		out[0]=rad.Planck(Tcool,lambda,1)/Math.PI;
		
		for(int i = 1;i<nodecount+1;i++){
			out[i]=temp[i-1];
		}
		
		out[nodecount+1]=rad.Planck(Thot,lambda,1)/Math.PI;
		
		
		return out;
	}
	
	*/
	 
	/////////////////////////////////////////////////////
	////////////////////Internal Methods/////////////////
	/////////////////////////////////////////////////////
	
	/*
	//////----Successive Over-Relaxation Matrix Solver----/////////
	public double[] SOR(double[][] AA, double[] bb, double tol, double guess, double w){
		//////Uses the Successive Over Relaxation Method to Solve the system Ax=b

		int n = AA.length;
		double[] xk_1 = new double[n];		
		double[] xk = new double[n];		
		
		
		//Initializes output vector with initial guess
		for (int i=0; i<n; i++) {
			xk_1[i]=guess;
			xk[i]=guess*2;
		}
		
		margin = 100;//arbitrarily large initial percentile error
		
		while (margin>tol) {//checks convergence
			
			for (int i =0; i<n; i++) {
				xk_1[i]=xk[i]; //assigns kth-iteration of x to xk_1 before next iteration
			}
			
			for (int i=0; i<n; i++) {
				
				sigma = 0; 
				
				for (int j=0; j<n; j++) {

					if (j!=i){sigma = sigma + AA[i][j]*xk[j];}
				
				}//end of j-loop
			
				xk[i]=(1-w)*xk[i]+(w/AA[i][i])*(bb[i]-sigma);//updates xk with new values
				
			}//end of i-loop: kth iteration has completed
		
			margin = maxerror(xk_1,xk);
			
		}//end of while loop
		
		return xk;	
	}
	*/
	
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
	
	///////----Maximum percentile difference between elements----//////////
	public double maxerror(double[] x1, double[] x2){
		double out=0; 
		double diff = 0;

		for(int i=0; i<x1.length; i++){
			diff=100*Math.abs(((x1[i]-x2[i])/x1[i]));
			
			if(diff>out){out=diff;}//replaces lower error with greater error
			
		}
		return out;//returns greatest error
	}
	
	
}


