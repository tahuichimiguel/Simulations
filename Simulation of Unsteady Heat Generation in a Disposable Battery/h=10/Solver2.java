import java.lang.Math;
import java.io.*;

class Solver2 {

	//Domain Parameters
	private double r_max = 0.02;						// cylinder radius - in meters - [fixed constant]
	private double theta_max = Math.PI;				    // angular span of domain {equal to pi for a half cylinder} - [fixed constant]
	private double theta_secondary_span = Math.PI/8;	// angular span of smaller section {angular space of [0 , theta_secondary_span]} - [fixed constant]
		
	private double kmajor = 3;							// conductivity of larger section - in W/(m*K) - [fixed constant]
	private double kminor = 1;							// conductivity of smaller section - in W/(m*K) - [fixed constant]
	private double To = 298.15;							// initial temperature of domain - in Kelvin - [fixed constant]
	private double Cp = 950;							// specific heat - in J/(kg*K) - [fixed constant]
	private double rho = 2618;							// density - in kg/m^3 - [fixed constant]
	private double Tinf = 298.15;						// far-field temperature - in Kelvin - specified [fixed constant]	
	
	private double h;									// heat transfer coefficient - in W/(m^2*K) - [specified constant]

	//Mesh Parameters
	private int r_count;								// # of cells in the radial direction - [specified constant]
	private int theta_count;							// # of cells in angular direction - [specified constant]
	private double del_r;								// radial spacing between cell centers - calculated @ construction [calculated constant]
	private double del_theta;							// angular spacing between cell centers - calculated @ construction [calculated constant]
	private int cellcount;								// total number of control volumes - [calculated constant]

	
	//Simulation Parameters
	private double current_timestep;					//  current timestep
	private double duration = 630;						//	simulation duration - in seconds - [fixed constant]
	private int timesteps;								//	# of timesteps - [specified constant]
	private double del_t;								//	timestep increment - in seconds - [calculated constant]
	
	//Calculation Variables
	private double LocalCoeff;							// Matrix Coefficient for center of cell of interest - calculated @ runtime [calculated constant]
	private double EastCoeff;							// Matrix Coefficient for cell right of cell of interest - calculated @ construction [calculated constant]
	private double WestCoeff;							// Matrix Coefficient for cell left of cell of interest  - calculated @ construction [calculated constant]
	private double NorthCoeff;						// Matrix Coefficient for cell above of cell of interest - calculated @ construction [calculated constant]
	private double SouthCoeff;							// Matrix Coefficient for cell below of cell of interest - calculated @ construction [calculated constant]
	
	private double K[][];								// Conductivity Matrix - [calculated constant matrix]
	private double A[][];								// Stiffness Matrix - [calculated constant matrix]
	private double b[];									//  Variable used repeatedly to construct the b-matrix - calculated @ runtime for each point of interest	
	
	//Output Variables		
	private double[][] T;								//Temperature field - in Kelvin - Dimensions specified in construction
														//T[timestep][mesh_index]
		
	//////////////////////////////////////////
	//////////////Constructor/////////////////
	//////////////////////////////////////////
	
	public Solver2(double h_convect, int radial_bands, int angular_bands, int time_steps){
		//variable assignments
		h = h_convect;
		
		r_count = radial_bands;
		theta_count = angular_bands;
		del_r = r_max/(r_count+1);
		del_theta = theta_max/theta_count;
		cellcount = r_count*theta_count+1;
		
		timesteps = time_steps;
		del_t = (double)(duration/timesteps);
	
		
		T = new double[cellcount][timesteps];
		A = new double[cellcount][cellcount];
		K = new double[r_count][theta_count];
		b = new double[cellcount]; 
	}
	
	//////////////////////////////////////////////
	////////////////SOLVER METHODS////////////////
	//////////////////////////////////////////////
	
	public void ConstructConductivity(){
		double angle;
		
		for (int theta =0; theta<theta_count; theta++) {
			angle = theta*del_theta;
			
			for (int r = 0; r<r_count; r++) {
			
				if(angle<=22.5*Math.PI/180){//minor region
					K[r][theta]=kminor;
					
				}else{
					K[r][theta]=kmajor;//major region
				}
			}
		}
	}
	
	
	public void AssembleStiffnessMatrix(){
		//central cell
		int i = 0;
		double r,rl,t,rc;
		
		rc = del_r;
		r = del_r;
		
		LocalCoeff = Math.PI*rho*Cp/del_t;
		
		for (int j=1; j<cellcount; j+=r_count) {
			//rl = r+del_r/2;
			//t = ((j-r)%theta_count+1)*del_theta;
			
			LocalCoeff+=2*(del_theta/rc)*ksouth(j)/(1.5*del_r);
			
			
			A[0][j] = -2*(del_theta/rc)*ksouth(j)/(1.5*del_r);
		}
		
		A[0][0] = LocalCoeff;
		
		
		////internal cells and boundary cells
		for (i=1; i<cellcount; i++) {
			r = (i%r_count)*del_r; ///radial position of inner boundary of cell
			rl = r+del_r/2;//radial position of central point
			t = ((i-(i%r_count))/r_count+1)*del_theta; //angular position of west cell wall
			
			///boundary cells	
			if(i%r_count==0){
				LocalCoeff = (rho*Cp/del_t) + keast(i)/((r_max-del_r)*rl*del_theta*del_theta) + kwest(i)/((r_max-del_r)*rl*del_theta*del_theta)+
																										ksouth(i)/(del_r*del_r)+ h*r_max/((r_max-del_r)*del_r);
				SouthCoeff = -ksouth(i)/(del_r*del_r);
				EastCoeff = -keast(i)/((r_max-del_r)*rl*del_theta*del_theta);
				WestCoeff = -kwest(i)/((r_max-del_r)*rl*del_theta*del_theta);
				
				A[i][i]= LocalCoeff; 			
				A[i][i-1]= SouthCoeff;			
				//A[i][i+1]= NorthCoeff(i);		
				
				if(i>=r_count+1){	
					A[i][i-r_count]= EastCoeff;
				}
				
				if(i<=theta_count*r_count-r_count){
					A[i][i+r_count]= WestCoeff;
				}
									
				continue;
			}
			
			//interior cells
			NorthCoeff = -knorth(i)*(r+del_r)/(r*del_r*del_r);
			A[i][i+1]= NorthCoeff;		
			
			if (r>del_r) {//not in contact with center
				LocalCoeff = (rho*Cp/del_t) + keast(i)/(r*rl*del_theta*del_theta)+kwest(i)/(r*rl*del_theta*del_theta)+ksouth(i)/(del_r*del_r)+ knorth(i)*(r+del_r)/(r*del_r*del_r);
				SouthCoeff = -ksouth(i)/(del_r*del_r);
				A[i][i]= LocalCoeff;
				A[i][i-1]= SouthCoeff;
			}else{//in contact with center
				LocalCoeff = (rho*Cp/del_t) + keast(i)/(r*rl*del_theta*del_theta)+kwest(i)/(r*rl*del_theta*del_theta)+ksouth(i)/(del_r*1.5*del_r)+knorth(i)*(r+del_r)/(r*del_r*del_r);
				SouthCoeff = -ksouth(i)/(del_r*1.5*del_r);
				A[i][i]= LocalCoeff;
				A[i][0]= SouthCoeff;
			}
			
				
			
			
			if(i>=r_count+1){
				EastCoeff = -keast(i)/(r*rl*del_theta*del_theta);
				A[i][i-r_count] = EastCoeff;
			}
			
			if(i<=theta_count*r_count-r_count){
				WestCoeff = -kwest(i)/(r*rl*del_theta*del_theta);
				A[i][i+r_count] = WestCoeff;
			}
		}
	}
	
	public void Initialize(){
		
		for(int i =0;i<cellcount;i++){
			T[i][0]=To;	
		}
		
	}
	
	public void updateB(double delt,int next_timestep){
		b[0]= T[0][next_timestep-1]*Math.PI*rho*Cp/del_t + qgen(next_timestep*delt)*Math.PI;
		
		for (int i = 1; i<cellcount; i++) {
			if(i%r_count==0)
			{//boundary cell
				b[i]=T[i][next_timestep-1]*rho*Cp/delt + qgen(next_timestep*delt)+h*Tinf*r_max/((r_max-del_r)*del_r);
				continue;
			}
			
			//interior cells
			b[i]=T[i][next_timestep-1]*rho*Cp/del_t+qgen(next_timestep*delt);
		
		}
	}
	
	
	public void Solve(){
		double[] temp;
		
		for (int i=1; i<timesteps; i++) {
			System.out.println( i);
			updateB(del_t, i);
			temp=SOR(A,b,0.05,To,1);//forces convergence to within 0.1%
			
			for (int j =0; j<cellcount; j++){//copy solution to solution matrix
				T[j][i]=temp[j];
			}

		}
	
	}
	
	public double[] getT(int time_step){
		double[] out = 	new double[cellcount];
		
		for (int j =0; j<cellcount; j++){//copy solution to solution matrix
			out[j]=T[j][time_step];
		}
		
		return out;
		
	}
	
	public double keast(int i){
		//i is greater or equal to 1
		int r1,t1,r2,t2;
		r1 = (i-1)%r_count;
		t1 = (i-r1)/r_count;
		r2 = r1;
		t2 = t1-1;		
		
		if (t1>0) {
			return 2*K[r1][t1]*K[r2][t2]/(K[r1][t1]+K[r2][t2]);	
		}else{
			return 0.0;
		}
				
	}
	
	public double kwest(int i){
		//i is greater or equal to 1
		int r1,t1,r2,t2;
		r1 = (i-1)%r_count;
		t1 = (i-r1)/r_count;
		r2 = r1;
		t2 = t1+1;		
		
		if (t1<theta_count-1) {
			return 2*K[r1][t1]*K[r2][t2]/(K[r1][t1]+K[r2][t2]);	
			
		}else{
			return 0.0;
		}
	}
	
	public double ksouth(int i){
		//i is greater or equal to 1
		int r1,t1;
		r1 = (i-1)%r_count;
		t1 = (i-r1)/r_count;
		return K[r1][t1];
	}
	
	public double knorth(int i){//for internal elements
		return ksouth(i);
	}
	
	public void writetoFile(int t_step, double h, int radd, int thetadd){
		double[][] out = new double[cellcount][3];
		
		for (int i = 0; i<cellcount; i++) {
			int r,t;
			double rr, tt;
			double x,y;
			
			r = i%r_count;
			t = (i-r)/r_count;
			
			rr=del_r/2+r*del_r;
			tt=(t+0.5)*del_theta;
			
			x=rr*Math.cos((t+0.5)*del_theta);
			y=rr*Math.sin((t+0.5)*del_theta);
			
			if(i==0){
				x=0.0;
				y=0.0;
			}
			
			
			out[i][0]=x;//x-coordinate
			out[i][1]=y;//y-coordinate
			out[i][2]=T[i][t_step];//temperature
		}
		
		
		Writer pen = new Writer();
		String filename = new String("T-Cartesian, Time_Step = " + t_step + ", h = " + h + ", radd = " +radd+ ", thetadd = " + thetadd);
		pen.writeVector(out,filename);
		
	}
	
	
	
	
	
	////-----Internal Heat Generation Function------////////
	public double qgen(double time){
		double c1,c2,c3,c4,c5,c6,c7;
		c1 = -1.7031*Math.pow(10,-11);
		c2 =  1.3146*Math.pow(10,-7);
		c3 = -1.8011*Math.pow(10,-4);
		c4 = 1.024*Math.pow(10,-1);
		c5 = -2.8133*Math.pow(10,1);
		c6 = 3.6444*Math.pow(10,3);
		c7 = -2.7545*Math.pow(10,3);
		
		double out = c1*Math.pow(time,6) + c2*Math.pow(time,5) + c3*Math.pow(time,4) + 
					 c4*Math.pow(time,3) + c5*Math.pow(time,2) + c6*Math.pow(time,1)+ c7;		
		return out;
	
	}
	
	///Test Function
	public double[][] getK(){
		double[][] out =new double[r_count][theta_count];
		
		for (int theta = 1; theta<theta_count+1; theta++) {			
			for (int r = 1; r<r_count+1; r++) {
				out[r-1][theta-1]=K[r-1][theta-1];	
			}
		}
		
		return out;
		
	}
	
	//////----Successive Over-Relaxation Matrix Solver----/////////
	public double[] SOR(double[][] AA, double[] bb, double tol, double guess, double w){
		//////Uses the Successive Over Relaxation Method to Solve the system Ax=b
		//w>1 accelerates convergence but risks divergence
		//w<1 forces an unstable system to converge
		
		int n = AA.length;
		double[] xk_1 = new double[n];		
		double[] xk = new double[n];		
		
		double sigma = 0.0;
		double margin = 0.0;
		double[] temp;
		
		
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
	
	///////----Maximum percentile difference between elements----//////////
	private double maxerror(double[] x1, double[] x2){
		double out=0; 
		double diff = 0;
		
		for(int i=0; i<x1.length; i++){
			diff=100*Math.abs(((x1[i]-x2[i])/x1[i]));
			
			if(diff>out){out=diff;}//replaces lower error with greater error
			
		}
		return out;//returns greatest error
	}

	
	
}


