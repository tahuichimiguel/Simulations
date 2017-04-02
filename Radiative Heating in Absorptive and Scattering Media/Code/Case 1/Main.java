
import java.io.*;
import java.lang.Math;

class Main{


	public static void main(String[] args){

		double thickness = .01;//in meters
		String filename;
		
		//int nodecount = 48;  //Coarse Mesh
		//int nodecount = 98;  //Normal Mesh
			
		//int nodecount = 152; //Fine Mesh
		int nodecount = 198; //Finest Mesh

		//-----------------------------------//-
		
		double wavelength;//in microns
		double kappa;
		double beta;
		
		//beta =  Math.pow(10,-19); //in m^-1
		//kappa = Math.pow(10,-19);//in m^-1
		
		//beta =  15; //in m^-1		//for band 2
		//kappa = 15;//in m^-1		//for band 2
		
		//beta =  600; //in m^-1	//for band 3
		//kappa = 600;//in m^-1		//for band 3
		
		beta =  10000; //in m^-1	//for band 4
		kappa = 10000;//in m^-1	//for band 4
		
		//for(double i = 0.1; i<0.5;i+=0.1){
		//for(int i = 0; i<10; i++){ //for band 2
		//for(int i = 0; i<9; i++){ //for band 3
		for(int i = 5; i<550; i=i+1){ //for band 4	
		
		//wavelength = i; //band 1	
		//wavelength = 0.8+0.2*i; //band 2
		//wavelength = 2.8+0.2*i; //band 3
		wavelength = i; //band 4
			
		double[] mesh;
		
		//Spectral Intensities
		double[] IMinus;
		double[] IPlus;
		//double[] I;
		
		MeshGen generator;
		Solver solver = new Solver(nodecount);
		Radiation rad = new Radiation();
		Writer pen = new Writer();
		
		//Specify the Mesh Dimensions and Number of Computational Nodes
		generator = new MeshGen(0,thickness,nodecount);
		generator.generate();
		mesh = generator.getMesh();
		
		//Specify the Number of Computational Nodes
		solver = new Solver(nodecount);
		
		//Input Beta and the slab thickness
		solver.DiffEqConstants(beta,thickness);
		
		//Input the wavelength of interest and the absorption coefficient
		solver.BlackBodyConstants(wavelength,kappa);
		
		solver.AssembleMatrices();
		//System.out.println("");
		
		//System.out.println("");
		IPlus = solver.SolveIPlus();
		IMinus = solver.SolveIMinus();
		//I = solver.SolveI();
		
		
		//filename = new String(wavelength +  " Intensity");
		//pen.writeVector(I,filename);
		
		filename = new String(wavelength +  " I+");
		pen.writeVector(IPlus,filename);

		filename = new String(wavelength +  " I-");
		pen.writeVector(IMinus,filename);
		
		 
		}//End of Wavelength Loop
		
		//pen.writeVector(mesh, new String("mesh"));
		
		
		
		
	}



}