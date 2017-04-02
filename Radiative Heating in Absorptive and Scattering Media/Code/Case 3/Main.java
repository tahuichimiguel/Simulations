
import java.io.*;
import java.lang.Math;

class Main{


	public static void main(String[] args){

		double thickness = .01;//in meters
		
		//int nodecount = 48;  //Coarse Mesh
		//int nodecount = 98;  //Normal Mesh
		//int nodecount = 152; //Fine Mesh
		//int nodecount = 198; //Finer Mesh
		
		int nodecount = 500; //Finest Mesh
		
		double kappa=0;
		
		//double sigma=10;
		
		for(double i = 100;i<1010;i+=100){

		double sigma=i;
		
		double[] mesh;
		double[] heatflux = new double[nodecount];
		
		//Spectral Intensities
		double[] IMinus;
		double[] IPlus;
		double[] I;
		
		String filename;
		
		MeshGen generator;
		Solver solver = new Solver(nodecount);
		Writer pen = new Writer();
		
		//Specify the Mesh Dimensions and Number of Computational Nodes
		generator = new MeshGen(0,thickness,nodecount);
		generator.generate();
		mesh = generator.getMesh();
		
		//Specify the Number of Computational Nodes
		solver = new Solver(nodecount);
		
		//Input Beta and the slab thickness
		solver.DiffEqConstants(sigma, kappa ,thickness);
		
		System.out.println("Initializing I+ and I-");
		
		solver.InitializeIPlusIMinus(1000);//Max Iterations = 100
		//solver.AssembleMatrices();
		
		//System.out.println("");
		
		//System.out.println("");
		IPlus = solver.getIPlus();
		IMinus = solver.getIMinus();
		//I = solver.SolveI();
		
		
		 //System.out.println("\n Heat Flux \n");
		 for(int ii = 0; ii<nodecount; ii++){
			 heatflux[ii]=(Math.PI*(IPlus[ii]-IMinus[ii]));
		 }
		
		double[] divheat = new double[nodecount];
		//System.out.println("\n Heat Flux \n");
		for(int ii = 1; ii<nodecount-2; ii++){
			divheat[ii]=(heatflux[ii+1]-heatflux[ii-1])/(.01/(2*(nodecount+2)));
		}
		
		//filename = new String("Intensity -> [S] = " +  sigma  );
		//pen.writeVector(I,filename);
		
		filename = new String("I+ -> [S] = " + sigma );
		pen.writeVector(IPlus,filename);

		filename = new String("I- -> [S] = " + sigma );
		pen.writeVector(IMinus,filename);
		
		filename = new String("Heat Flux -> [S] = " + sigma );
		pen.writeVector(heatflux,filename);
			
		filename = new String("DivHeat -> [K,S] = " + "[" + kappa + " , " + sigma + "]");
		pen.writeVector(divheat,filename);
			
		//pen.writeVector(mesh, new String("mesh"));
		}
	}



}