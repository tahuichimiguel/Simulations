
import java.io.*;
import java.lang.Math;

class Main{


	public static void main(String[] args){
		//int radd=10; //-- 641 cells for circle
		//int radd=15; //-- 961
		//int radd=20; //-- 1281
		//int radd=25; //-- 1501
		int radd = 30;
		
		//int radd = 50;
		//int thetad = 16;
		
		int thetad=16;
		double h = 2;

		//int t_steps = 700; //0.9 seconds
		//int t_steps = 800; //0.78 seconds
		//int t_steps = 900; //0.7 seconds
		//int t_steps = 1000; //0.63 seconds
		int t_steps = 1100; //0.63 seconds
		
		//for(int radd = 10; radd<26; radd +=5){
		
		//thetad MUST be integral multiples of 8
		
		Solver2 solver = new Solver2(h,radd,thetad,t_steps);
		//Solver2 solver = new Solver2(10,radd,thetad,400);
		//double[][] A = new double[4][4];
		//double[] b = new double[4];
		double[] out;
		
		solver.ConstructConductivity();
		solver.AssembleStiffnessMatrix();
		solver.Initialize();
		solver.Solve();
		
		//out=solver.getT(t_steps-1);
		//out= solver.getK();
		
		//System.out.println(	"Solution");	
		//for (int i = 0; i<radd*thetad; i++) {
		//	System.out.println(	out[i]);
		//}
		
		solver.writetoFile(275,h,radd,thetad);
		solver.writetoFile(550,h,radd,thetad);
		solver.writetoFile(825,h,radd,thetad);
		solver.writetoFile(t_steps-1,h,radd,thetad);
		/*for (int theta =0; theta<thetad; theta++) {
			
			for (int r = 0; r<radd; r++) {
				System.out.println(out[r][theta]);
				
			}
		}*/
		//}
	}



}