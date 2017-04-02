import java.lang.Math;

class MeshGen{
	
	/////////
	/////////
	////NOTE: This class makes only the mesh that will be used during post-processing. 
	/////////  The Solver class assembles and solves the A-Matrix.
	/////////
	
	//////NOTE: This mesh generator makes ONLY uniform linear meshes
	
	double xmin, xmax;
	int node_count;
	
	double[] mesh;
	double[] positions; //The spatial locations of the nodes defined by xmin, xmax, and 
	
	///NOTE: The dimension of the solution vector is equal to the number of nodes(node_count)
	///NOTE: The A-matrix is a square matrix with rank = node_count

	
	public MeshGen(double x1, double x2, int nodes){
		
		xmin = x1;
		xmax = x2;

		node_count = nodes+2;
		
		////Initializes mesh with node_count rows
		////mesh[][1] = positions
		////mesh[][2] = empty matrix elements for post-processing use (storage of the solutions)
		mesh = new double[node_count];
	}

	////Generates positions for all nodes and stores them in the mesh vector
	public void generate(){
		mesh[0]=0;
		
			for (int i = 0; i<node_count-2; i++) {
			mesh[i+1] = xmin+(i+1)*(xmax-xmin)/(node_count);
		}
		mesh[node_count-1]=xmax;
	}

	public double[] getMesh(){return mesh;}//returns the vector with all the 
	
}