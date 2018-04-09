package cavityFETD;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;








/** Run the Simulation Solver loop for an Eigenmode Analysis
 * @author Darwin Li
 */
class SolverEigenmode extends Solver {

	SolverEigenmode(double fMin, double fMax, MatrixCalc matCalc, String outDir) {
		if (matCalc.getVerbose()) {
			_verbose = true;
		}
		_outDir = outDir;
		_fMin = fMin;
		_fMax = fMax;
		S = matCalc.getS().scalarMultiply(FastMath.pow(PhysicsConstants.u_r, -1));
		T = matCalc.getT().scalarMultiply(PhysicsConstants.e_r);
		_numFreeNodes = matCalc.getNodeList().size() - matCalc.getDirichiletNodesCount();
	}
	
	/** Run the Solver loop for the Cavity Eigenmode analysis problem */
	@Override
	void solverLoop() {
		PrintWriter writer = null;
    	try {
			writer = new PrintWriter(_outDir + "\\" + "ResonanceFrequencies.txt", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		LUDecomposition solverT = new LUDecomposition(T);
		DecompositionSolver decomp = solverT.getSolver();
		RealMatrix TS = decomp.solve(S);
		EigenDecomposition eigDecomp = new EigenDecomposition(TS);
		double[] eigValues = eigDecomp.getRealEigenvalues();
		Arrays.sort(eigValues);
		eigValues = Arrays.copyOfRange(eigValues, _numFreeNodes, eigValues.length);
		for (int i = 0; i < eigValues.length; i ++) {
			double evNonAngular = FastMath.sqrt(eigValues[i])*(FastMath.sqrt(PhysicsConstants.u_f*PhysicsConstants.e_f))/(2*FastMath.PI);
			double kc = FastMath.sqrt(eigValues[i]);
			if (evNonAngular >= _fMin) {
				System.out.println("Resonance at " + evNonAngular + " Hz");
				writer.println("Resonance at " + evNonAngular + " Hz");
				System.out.println(kc + " rad/m");
				writer.println(kc + " rad/m");
			}
			if (evNonAngular > _fMax) {
				break;
			}
		}
		writer.close();
	}
	
	private void printMATLABfile(double[][] Smat, double[][] Tmat) {
		PrintWriter writer = null;
    	try {
			writer = new PrintWriter(_outDir + "\\" + "RunMeInMatLabForEigenValues.m", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		String row = "[";
		for (int i = 0; i<Smat.length; i++) {
			for (int j = 0; j<Smat.length; j++){
		    	if (i == Smat.length-1 && j==Smat.length-1) {
		    		row = row + Smat[i][j] + "];";
		    		writer.println("S = " + row);
	    		} else if (j == Smat.length-1) {
	    			row = row + Smat[i][j] + "; ";
    			} else {
    		    	row = row + Smat[i][j] + " ";
    			}
		     }
		}
		row = "[";
		for (int i = 0; i<Tmat.length; i++) {
			for (int j = 0; j<Tmat.length; j++){
		    	if (i == Tmat.length-1 && j==Tmat.length-1) {
		    		row = row + Tmat[i][j] + "];";
		    		writer.println("T = " + row);
	    		} else if (j == Tmat.length-1) {
	    			row = row + Tmat[i][j] + "; ";
    			} else {
    		    	row = row + Tmat[i][j] + " ";
    			}
		     }
		}
		writer.println("[eigVectors, eigValues] = + eig(S, T)");
		writer.println("eigValues = diag(eigValues);");
		writer.println("eigValues = sort(eigValues);");
		writer.println("eigValues(" + _numFreeNodes + "+1,:)");
	}
	
	/** Computes the first 8 analytical eigenvalues based on a rectangular cavity. Need to fix this, but a
	 * general iteration for the first 8 modes is not that easy*/
	private void getAnalyticalEigenValues(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax) {
		double[] tempTE = new double[8];
		double[] tempTM = new double[8];
		double[] tempTETM = new double[16];
		int n = 0;
		int m = 0;
		int l = 0;
		int count = 0;
		int countTETM = 0;
		while (countTETM < 16) {
			double dTE = FastMath.sqrt(FastMath.pow(n*FastMath.PI/(xMax-xMin), 2) + FastMath.pow(m*FastMath.PI/(yMax-yMin), 2) + FastMath.pow(l*FastMath.PI/(zMax-zMin), 2));
			double dTM = FastMath.sqrt(FastMath.pow(l*FastMath.PI/(xMax-xMin), 2) + FastMath.pow(m*FastMath.PI/(yMax-yMin), 2) + FastMath.pow(n*FastMath.PI/(zMax-zMin), 2));
			tempTE[count] = dTE;
			tempTM[count] = dTM;
			count++;
			tempTETM[countTETM] = dTE;
			tempTETM[countTETM + 1] = dTM;
			count = count + 2;
		}
		Arrays.sort(tempTETM);
		correctAnalyticalValues = Arrays.copyOfRange(tempTETM, 0, 9);
	}
	
	/** The count of the number of free nodes I have. */
	private int _numFreeNodes;
	
	/** The correct EigenValues */
	private double[] correctAnalyticalValues;
}