package cavityFETD;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.io.File;


/** Class for running the simulation solver loop
 * @author Darwin Li
 */
abstract class Solver {

	Solver() {
	}
	
	/** Runs the solver loop for the specified problem*/
	void solverLoop() {
	}
	
	/** Method for printing a double[][] for debugging into MATLAB*/
	protected void printMatrix(double[][] mat, PrintWriter out) {
		String row = "[";
		for (int i = 0; i<mat.length; i++) {
			for (int j = 0; j<mat.length; j++){
		    	if (i == mat.length-1 && j==mat.length-1) {
		    		row = row + mat[i][j] + "]";
		    		out.println(row);
	    		} else if (j == mat.length-1) {
	    			row = row + mat[i][j] + "; ";
    			} else {
    		    	row = row + mat[i][j] + " ";
    			}
		     }
		}
	}
	
    /** The maximum number of time steps*/
	protected int _maxSteps;
	
    /** The time step in seconds*/
	protected double _dt;
	
    /** The minimum frequency in Hz*/
	protected double _fMin;
	
    /** The maximum frequency in Hz*/
	protected double _fMax;
	
	
    /** The time step in seconds squared*/
	protected double _dtSquare;
	
    /** The distance between the two PEC planes*/
	protected double _plateDistance;
	
    /** The beta factor as part of Newmark Beta*/
	protected double beta;
	
    /** The E field unknowns at the current timestep */
	protected RealVector _En;
	
    /** The E field unknowns at the next timestep */
	protected RealVector _EnPlus1;
	
    /** The E field unknowns at the previous timestep */
	protected RealVector _EnMinus1;
	
    /** The f vector for the Forcing Vector. At the current timestep */
	protected RealVector _fn;
	
    /** The f vector for the Forcing Vector. At the next timestep */
	protected RealVector _fnPlus1;
	
    /** The f vector for the Forcing Vector. At the previous timestep */
	protected RealVector _fnMinus1;
	
    /** The normalized f vector for the Forcing Vector. Multiply this by the signal */
	protected RealVector f;
	
    /** ArrayList consisting of a Voltage Monitor for the source 
     * Elements are ArrayLists that contain the edges part of the monitor */
	protected Integer[] _voltageSourceMonitor;
	
    /** ArrayList consisting of a Voltage Monitor as defined by user 
     * Elements are ArrayLists that contain the edges part of the monitor */
	protected Integer[] _voltageMonitor1;
	
	/** The S Matrix as defined by Matrix Calc*/
	protected RealMatrix S;
	
	/** The T Matrix as defined by Matrix Calc*/
	protected RealMatrix T;
	
	/** The S Matrix as defined by Matrix Calc*/
	protected RealMatrix R;
	
	/** The A Matrix which combines the T, R, and S matrices with coefficients
	 * goes before the basis vector at next time step*/
	protected RealMatrix A;
	
    /** The B Matrix which combines the T and S matrices with coefficients
     * goes before the basis vector at the current time step*/
	protected RealMatrix B;
	
    /** The C Matrix which combines the T, R, and S matrices with coefficients
     * goes before the basis vector at the last time step*/
	protected RealMatrix C;
	
    /** The signal to pull the excitation from*/
	protected ExcitationSignal _signal;

    /** My output director */
	protected String _outDir;
	
    /** Verbose, for debugging.*/
	protected boolean _verbose;

}