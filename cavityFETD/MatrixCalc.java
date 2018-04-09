package cavityFETD;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;

import java.util.HashMap;
import java.util.HashSet;
import java.io.PrintWriter;
import java.lang.Integer;
import java.lang.Double;
import java.util.ArrayList;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** The Matrix calc computations are initiated here.
 * @author Darwin Li
 */
abstract class MatrixCalc {

	/** Abstract class for Matrix Calculations*/
	MatrixCalc() {
		nodeList = new ArrayList<Node>();
		eleList = new ArrayList<Element>();
    	internalEdgeMap = new HashMap<ArrayList<Integer>, Integer>();
    	boundaryEdgeMap = new HashMap<ArrayList<Integer>, Integer>();
    	backwardsEdgeSet = new HashSet<Integer>();
    	_dirichiletNodesCount = 0;
	}

    /** Reads in the node file and starts processing the nodes. Also initializes the matrices, fills in the boundary locations*/
	void meshFileParser(String e) throws FileNotFoundException {
	}
	
	/** Fills the T and S Matrices and fills the forcing vector F**/
	protected void fillMatrices(){
		if (_verbose) {
			writer.println("boundary start");
			for (ArrayList<Integer> test: boundaryEdgeMap.keySet()){
	            String key = test.get(0) + "," + test.get(1);
	            String value = boundaryEdgeMap.get(test).toString();  
	            writer.println(key + " " + value);  
	        }
			writer.println("internal start");
			for (ArrayList<Integer> test: internalEdgeMap.keySet()){
	            String key = test.get(0) + "," + test.get(1);
	            String value = internalEdgeMap.get(test).toString();  
	            writer.println(key + " " + value);  
	        }
		}
    	int n = internalEdgeMap.size();
    	tMat = MatrixUtils.createRealMatrix(n, n);
    	sMat = MatrixUtils.createRealMatrix(n, n);
    	rMat = MatrixUtils.createRealMatrix(n, n);
    	voltageSourceMonitor = new HashSet<Integer>();
    	voltageMonitor1 = new HashSet<Integer>();
    	double[] empty = new double[n];
		fVec = MatrixUtils.createRealVector(empty);
		for (Element e: eleList) {
			e.addMatrixTerms(internalEdgeMap, boundaryEdgeMap, rMat, tMat, 
					sMat, fVec, voltageSourceMonitor, voltageMonitor1, backwardsEdgeSet);
		}
	}
	
	/** Returns my complete list of nodes.*/
	ArrayList<Node> getNodeList() {
		return nodeList;
	}
	
	
	/** Returns the list of edges that should be monitored for the source's 
	 * Voltage Monitor*/
	HashSet<Integer> getSourceVoltageMonitor() {
		return voltageSourceMonitor;
	}
	
	/** Returns the list of edges that should be monitored for the Voltage Monitor*/
	HashSet<Integer> getVoltageMonitor1() {
		return voltageMonitor1;
	}
	
	/** Returns the T Matrix*/
	RealMatrix getT() {
		return tMat;
	}
	
	/** Returns the S Matrix*/
	RealMatrix getS() {
		return sMat;
	}
	
	/** Returns the R Matrix*/
	RealMatrix getR() {
		return rMat;
	}
	
	/** Returns the forcing vector F*/
	RealVector getF() {
		return fVec;
	}
	
	/** Returns the number of Dirichilet Boundary Nodes */
	int getDirichiletNodesCount() {
		return _dirichiletNodesCount;
	}
	
	
    /** The distance between the two PEC planes*/
	double getPlateDistance() {
		return _plateDistance;
	}
	
	HashSet<Integer> getBackwardsSet(){
		return backwardsEdgeSet;
	}
	
    /** The min max X, Y, Z boundaries of the domain*/
	double getXMin() {
		return xMin;
	}
	double getXMax() {
		return xMax;
	}
	double getYMin() {
		return yMin;
	}
	double getYMax() {
		return yMax;
	}
	double getZMin() {
		return zMin;
	}
	double getZMax() {
		return zMax;
	}
	
    /** Returns the boolean for verbose.*/
	boolean getVerbose() {
		return _verbose;
	}
	
    /** Returns the PrintWriter for the logfile.*/
	PrintWriter getWriter() {
		return writer;
	}

	
    /** The List of Nodes in the Computational Domain.*/
	protected ArrayList<Node> nodeList;
	
    /** The List of Elements in the Computational Domain.*/
	protected ArrayList<Element> eleList;
	
    /** The List of Edges to record data in the Computational Domain.
     * These edges lie on the current source*/
	protected HashSet<Integer> voltageSourceMonitor;
	
    /** The List of Edges to record data in the Computational Domain.
     * These edges lie on the user specified Voltage Monitor*/
	protected HashSet<Integer> voltageMonitor1;
	
	/** The number of Dirichilet nodes ie. nodes thatlie on a Dirichilet Boundary. Neumann
	 * Boundaries are still counted as a free node.
	 */
	protected int _dirichiletNodesCount;
	
	/** The X,Y coordinate of my current Source as a user input in an array of 2 items*/
	protected double jPosition[];
	
	/** The X,Y coordinate of my current Source mapped to the mesh in an array of 2 items*/
	protected double jPositionMapped[];
	
	/** The X,Y coordinate of my voltageMonitor1 in an array of 2 items*/
	protected double vMonitorPosition[];
	
	
	/** The X,Y coordinate of my voltage monitor mapped to the mesh in an array of 2 items*/
	protected double vMonitorPositionMapped[];
	
	/** My +-x,y,z boundaries*/
	protected double xMin;
	protected double yMin;
	protected double zMin;
	protected double xMax;
	protected double yMax;
	protected double zMax;
	
    /** The distance between the two PEC planes*/
	protected double _plateDistance;
	
    /** The mapping of Edges from a pair of nodes. The pair of nodes
     * will always be in increasing global order, ie. there will never be
     * an edge with nodes {2, 1}, it must be {1, 2}. Neumann boundaries
     * edges and sources are included in this map among internal edges*/
	protected HashMap<ArrayList<Integer>, Integer> internalEdgeMap;
	
    /** The mapping of Edges from a pair of nodes. The pair of nodes
     * will always be in increasing global order, ie. there will never be
     * an edge with nodes {2, 1}, it must be {1, 2}. This handles only
     * the Dirichilet boundary case for now*/
	protected HashMap<ArrayList<Integer>, Integer> boundaryEdgeMap;
	
	/** Handles the case where an Edge is backwards along a Z Axis if we want edges to always
	 * point in direction of lower Z Value to Higher z value along its two nodes
	 */
	protected HashSet<Integer> backwardsEdgeSet;
	
    /** The forcing vector for the dirichilet boundaries and source.*/
	protected RealVector fVec;
	
    /** The S matrix. Has the dotted Curl integrals by mu*/
	protected RealMatrix sMat;
	
	/** The T matrix. Has the dotted integral times epsilon*/
	protected RealMatrix tMat;
	
	/** The R matrix. Has the dotted integral times sigma*/
	protected RealMatrix rMat;
	
    /** My PrintWriter*/
	protected PrintWriter writer;
	
    /** Verbose, for debugging.*/
	protected boolean _verbose;
	
	/** The mesh file I need.*/
	protected String meshFilePath;
	
}
