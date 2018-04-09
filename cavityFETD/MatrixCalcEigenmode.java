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
class MatrixCalcEigenmode extends MatrixCalc {

	/** The constructor for the cavity Eignemode MatrixCalc. VERBOSE will toggle outputs to user.
	 * FILENAME designates the path to the mesh file. JPOS designates the position
	 * of the current source. VMONPOS designates the position of the voltage monitor
	 * This constructor also constructs the Matrices*/
	MatrixCalcEigenmode(boolean verbose, String fileName, PrintWriter writerIn) {
		writer = writerIn;
		meshFilePath = fileName;
		_verbose = true;
		nodeList = new ArrayList<Node>();
		eleList = new ArrayList<Element>();
		_dirichiletNodesCount = 0;
		try {
			meshFileParser(meshFilePath);
		} catch (FileNotFoundException e1) {
	        System.out.println(".node file not found");
	        writer.println(".node file not found");
			e1.printStackTrace();
		}
		fillMatrices();
	}

    /** Reads in the node file and starts processing the nodes. Also initializes the matrices, fills in the boundary locations*/
	@Override
	void meshFileParser(String e) throws FileNotFoundException {
		//String nodeStartPattern = "$Nodes";
		//String nodeEndPattern = "$Nodes";
		//String elementStartPattern = "$Elements";
		//String elementEndPattern = "$EndElements";
		//String numericPattern = "([0-9\\.]+)";
		//Pattern nodeStart = Pattern.compile(nodeStartPattern);
    	internalEdgeMap = new HashMap<ArrayList<Integer>, Integer>();
    	boundaryEdgeMap = new HashMap<ArrayList<Integer>, Integer>();
		xMin = 1.0e12;
		xMax = 0.0;
		yMin = 1.0e12;
		yMax = 0.0;
		zMin = 1.0e12;
		zMax = 0.0;
		String sNode = "([0-9]+)\\s+([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)";
		//Need to handle gmesh ascii mesh file versions 2.0 and 2.2
		String sEle22 = "[0-9]+\\s+4\\s+[0-9]+\\s+[0-9]+\\s+[0-9]+\\s+([0-9]+)\\s+([0-9]+)\\s+([0-9]+)\\s+([0-9]+)";
		String sEle20 = "[0-9]+\\s+4\\s+[0-9]+\\s+[0-9]+\\s+[0-9]+\\s+[0-9]+\\s+([0-9]+)\\s+([0-9]+)\\s+([0-9]+)\\s+([0-9]+)";
		Pattern pattNode = Pattern.compile(sNode);
		Pattern pattEle20 = Pattern.compile(sEle20);
		Pattern pattEle22 = Pattern.compile(sEle22);
		if (_verbose) {
    	   System.out.println("Looking for mesh file in " + e);
    	   writer.println("Looking for mesh file in " + e);
        }
		File fil = new File(meshFilePath);
	    Scanner scan = new Scanner(fil);
	    //Scan until we get to the nodes
	    boolean nodesUnfinished = true;
	    while(scan.hasNextLine()){
	    	String line = scan.nextLine();
	        Matcher mNode = pattNode.matcher(line);
	        Matcher mEle22 = pattEle22.matcher(line);
	        Matcher mEle20 = pattEle20.matcher(line);
	        if (nodesUnfinished &&mNode.find()) {
	        	int n = Integer.valueOf(mNode.group(1));
	        	double x = Double.valueOf(mNode.group(2));
	        	double y = Double.valueOf(mNode.group(3));
	        	double z = Double.valueOf(mNode.group(4));
	        	writer.println("Mesh file node number " + n + " internal node number " + nodeList.size());
	        	writer.println(x + " " + y + " " + z);
	        	nodeList.add(new Node(n, x, y, z, _verbose));
	        	if (x < xMin) {
	        		xMin = x;
	        	} else if (x > xMax) {
	        		xMax = x;
	        	}
	        	if (y < yMin) {
	        		yMin = y;
	        	} else if (y > yMax) {
	        		yMax = y;
	        	}
	        	if (z < zMin) {
	        		zMin = z;
	        	} else if (z > zMax) {
	        		zMax = z;
	        	}
	        // The meshfile has element lines only after all the node lines
	        } else if (line.equals("$Elements")) {
	        	for (Node i: nodeList) {
	        		_dirichiletNodesCount = _dirichiletNodesCount + i.setBound(xMin, xMax, yMin, yMax, zMin, zMax);
	        	}
	        	nodesUnfinished = false;
	        } else if (mEle20.find()) {
	        	int node0 = Integer.valueOf(mEle20.group(1))-1;
	        	int node1 = Integer.valueOf(mEle20.group(2))-1;
	        	int node2 = Integer.valueOf(mEle20.group(3))-1;
	        	int node3 = Integer.valueOf(mEle20.group(4))-1;
	        	Element tempEle = new Element(eleList.size(), nodeList.get(node0), nodeList.get(node1), nodeList.get(node2), nodeList.get(node3), _verbose, writer);
	        	eleList.add(tempEle);
	        	tempEle.mapEdges(internalEdgeMap, boundaryEdgeMap);
	        } else if(mEle22.find()) {
	        	int node0 = Integer.valueOf(mEle22.group(1))-1;
	        	int node1 = Integer.valueOf(mEle22.group(2))-1;
	        	int node2 = Integer.valueOf(mEle22.group(3))-1;
	        	int node3 = Integer.valueOf(mEle22.group(4))-1;
	        	Element tempEle = new Element(eleList.size(), nodeList.get(node0), nodeList.get(node1), nodeList.get(node2), nodeList.get(node3), _verbose, writer);
	        	eleList.add(tempEle);
	        	tempEle.mapEdges(internalEdgeMap, boundaryEdgeMap);
	        } else if(_verbose) {
        		System.out.println("Skipped line with " + line);
        		writer.println("Skipped line with " + line);
	        }
	    }
	    scan.close();
	}
}