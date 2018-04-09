package cavityFETD;

import cavityFETD.PhysicsConstants;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.io.PrintWriter;

/** A Tetrahedron is an element with 4 points
 * Each Tetrahedron is defined locally in the class
 *  @author Darwin Li
 */
public class Element {
	
	/** The constructor that designates a number N
	 * as a name and sets the nodes attached to the 
	 * element as A B C and D. */
    public Element(int n, Node a, Node b, Node c, Node d, boolean verbose, PrintWriter writeR) {
        _verbose = verbose;
        _n = n;
        writer = writeR;
        nodes = nodeOrder(a, b, c, d);
        _a = nodes[0].getNumber();
        _b = nodes[1].getNumber();
        _c = nodes[2].getNumber();
        _d = nodes[3].getNumber();
        RealMatrix volMat = MatrixUtils.createRealMatrix(4, 4);
        coeffMat = MatrixUtils.createRealMatrix(4,4);
        for (int i = 0; i < 4; i ++) {
        	for (int j = 0; j < 4; j++) {
            	if (i == 0) {
                   	coeffMat.addToEntry(i, j, nodes[j].getX());
                   	volMat.addToEntry(j, i, 1.0);
            	} else if (i == 1) {
                   	coeffMat.addToEntry(i, j, nodes[j].getY());
            		volMat.addToEntry(j, i, nodes[j].getX());
            	} else if (i == 2) {
                   	coeffMat.addToEntry(i, j, nodes[j].getZ());
            		volMat.addToEntry(j, i, nodes[j].getY());
            	} else if (i == 3) {
                   	coeffMat.addToEntry(i, j, 1.0);
            		volMat.addToEntry(j, i, nodes[j].getZ());
            	}
        	}
        }
        coeffMat = MatrixUtils.inverse(coeffMat);
        _volume = new LUDecomposition(volMat).getDeterminant()/6;
        _volume = FastMath.abs(_volume);
 	   if (_verbose) {
    	   System.out.println(this.getName()+" has Nodes " 
 	   + _a + "," + _b + "," + _c + "," + _d + " and volume of " + _volume);
    	   writer.println(this.getName()+" has Nodes " 
 	   + _a + "," + _b + "," + _c + "," + _d + " and volume of " + _volume);
    	   writer.println(this.getName()+" has coefficient matrix of");
    	   writer.println(coeffMat);
 	   }
    }

    /** Return my name.*/
    int getName() {
        return _n;
    }
    
    /** Return my first node name.*/
    int getA() {
        return _a;
    }
    
    /** Return my second node name.*/
    int getB() {
        return _b;
    }
    
    /** Return my third node name.*/
    int getC() {
        return _c;
    }
    
    /** Return my fourth node name.*/
    int getD() {
        return _d;
    }
    
    /** Return the local edge number. The convention used is
     * Edge 1 is between first and second node, Edge 2
     * is between first and third. Edge 3 is between
     * first and fourth. Edge 4 is between second and third
     * Edge 5 is between second and fourth. Edge 6 is
     * between third and fourth. Index starts at 0, so
     * edge 1 is index 0. */
    int getEdge(int i, int j) {
    	int temp;
    	temp = 100;
        if (i == 0 && j == 1) {
        	temp = 0;
        } else if (i ==0 && j == 2) {
        	temp = 1;
        } else if (i ==0 && j == 3) {
        	temp = 2;
        } else if (i ==1 && j == 2) {
        	temp = 3;
        } else if (i ==1 && j == 3) {
        	temp = 4;
        } else if (i ==2 && j == 3) {
        	temp = 5;
        }
        return temp;
    }
    
    /** Return local index edge E's pair of Nodes in an array with
     * the convention defined as: 
     * Edge 1 is between first and second node, Edge 2
     * is between first and third. Edge 3 is between
     * first and fourth. Edge 4 is between second and third
     * Edge 5 is between second and fourth. Edge 6 is
     * between third and fourth. Index starts at 0, so
     * edge 1 is index 0.*/
    private Node[] getNodesOfEdge(int e) {
        if (e == 0) {
        	Node[] temp = {nodes[0], nodes[1]};
        	return temp;
        } else if (e == 1) {
        	Node[] temp = {nodes[0], nodes[2]};
        	return temp;
        } else if (e == 2) {
          	Node[] temp = {nodes[0], nodes[3]};
          	return temp;
        } else if (e == 3) {
         	Node[] temp = {nodes[1], nodes[2]};
         	return temp;
        } else if (e == 4) {
        	Node[] temp = {nodes[1], nodes[3]};
        	return temp;
        } else if (e == 5) {
        	Node[] temp = {nodes[2], nodes[3]};
        	return temp;
        }
        return new Node[1];
    }
    
    /** Return local index edge E's pair of Node's global indices in an array with
     * the convention defined as: 
     * Edge 1 is between first and second node, Edge 2
     * is between first and third. Edge 3 is between
     * first and fourth. Edge 4 is between second and third
     * Edge 5 is between second and fourth. Edge 6 is
     * between third and fourth. Index starts at 1, so
     * edge 1 is index 1, This is because of .getNumber() */
    private ArrayList<Integer> getGlobalNodeNumbersOfEdge(int e) {
    	ArrayList<Integer> tempArr = new ArrayList<Integer>();
        if (e == 0) {
        	tempArr.add(nodes[0].getNumber());
        	tempArr.add(nodes[1].getNumber());
        } else if (e == 1) {
        	tempArr.add(nodes[0].getNumber());
        	tempArr.add(nodes[2].getNumber());
        } else if (e == 2) {
        	tempArr.add(nodes[0].getNumber());
        	tempArr.add(nodes[3].getNumber());
        } else if (e == 3) {
        	tempArr.add(nodes[1].getNumber());
        	tempArr.add(nodes[2].getNumber());
        } else if (e == 4) {
        	tempArr.add(nodes[1].getNumber());
        	tempArr.add(nodes[3].getNumber());
        } else if (e == 5) {
        	tempArr.add(nodes[2].getNumber());
        	tempArr.add(nodes[3].getNumber());
        }
        return tempArr;
    }
    
    /** Return E's pair of Node's local indices in an array with
     * the convention defined as: 
     * Edge 1 is between first and second node, Edge 2
     * is between first and third. Edge 3 is between
     * first and fourth. Edge 4 is between second and third
     * Edge 5 is between second and fourth. Edge 6 is
     * between third and fourth. Index starts at 0, so
     * edge 1 is index 0.*/
    private int[] getLocalNodeNumbersOfEdge(int e) {
        if (e == 0) {
        	int[] temp = {0, 1};
        	return temp;
        } else if (e == 1) {
        	int[] temp = {0, 2};
        	return temp;
        } else if (e == 2) {
          	int[] temp = {0, 3};
          	return temp;
        } else if (e == 3) {
         	int[] temp = {1, 2};
         	return temp;
        } else if (e == 4) {
        	int[] temp = {1, 3};
        	return temp;
        } else if (e == 5) {
        	int[] temp = {2, 3};
        	return temp;
        }
        return new int[1];
    }
    
    /** Return an array of nodes ordered by their global
     * numbering scheme from smallest to largest*/
    private Node[] nodeOrder(Node aa, Node bb, Node cc, Node dd) {
    	Node low1;
    	Node low2;
    	Node high1;
    	Node high2;
    	Node lowest;
    	Node middle1;
    	Node middle2;
    	Node highest;
    	if (aa.getNumber() < bb.getNumber()) {
    		low1 = aa;
    		high1 = bb;
    	} else {
    		low1 = bb;
    		high1 = aa;
    	}
    	if (cc.getNumber() < dd.getNumber()) {
    		low2 = cc;
    		high2 = dd;
    	} else {
    		low2 = dd;
    		high2 = cc;
    	}
    	if (low1.getNumber() < low2.getNumber()) {
    		lowest = low1;
    		middle1 = low2;
    	} else {
    		lowest = low2;
    		middle1 = low1;
    	}
    	if (high1.getNumber() > high2.getNumber()) {
    		highest = high1;
    		middle2 = high2;
    	} else {
    		highest = high2;
    		middle2 = high1;
    	}
    	if ( middle1.getNumber() < middle2.getNumber()) {
    		return new Node[] {lowest, middle1, middle2, highest};
    	}
        return new Node[]{lowest, middle2, middle1, highest};
    }
    
    /** Computes the phi term between nodes locally numbered I and J. Phi is defined as
     * the dot product of grad(lambda_i) and grad(lambda_j)*/
    private double phi(int i, int j) {
    	double bi = coeffMat.getEntry(i, 0);
    	double bj = coeffMat.getEntry(j, 0);
    	double ci = coeffMat.getEntry(i, 1);
    	double cj = coeffMat.getEntry(j, 1);
    	double di = coeffMat.getEntry(i, 2);
    	double dj = coeffMat.getEntry(j, 2);
        return bi*bj + ci*cj + di*dj;
    }
    
    /** Computes the RealVector v defined on the edge E with end nodes numbered I and J.
     * v is defined as the cross product between grad(lambda_i) and grad(lambda_j)*/
    private RealVector vVector(int e) {
    	int i = getLocalNodeNumbersOfEdge(e)[0];
    	int j = getLocalNodeNumbersOfEdge(e)[1];
    	double bi = coeffMat.getEntry(i, 0);
    	double bj = coeffMat.getEntry(j, 0);
    	double ci = coeffMat.getEntry(i, 1);
    	double cj = coeffMat.getEntry(j, 1);
    	double di = coeffMat.getEntry(i, 2);
    	double dj = coeffMat.getEntry(j, 2);
    	double[] tempArr = new double[3];
    	tempArr[0] = ci*dj - cj*di;
    	tempArr[1] = bj*di - bi*dj;
    	tempArr[2] = bi*cj - bj*ci;
        return MatrixUtils.createRealVector(tempArr);
    }
    
    /** Computes the S matrix term between edges numbered I and J*/
    private double sIJ(int eI, int eJ) {
    	RealVector vI = vVector(eI);
    	RealVector vJ = vVector(eJ);
    	return 4*getAbsVolume()*vI.dotProduct(vJ);
    }
    
    /** Computes the T matrix term between edges numbered I and J*/
    private double tIJ(int eI, int eJ) {
    	int i1 = getLocalNodeNumbersOfEdge(eI)[0];
    	int i2 = getLocalNodeNumbersOfEdge(eI)[1];
    	int j1 = getLocalNodeNumbersOfEdge(eJ)[0];
    	int j2 = getLocalNodeNumbersOfEdge(eJ)[1];
    	double phiI1J1 = phi(i1, j1);
    	double phiI1J2 = phi(i1, j2);
    	double phiI2J1 = phi(i2, j1);
    	double phiI2J2 = phi(i2, j2);
    	double MI1J1=1.0;
    	double MI1J2=1.0;
    	double MI2J1=1.0;
    	double MI2J2=1.0;
    	if (i1 == j1) {
    		MI1J1 = 2.0;
    	} else if (i1 == j2) {
    		MI1J2 = 2.0;
    	}
    	// should never have the case where i2==j1==i1==j2
    	if (i2 == j1) {
    		MI2J1 = 2.0;
    	} else if (i2 == j2) {
    		MI2J2 = 2.0;
    	}
    	return getAbsVolume()*(phiI2J2*MI1J1-phiI2J1*MI1J2-phiI1J2*MI2J1+phiI1J1*MI2J2)/20;
    }
    
    /** Maps my edges to the global HashMap of internal edges or dirichilet boundary edges. If the two nodes
     * of one of my edges is not already a key, then make a new entry mapping to a new integer globally. */
    void mapEdges(HashMap<ArrayList<Integer>, Integer> iMap, HashMap<ArrayList<Integer>, Integer> dMap) {
    	for (int i = 0; i < 6; i ++) {
    		ArrayList<Integer> tempArr= getGlobalNodeNumbersOfEdge(i);
    		if (getNodesOfEdge(i)[0].getBound().equals(getNodesOfEdge(i)[1]) && !getNodesOfEdge(i)[0].getBound().equals("none")) {
        		if (dMap.isEmpty()) {
        			dMap.put(tempArr, 0);
        		} else {
        			if (!dMap.containsKey(getGlobalNodeNumbersOfEdge(i))) {
                    	dMap.put(tempArr, dMap.size());
        			}
        		}
    		} else {
        		if (iMap.isEmpty()) {
        			iMap.put(tempArr, 0);
        		} else {
        			if (!iMap.containsKey(tempArr)) {
                    	iMap.put(tempArr, iMap.size());
        			}
        		}
    		}
    	}
    }
 
    
    /** Return my absolute Volume in meters^3.*/
    double getAbsVolume() {
        return FastMath.abs(_volume);
    }
    
    /** Contribute to the T and S matrices using IMAP to find the Global Index of the interior
     * or Neumann boundary edge. Also adds the source factor to the forcing vector F. The R matrix
     * entries are just the T Matrix entries multiplied by a ratio of epsilon and sigma since there
     * is no handling of Surface Admittances */
    void addMatrixTerms(HashMap<ArrayList<Integer>, Integer> iMap, HashMap<ArrayList<Integer>, Integer> dMap, RealMatrix r, 
    		RealMatrix t, RealMatrix s, RealVector f, HashSet<Integer> voltageSourceMonitor, 
    		HashSet<Integer> voltageMonitor1, HashSet<Integer> backwardsEdgeSet) {
    	for (int e1 = 0; e1 < 6; e1++) {
    		ArrayList<Integer> tempArr = getGlobalNodeNumbersOfEdge(e1);
    		Integer eI = iMap.get(tempArr);
    		if (eI!=null) {
        		for (int e2 = 0; e2 < 6; e2++) {
        			ArrayList<Integer> tempArr2 = getGlobalNodeNumbersOfEdge(e2);
            		Integer eJ = iMap.get(tempArr2);
    				HashSet<Integer> debugger = new HashSet<Integer>();
    				for (Integer i = 1; i < 17; i ++) {
    					debugger.add(i);
    				}
    				if (debugger.contains(tempArr2.get(0)) && debugger.contains(tempArr2.get(1))) {
    					System.out.println("Look at edge " + eI + ", which has nodes " + tempArr2.get(0) + " " + tempArr2.get(1));
    				}
        			if (eJ!=null) {
        				double tempIJ = tIJ(e1,e2);
        				double tempSIJ = sIJ(e1, e2);
        				s.addToEntry(eI, eJ, tempSIJ);
        				t.addToEntry(eI, eJ, tempIJ);
        				r.addToEntry(eI, eJ, tempIJ);
        			}
        		}
        		Node[] tempNodeArr = getNodesOfEdge(e1);
        		if (tempNodeArr[0].getMonitor()==3 && tempNodeArr[1].getMonitor()==3) {
        			double toAdd = 0.25*getAbsVolume();
        			if (tempNodeArr[0].getZ() > tempNodeArr[1].getZ()) {
        				backwardsEdgeSet.add(eI);
        			}
        			f.addToEntry(eI, 0.25*getAbsVolume());
        			voltageSourceMonitor.add((Integer) eI);
        			if (_verbose) {
        				System.out.println("Element " + _n + " has source nodes " + getNodesOfEdge(e1)[0].getNumber() + " " + getNodesOfEdge(e1)[1].getNumber());
        				writer.println("Element " + _n + " has source nodes " + getNodesOfEdge(e1)[0].getNumber() + " " + getNodesOfEdge(e1)[1].getNumber());
        			}
        		} else if (tempNodeArr[0].getMonitor()==4 && tempNodeArr[1].getMonitor()==4) {
        			voltageMonitor1.add((Integer) eI);
        			if (tempNodeArr[0].getZ() > tempNodeArr[1].getZ()) {
        				backwardsEdgeSet.add(eI);
        			}
        		}
    		} else {
    			if (_verbose) {
    				if (dMap.get(tempArr) == null) {
    					System.out.println("Uh oh. The nodes " + tempArr.get(0) + ", " + tempArr.get(1) + 
    							" are mapped to Z locations " + getNodesOfEdge(e1)[0].getZ() + ", " + getNodesOfEdge(e1)[1].getZ());
    					writer.println("Uh oh. The nodes " + tempArr.get(0) + ", " + tempArr.get(1) + 
    							" are mapped to Z locations " + getNodesOfEdge(e1)[0].getZ() + ", " + getNodesOfEdge(e1)[1].getZ());
    				}
    				System.out.println("PEC found on edge " + dMap.get(tempArr) + " connecting nodes " + tempArr.get(0) + ", " + tempArr.get(1));
    				writer.println("PEC found on edge " + dMap.get(tempArr) + " connecting nodes " + tempArr.get(0) + ", " + tempArr.get(1));
    			}
    		}
    	}    	       
    }
    
    
    /** My name as an integer.*/
    private int _n;
    /** My first node's name.*/
    private int _a;
    /** My second node's name.*/
    private int _b;
    /** My third node's name.*/
    private int _c;
    /** My fourth node's name.*/
    private int _d;
    /** My Volume in m^3 */
    private double _volume;
    
    /** My Interpolating function Coefficient Matrix.
     * a + bx +cy + dz*/
    private RealMatrix coeffMat;
    
    /** My array of nodes.*/
    private Node[] nodes;
    
    /** My PrintWriter.*/
    private PrintWriter writer;
    
    /** My verbose toggle.*/
    private boolean _verbose;
    
}
