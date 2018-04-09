package cavityFETD;

import org.apache.commons.math3.util.FastMath;

/** A Node is a point in space with x, y, and z coordinates
 *  @author Darwin Li
 */
public class Node {

	/** The constructor that gives this node its
	 *name N as integer and X and Y Z values as double. 
	 */
    public Node(int n, double x, double y, double z, boolean verbose) {
    	_n = n;
        _x = x;
        _y = y;
        _z = z;
        _bound = "none";
        _monitor = 0;
        _verbose = verbose;
        if (_verbose) {
        	System.out.println("Node" + _n + ": " + _x + "," + _y + "," + _z + " bound type " + _bound);
        }
    }
    
    /** Return my name as an integer.*/
    int getNumber() {
        return _n;
    }
    
    /** Return my X value in meters.*/
    double getX() {
        return _x;
    }
    
    /** Return my Y value in meters.*/
    double getY() {
        return _y;
    }
    
    /** Return my Y value in meters.*/
    double getZ() {
        return _z;
    }
    
    /** Return my the node I belong to. A value
     * of 0 means that I am an inner node*/
    String getBound() {
        return _bound;
    }
    
    /** Return my the node I belong to. A value
     * of 0 means that I am an inner node*/
    int getMonitor() {
        return _monitor;
    }
    
	 /** Takes in PEC as a set of planes that are PEC boundaries. Takes the
	  * x-y position of the current source as an array and changes the
	  * bound identification accordingly. Also takes in the x-y position
	  * of a voltage monitor and changes the boundary identification
	  * accordingly. */
    void setBound(double pecZMin, double pecZMax, double[] jPosition, double [] vMonPos) {
        if (pecZMin == _z) {
        	_bound = "zMin";
        } else if (pecZMax == _z) {
        	_bound = "zMax";
        }
        if (FastMath.abs(jPosition[0] - _x)/jPosition[0] <=.00001 && FastMath.abs(jPosition[1] -_y)/jPosition[1] <=.00001) {
        	_monitor = 3;
        } else if (FastMath.abs(vMonPos[0] - _x)/vMonPos[0] <=.00001 && FastMath.abs(vMonPos[1] -_y)/vMonPos[1] <=.00001) {
        	_monitor = 4;
        }
    }
    
	 /** Takes in PEC as a set of planes that are PEC boundaries. Then, applies the boundaries accordingly and
	  * returns 1 if I am on a PEC boundary and 0 if I am a free Node. */
    int setBound(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax) {  	
    	if (zMin == _z) {
    		_bound = _bound + "zMin";
    	} else if (zMax == _z) {
    	  	_bound = _bound + "zMax";
    	  }
    	if (yMin == _y) {
    	  	_bound = _bound + "yMin";
    	} else if (yMax == _y) {
    	    	_bound = _bound + "yMax";
    	}
    	if (xMin == _x) {
    	  	_bound = _bound + "xMin";
    	} else if (xMax == _x) {
    	 	_bound = _bound + "xMax";
    	}
    	if (_bound.equals("none")) {
    		return 0;
    	} else {
    		return 1;
    	}
    }
    
    /** My name.*/
    private int _n;
    
    /** My x value.*/
    private double _x;
    
    /** My y value.*/
    private double _y;
    
    /** My z value.*/
    private double _z;
    
    /** The boundary I belong to if I belong to one or more Dirichilet Boundaries*/
    private String _bound;
    
    /** A value of 3 means that I am located
     * at a source. A value of 4 means that I am located at the
     * user defined voltage monitor. */
    private int _monitor;

    /** My verbose toggle*/
    private boolean _verbose;
    
}
