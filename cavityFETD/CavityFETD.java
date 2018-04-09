package cavityFETD;

import static cavityFETD.FormatException.getTotalErrors;
import static cavityFETD.FormatException.reportError;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** Compute the PUL Capacitance of a simple coax given a .poly file.
 * @author Darwin Li
 */
public class CavityFETD {
	
    /** ARGS must have length 3. The ARGS[0] gives the verbose command, ARGS[1] gives the
     *  .ele file, ARGS[2] gives the .node file Print a usage message 
     *  otherwise or if the files are unreadable or unwritable, respectively. */
	public static void main(String[] args) {
		double fMin = 0.0;
		double fMax = 0.0;
		double jXPos = 0.0;
		double jYPos = 0.0;
		double vMonXPos = 0.0;
		double vMonYPos = 0.0;
	    int timeStep = 0;
		Scanner scanner = new Scanner(System.in);
		System.out.println("Enter gmsh .mesh File Path name ie. C:\\Work\\mesh.msh");
		String meshfilepathIn = scanner.next();
		System.out.println("Specify results output directory ie. C:\\Work\\");
		Pattern resultsDirPatt = Pattern.compile("(.)+");
		String resultsDirIn = scanner.next();
		if (!resultsDirIn.endsWith("\\")) {
			resultsDirIn = resultsDirIn + "\\";
		}
		Pattern freqRangePatt = Pattern.compile("([0-9.eE]+)");
		boolean freqRangeMatBool = false;
		Matcher freqRangeMat;
		String freqRangeIn;
		while (!freqRangeMatBool) {
			System.out.println("Please Enter fMin in Hz");
			freqRangeIn = scanner.next();
			freqRangeMat = freqRangePatt.matcher(freqRangeIn);
			freqRangeMatBool = freqRangeMat.find();
			if (freqRangeMatBool) {
				fMin = Double.valueOf(freqRangeMat.group(1));
			}
		}
		freqRangeMatBool = false;
		while (!freqRangeMatBool) {
			System.out.println("Please Enter fMax in Hz");
			freqRangeIn = scanner.next();
			freqRangeMat = freqRangePatt.matcher(freqRangeIn);
			freqRangeMatBool = freqRangeMat.find();
			if (freqRangeMatBool) {
				fMax = Double.valueOf(freqRangeMat.group(1));
			}
		}
		if (fMin > fMax) {
			System.out.println("Oops, fMin is greater than fMax. Please re-run");
			System.exit(0);
		}
		System.out.println("Please choose from one of the following simulation types by entering 1 or 2");
		System.out.println("1) PEC Cavity Eigenmode Analysis with FEM");
		System.out.println("2) PDN Dual Plane Impedance with Finite Element Time Domain (beta)");
		Pattern problemTypePatt = Pattern.compile("1|2");
		String problemTypeIn = scanner.next();
		Matcher problemTypeMat = problemTypePatt.matcher(problemTypeIn);
		while (!problemTypeMat.find()) {
			System.out.println("Incorrect input");
			problemTypeIn = scanner.next();
			problemTypeMat = problemTypePatt.matcher(problemTypeIn);
		}
		int problemType = Integer.valueOf(problemTypeMat.group(0));
		if (problemType == 2) {
			System.out.println("Enter where you would like the current source to be in meters for x and y ie. 50.0,50.0");
			System.out.println("If the location is not valid in the problem domain, then no results will return!");
			System.out.println("In that case, look in the logfile for a clue to a valid position");
			Pattern jPosPatt = Pattern.compile("([0-9\\.eE]+),([0-9\\.eE]+)");
			String jPosIn = scanner.next();
			Matcher jPosMat = jPosPatt.matcher(jPosIn);
		    boolean jPosBool = jPosMat.find();
			while (!jPosBool) {
				System.out.println("Incorrect input!");
				jPosIn = scanner.next();
				jPosMat = jPosPatt.matcher(jPosIn);
				jPosBool = jPosMat.find();
			}
			jXPos = Double.valueOf(jPosMat.group(1));
			jYPos = Double.valueOf(jPosMat.group(2));
			System.out.println("Enter where you would like to monitor Voltage in meters for x and y ie. 50.0,50.0");
			System.out.println("If the location is not valid in the problem domain, then no results will return!");
			System.out.println("In that case, look in the logfile for a clue to a valid position");
			Pattern vMonPosPatt = Pattern.compile("([0-9\\.eE]+),([0-9\\.eE]+)");
			String vMonPosIn = scanner.next();
			Matcher vMonPosMat = vMonPosPatt.matcher(vMonPosIn);
		    boolean vMonPosBool = vMonPosMat.find();

			while (!vMonPosBool) {
				System.out.println("Incorrect input!");
				vMonPosIn = scanner.next();
				vMonPosMat = vMonPosPatt.matcher(vMonPosIn);
				vMonPosBool = vMonPosMat.find();
			}
			vMonXPos = Double.valueOf(vMonPosMat.group(1));
			vMonYPos = Double.valueOf(vMonPosMat.group(2));
			System.out.println("Enter maximum number of time steps in integer form ie. 1000");
			Pattern timeStepPatt = Pattern.compile("([0-9\\.eE]+)");
			String timeStepIn = scanner.next();
			Matcher timeStepMat = timeStepPatt.matcher(timeStepIn);
		    boolean timeStepBool = timeStepMat.find();
			while (!timeStepBool) {
				System.out.println("Incorrect input!");
				timeStepIn = scanner.next();
				timeStepMat = timeStepPatt.matcher(timeStepIn);
				timeStepBool = timeStepMat.find();
			}
			timeStep = Integer.valueOf(timeStepMat.group(0));
		}
		System.out.println("Would you like to keep a log? Enter anything except \"no\" to keep a log and have outputs in console");
		System.out.println("Log Files are kept in the java directory");
		String verboseIn = scanner.next();
		boolean verbose = true;
		if (verboseIn.equals("no")) {
			verbose = false;
		}
		System.out.println("Simulation Started");
		try {
			writerImport = new PrintWriter(resultsDirIn+"LogImport.txt", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		try {
			solverLogfile = new PrintWriter(resultsDirIn+"VoltageResults.txt", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
        MatrixCalc matCalc;
        Solver solve;
        if (problemType == 2) {
            matCalc = new MatrixCalcFETD(verbose, meshfilepathIn, jXPos, jYPos, vMonXPos, vMonYPos, writerImport);
            ExcitationSignal excit;
            excit = new ExcitationSignal(fMin, fMax, resultsDirIn);
            solve = new SolverFETD(timeStep, excit, matCalc, resultsDirIn);
            solve.solverLoop();
        } else {
        	matCalc = new MatrixCalcEigenmode(verbose, meshfilepathIn, writerImport);
            solve = new SolverEigenmode(fMin, fMax, matCalc, resultsDirIn);
            solve.solverLoop();
        }
        writerImport.close();
        solverLogfile.close();
        System.exit(getTotalErrors() == 0 ? 0 : 1);
    } 
	
    /** Print usage message. */
    private static void usage() {
        System.out.printf("Usage: java format.Main INFILE [OUTFILE]%n"
                          + "   Format INFILE, sending output to OUTFILE "
                          + "(default: standard output).%n");
    }
	
    /** PrintWriter for the main log file.*/
    private static PrintWriter writer;
    
    /** PrintWriter for the mesh file import log file.*/
    private static PrintWriter writerImport;
    
    /** PrintWriter for the Solver log file.*/
    private static PrintWriter solverLogfile;
	
}
