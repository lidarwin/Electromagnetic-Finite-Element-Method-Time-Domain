package cavityFETD;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.util.FastMath;


/** This class is used to compute the
 * Gaussian Excitation signal. Currently,
 * arbitrary excitation signals are not supported,
 * but can easily be added by parsing ASCII files
 *  @author Darwin Li
 */
public class ExcitationSignal {

	/** The constructor that takes in FMIN as the minimum
	 * frequency in hertz and FMAX as the maximum frequency in hertz.
	 * TMAX is in units of seconds and makes sure that the simulation
	 * truncates. */
	public ExcitationSignal(double fMinn, double fMaxx, String outDir) {
		fMin = fMinn;
		fMax = fMaxx;
    	try {
			writer = new PrintWriter(outDir + "\\" + "ExcitationSignalLog.txt", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}    
		bbw=0.0001;
		bt=0.001;
		bf=0.001;
		fc = (fMax+fMin)/2;
		t0=-FastMath.log10(bt);
		t0 = FastMath.sqrt(t0);
		alpha = FastMath.PI*(fMax-fMin)/FastMath.sqrt(-FastMath.log(bbw));
		t0=FastMath.sqrt(-FastMath.log10(bt))/alpha;

		if (fc < alpha*FastMath.PI/FastMath.sqrt(-FastMath.log(bf))) {
			writer.println("DC component is too large. Please consider using multiple frequency bands.");
		}
		dt = 1/(2*fMax);
	}
	
	/** This returns the excitation signal in a 2 column array with
	 * the first column as time in seconds and the second as the value */
	double[][] computeSignal(int maxStep) {
		double[][] tempSignal = new double[maxStep][2];
		double time;
		for (int t = 0; t < maxStep; t++) {
			time = t * dt;
			tempSignal[t][0]= time;
			double expTerm = t*dt - t0;
			expTerm = expTerm*alpha;
			expTerm = -FastMath.pow(expTerm, 2.0);
			double cosTerm = t*dt-t0;
			cosTerm = cosTerm * 2.0 * FastMath.PI*fc;
			cosTerm = FastMath.cos(cosTerm);
			tempSignal[t][1] = FastMath.exp(-FastMath.pow(((t)*dt-t0)*alpha, 2))*FastMath.cos(2*FastMath.PI*fc*((t)*dt-t0));
			writer.println(tempSignal[t][0]+ "	" + tempSignal[t][1]);
		}
		writer.close();
		return tempSignal;
	}
	
	/** This returns the analytical first order derivative of the excitation signal in a 2 column array with
	 * the first column as time in seconds and the second as the value over seconds */
	double[][] computeSignalDerivativeAnalytical(int maxStep) {
		double[][] tempSignal = new double[maxStep][2];
		double time;
		for (int t = 0; t < maxStep; t++) {
			time = t * dt;
			tempSignal[t][0]= time;
			tempSignal[t][1] = 2*FastMath.pow(alpha, 2)*dt*FastMath.exp(-FastMath.pow(alpha,2)*
					FastMath.pow(t0 - dt*t, 2))*FastMath.cos(2*fc*FastMath.PI*(t0 - dt*t))*(t0 - dt*t) + 
					2*dt*fc*FastMath.PI*FastMath.exp(-FastMath.pow(alpha,2)*FastMath.pow(t0 - dt*t,2))*
							FastMath.sin(2*fc*FastMath.PI*(t0 - dt*t));
		}
		return tempSignal;
	}
	
	/** This returns the first order differentiation by data of the excitation signal in a 2 column array with
	 * the first column as time in seconds and the second as the value over seconds */
	double[][] computeSignalDerivativeHelper(double[][] signal) {
		int maxStep = signal.length - 1;
		double[][] tempSignal = new double[maxStep][2];
		double time;
		for (int t = 0; t < maxStep; t++) {
			time = t * dt;
			tempSignal[t][0]= time;
			tempSignal[t][1]= (signal[t+1][1] - signal[t][1])/dt;
		}
		return tempSignal;
	}
	
	/** This returns the first order differentiation by data of the excitation signal in a 2 column array with
	 * the first column as time in seconds and the second as the value over seconds */
	double[][] computeSignalDerivative(int maxStep) {
		return computeSignalDerivativeHelper(computeSignal(maxStep + 1));
	}

	/** Returns dt, time step in seconds */
	double getDT() {
		return dt;
	}
	
	/** This constant is the time step in seconds */
	private double dt;
	
	/** My PrintWriter */
	private PrintWriter writer;
	
	/** My minimum, maximum, and center frequencies. In units of Hertz */
	private double fMin;
	private double fMax;
	private double fc;
	
	/** These are constants used to compute the gaussian source. It is a cosine function with a ramp */
	private double bbw;
	private double bt;
	private double bf;
	private double alpha;
	private double t0;
	
}
