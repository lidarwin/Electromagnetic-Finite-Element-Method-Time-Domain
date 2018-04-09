package cavityFETD;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.HashSet;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.io.File;


/** Run the Simulation Solver loop for a Time Domain Finite Element treatment of dual ground plane
 * @author Darwin Li
 */
class SolverFETD extends Solver {

	SolverFETD(int tStepMax, ExcitationSignal signal, MatrixCalc matCalc, String outDir) {
		beta = 0.25;
		_outDir = outDir;
		_maxSteps = tStepMax;
		_signal = signal;
		_plateDistance = matCalc.getPlateDistance();
		_dt = _signal.getDT();
		_dtSquare = FastMath.pow(_dt, 2.0);
		if (matCalc.getVerbose()) {
			_verbose = true;
		}
		RealMatrix S = matCalc.getS().scalarMultiply(FastMath.pow(PhysicsConstants.u_f, -1));
		RealMatrix T = matCalc.getT().scalarMultiply(PhysicsConstants.e_f);
		RealMatrix R = matCalc.getR().scalarMultiply(PhysicsConstants.sigma);
		A = T.scalarMultiply(1/(_dtSquare));
		A = A.add(R.scalarMultiply(1/(2*_dt)));
		A = A.add(S.scalarMultiply(beta));
		B = T.scalarMultiply(2/(_dtSquare));
		B = B.subtract(S.scalarMultiply(1-2*beta));
		C = T.scalarMultiply(1/(_dtSquare));
		C = C.subtract(R.scalarMultiply(1/(2*_dt)));
		C = C.add(S.scalarMultiply(beta));
		C = C.scalarMultiply(-1.0);
		//_Asolver = new LUDecomposition(A).getSolver();
		A = MatrixUtils.inverse(A);
		f = matCalc.getF();
		_voltageSourceMonitor = matCalc.getSourceVoltageMonitor().toArray(new Integer[0]);
		_voltageMonitor1 = matCalc.getVoltageMonitor1().toArray(new Integer[0]);
		_backwardsEdgeSet = matCalc.getBackwardsSet();
	}
	
	/** Run the solver loop for FETD simulation. Outputs are saved in the directory specified in Main*/
	@Override
	void solverLoop() {
		XYSeries currentSig = new XYSeries("Current A vs seconds");
		XYSeries currentDerivSig = new XYSeries("Current Derivative A vs second^2");
		XYSeries voltageSourceSig = new XYSeries("Voltage at the source V vs second");
		XYSeries voltageMonitorSig = new XYSeries("Voltage at the monitor Vvs second");
		PrintWriter currWriter = null;
    	try {
			currWriter = new PrintWriter(_outDir + "\\" + "CurrentSignal.txt", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		PrintWriter currDerWriter = null;
    	try {
			currDerWriter = new PrintWriter(_outDir + "\\" + "CurrentDerivativeSignal.txt", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		PrintWriter vS_out = null;
		try {
			vS_out = new PrintWriter(_outDir + "\\" + "VsourceSignal.txt", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		PrintWriter v1_out = null;
		try {
			v1_out = new PrintWriter(_outDir + "\\" + "V1Signal.txt", "UTF-8");
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		} catch (UnsupportedEncodingException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		double[][] tempI = _signal.computeSignal(_maxSteps);
		Array2DRowRealMatrix currentSignal = new Array2DRowRealMatrix(tempI);
		Array2DRowRealMatrix currentDerivative = new Array2DRowRealMatrix(_signal.computeSignalDerivative(_maxSteps+1));
		Array2DRowRealMatrix voltageSource = new Array2DRowRealMatrix(_maxSteps, 2);
		voltageSource.setColumnMatrix(0, currentSignal.getColumnMatrix(0));
		Array2DRowRealMatrix voltage1 = (Array2DRowRealMatrix) voltageSource.copy();
		//Set initial conditions to be 0
		_fnMinus1 = f.mapMultiply(0.0);
		_En = f.mapMultiply(0.0);
		_EnMinus1 = f.mapMultiply(0.0);
		//Start the solver loop
		for (int i = 0; i < _maxSteps; i ++) {
			_fn = f.mapMultiply(currentDerivative.getEntry(i, 1));
			_fnPlus1 = f.mapMultiply(currentDerivative.getEntry(i+1, 1));
			
			_EnPlus1 = _fnPlus1.mapMultiply(beta);
			_EnPlus1 = _EnPlus1.add(_fn.mapMultiply(1.0-2.0*beta));
			_EnPlus1 = _EnPlus1.add(_fnMinus1.mapMultiply(beta));
			_EnPlus1 = B.operate(_En).subtract(_EnPlus1);
			_EnPlus1 = _EnPlus1.add(C.operate(_EnMinus1));
			//_EnPlus1 = _Asolver.solve(_EnPlus1);
			_EnPlus1 = A.operate(_EnPlus1);
			RealVector tempB = B.operate(_En);
			RealVector tempC = C.operate(_EnMinus1);
			RealVector tempFnPlus1 = _fnPlus1.mapMultiply(beta);
			RealVector tempFnMinus1 = _fnMinus1.mapMultiply(beta);
			RealVector tempFn = _fn.mapMultiply(1-(2*beta));
			_EnPlus1 = tempB;
			_EnPlus1 = _EnPlus1.add(tempC);
			_EnPlus1 = _EnPlus1.subtract(tempFnPlus1);
			_EnPlus1 = _EnPlus1.subtract(tempFn);
			_EnPlus1 = _EnPlus1.subtract(tempFnMinus1);
			_EnPlus1 = A.operate(_EnPlus1);
			
			for (Integer e: _voltageSourceMonitor) {
				double integ = _En.getEntry(e);
				integ = integ * _plateDistance/(_voltageSourceMonitor.length);
				if (_backwardsEdgeSet.contains(e)) {
					integ = -1.0*integ;
				}
				voltageSource.addToEntry(i, 1, integ);
			}
			for (Integer e: _voltageMonitor1) {
				double integ = _En.getEntry(e);
				integ = integ * _plateDistance/(_voltageMonitor1.length);
				if (_backwardsEdgeSet.contains(e)) {
					integ = -1.0*integ;
				}
				voltage1.addToEntry(i, 1, integ);
			}
			currWriter.println(currentSignal.getEntry(i, 0) + "    " + currentSignal.getEntry(i, 1));
			currDerWriter.println(currentDerivative.getEntry(i, 0) + "    " + currentDerivative.getEntry(i, 1));
			vS_out.println(voltageSource.getEntry(i, 0)+ "    " + voltageSource.getEntry(i, 1));
			v1_out.println(voltage1.getEntry(i, 0)+ "    " + voltage1.getEntry(i, 1));
			currentSig.add(currentSignal.getEntry(i, 0), currentSignal.getEntry(i, 1));
			currentDerivSig.add(currentDerivative.getEntry(i, 0), currentDerivative.getEntry(i, 1));
			voltageSourceSig.add(voltageSource.getEntry(i, 0), voltageSource.getEntry(i, 1));
			voltageMonitorSig.add(voltage1.getEntry(i, 0), voltage1.getEntry(i, 1));
			_fnMinus1 = _fn;
			_EnMinus1 = _En;
			_En = _EnPlus1;
		}
		currWriter.close();
		currDerWriter.close();
		vS_out.close();
		v1_out.close();
		XYSeriesCollection currentSigData = new XYSeriesCollection();
		XYSeriesCollection currentDerivData = new XYSeriesCollection();
		XYSeriesCollection voltageData = new XYSeriesCollection();
		currentSigData.addSeries(currentSig);
		currentDerivData.addSeries(currentDerivSig);
		voltageData.addSeries(voltageSourceSig);
		voltageData.addSeries(voltageMonitorSig);
		JFreeChart voltageChart = ChartFactory.createXYLineChart("Voltage vs time", "Seconds", "Volts", voltageData);
		JFreeChart currentChart = ChartFactory.createXYLineChart("Current vs time", "Seconds", "Amps", currentSigData);
		JFreeChart currentDerivChart = ChartFactory.createXYLineChart("Current Derivative vs time", "Seconds", "Amps/s", currentDerivData);
		try {
			ChartUtilities.saveChartAsJPEG(new File(_outDir + "\\" + "Voltages.jpg"), voltageChart, 1000, 600);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}
		try {
			ChartUtilities.saveChartAsJPEG(new File(_outDir + "\\" + "Current.jpg"), currentChart, 1000, 600);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}
		try {
			ChartUtilities.saveChartAsJPEG(new File(_outDir + "\\" + "CurrentDerivative.jpg"), currentDerivChart, 1000, 600);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}
		}
	
	/** My Decomposition Solver for the A matrix*/
	private DecompositionSolver _Asolver;
	
	/** The list of edges that are backwards in global numbering. */
	private HashSet<Integer> _backwardsEdgeSet;
}