package cavityFETD;

/** Relevant Physics Constants.
 * @author Darwin Li
 */
class PhysicsConstants {
	
	/** The Permittivity of Free Space */
	public static Double e_o = 8.85418782e-12;
	
	/** The Permeability of Free Space */
	public static Double u_o = 1.25663706e-6;
	
	/** The speed of light in meters per second */
	public static Double c_o = 299792458.00;
	
	/** The relative dieletric constant */
	public static double e_r = 1.0;
	
	/** The effective dieletric constant */
	public static double e_f = e_r*e_o;
	
	/** The relative permeability constant */
	public static double u_r = 1.0;
	
	/** The effective permeability constant */
	public static double u_f = u_r * u_o;
	
	/** The conductivity of the dieletric in S/m */
	public static double sigma = 0.005;
}
