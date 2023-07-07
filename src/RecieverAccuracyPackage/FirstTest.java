/*
 * This class (FirstTest) is the demonestration of orekit ISL orbit determination capability 
 * using OneWayGNSSRange measurements
 *  Further developements are under consideration.
 * 
 * 
 * 
 */
package RecieverAccuracyPackage;

import java.io.*;
import org.orekit.time.*;
import org.orekit.data.*;
import org.orekit.orbits.*;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.utils.*;
import org.orekit.frames.*;
import org.orekit.models.earth.ReferenceEllipsoid;
import org.orekit.attitudes.*;
import org.orekit.forces.gravity.potential.*;
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.propagation.conversion.*;
import org.hipparchus.linear.*;
import org.hipparchus.optim.nonlinear.vector.leastsquares.GaussNewtonOptimizer;
import org.orekit.estimation.leastsquares.*; 
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.estimation.measurements.gnss.OneWayGNSSRange;
import org.orekit.propagation.*;

public class FirstTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
    	File orekitData = new File("C:\\Program Files\\Java\\jdk-19\\lib"); // Path to orekit-data.zip
    	manager.addProvider(new DirectoryCrawler(orekitData));
    	TimeScale    utc   = TimeScalesFactory.getUTC();
    	AbsoluteDate startTime = new AbsoluteDate(2018, 4, 26, 8, 30, 0.0, utc);
    	final double RAD2DEG = 180 / Math.PI;
    	final double DEG2RAD = 1 / RAD2DEG;
    	KeplerianOrbit initialOrbitUn = new KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 700000.0,
    			0.0002654 ,    
                65.0*DEG2RAD, 
                0.0*DEG2RAD,  
                11.0*DEG2RAD ,
                30.0*DEG2RAD,
                PositionAngle.TRUE,
                FramesFactory.getEME2000(),
                startTime,
                Constants.WGS84_EARTH_MU);
    	KeplerianOrbit initialOrbitKn = new KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 700000.0,
    			0.0002654 ,    
                65.0*DEG2RAD, 
                0.0*DEG2RAD,  
                10.0*DEG2RAD ,
                30.0*DEG2RAD,
                PositionAngle.TRUE,
                FramesFactory.getEME2000(),
                startTime,
                Constants.WGS84_EARTH_MU);
    	
    	Frame eci = FramesFactory.getGCRF();
    	Frame ecef = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, false);
    	ReferenceEllipsoid wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(ecef);
    	CartesianOrbit initialEciUn = new CartesianOrbit(initialOrbitUn.getPVCoordinates(eci), FramesFactory.getGCRF(), wgs84Ellipsoid.getGM());
    	
    	NadirPointing nadirPointing = new NadirPointing(eci, wgs84Ellipsoid);
    	NormalizedSphericalHarmonicsProvider gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(1, 1);
    	HolmesFeatherstoneAttractionModel gravityAttractionModel = new HolmesFeatherstoneAttractionModel(ecef, gravityProvider);
    	DormandPrince853IntegratorBuilder integratorBuilder = new DormandPrince853IntegratorBuilder(0.01, 2.0, 3.0);
    	NumericalPropagatorBuilder odPropBuilderUn = new NumericalPropagatorBuilder(initialEciUn,
                integratorBuilder, PositionAngle.MEAN, 1.0);
    	
    	
    	odPropBuilderUn.setMass(100.0);
    	odPropBuilderUn.setAttitudeProvider(nadirPointing);
    	odPropBuilderUn.addForceModel(gravityAttractionModel);
    	
    	
    	
    	QRDecomposer matrixDecomposer = new QRDecomposer(1e-6);
    	GaussNewtonOptimizer optimizer = new GaussNewtonOptimizer(matrixDecomposer, false);
    	
    	BatchLSEstimator estimator1 = new BatchLSEstimator(optimizer, odPropBuilderUn);
    	ObservableSatellite TBOSat = new ObservableSatellite(0);
    	double estimator_convergence_thres = 0.01;
    	int estimator_max_iterations = 25;
    	int estimator_max_evaluations = 35;
    	
    	estimator1.setParametersConvergenceThreshold(estimator_convergence_thres);
    	estimator1.setMaxIterations(estimator_max_iterations);
    	estimator1.setMaxEvaluations(estimator_max_evaluations);
    	
    	KeplerianPropagator propagatorUn = new KeplerianPropagator(initialOrbitUn);
    	KeplerianPropagator propagatorKn = new KeplerianPropagator(initialOrbitKn);
    	
    	AbsoluteDate propDate = startTime;
    	AbsoluteDate endTime = startTime.shiftedBy(360.0);
    	
    	double SatsClockOffset = 0.0;//Math.pow(10,-9);
    	double sampleRate = 0.2; // Data Sampling frequency
    	double samplesTimeShift = 1.0 / sampleRate;
    	
    	double theoreticalStD = 1.0;
    	double baseWeight = 1.0;
    	
    	while (propDate.compareTo(endTime) <= 0) {
    		SpacecraftState pvInerUn = propagatorUn.propagate(propDate);
    		SpacecraftState pvInerKn = propagatorKn.propagate(propDate);
    		Vector3D pvInerKnVector = pvInerKn.getPVCoordinates().getPosition();
    	    double realRange = pvInerUn.getPVCoordinates().getPosition().distance(pvInerKnVector);
    	    double perturbationNoise = Math.random()*Math.pow(10,-2); // Noise 0.1-10 [mm]
    	    double perturbedRange = realRange + perturbationNoise;
    	    OneWayGNSSRange onewayMeas = new OneWayGNSSRange(propagatorKn, SatsClockOffset,
    	    								propDate,
    	    								perturbedRange,
    	    								theoreticalStD,
    	    								baseWeight,
    	                                       TBOSat);
    	    
    	    propDate = propDate.shiftedBy(samplesTimeShift);
    	    estimator1.addMeasurement(onewayMeas);
    	}
    	
    	Propagator estimatedPropagator = estimator1.estimate()[0];
    	SpacecraftState	estimatedInitialState = estimatedPropagator.getInitialState();
    	SpacecraftState	RealPosition = propagatorUn.getInitialState();
    	
    	
    	
    	System.out.println(estimatedInitialState.getPVCoordinates().getPosition());
    	System.out.println(RealPosition.getPVCoordinates().getPosition());

	}

}
