/**
 * \package simulator.detachment
 * 
 * \brief Package of classes that capture detachment of agents from the biomass
 * 
 * Package of classes that capture detachment of agents from the biomass. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.detachment;

import simulator.AgentContainer;
import simulator.Simulator;
import simulator.agent.LocatedGroup;
import utils.XMLParser;

/**
 * \Quadratic detachment method for granules based on DS_Quadratic. Function: kDet*r^2 or kDet*r^3, where r is distance from the domain center.
 * 
 * The mark-up defines the erosion forces that act on the biofilm surface. Detachment works by removing a layer of biomass based on the 
 * detachment speed and the timestep, with the detachment speed calculated via one of the given forms. This class captures the Quadratic 
 * detachment method for granules.
 * 
 * @author Bastiaan Cockx
 *
 */
public class DS_Granule extends LevelSet {

	/**
	 * Constant parameter used to determine the strength of detachment.
	 */
	private double kDet;
	
	/**
	 * Maximum thickness that the biofilm may reach
	 */
	private double maxTh;

	/**
	 * Domain X dimension in micrometers.
	 */
	private double length_X;
	
	/**
	 * Domain Y dimension in micrometers.
	 */
	private double length_Y;
	
	/**
	 * Domain Z dimension in micrometers.
	 */
	private double length_Z;
	
	/**
	 * Whether this computation domain is two or three dimensional.
	 */
	private boolean is3D;

	/**
	 * squared distance from the domain center (r^2), r^3 for 3D granules.
	 */
	private double r;	

	/**
	 * \brief Initialise the object by reading attributes from associated agent grid and XML protocol file
	 * 
	 * Initialise the object by reading attributes from associated agent grid and XML protocol file
	 * 
	 * @param anAgentGrid	Associated grid of agents
	 * @param root	XML tag containing information related to this detachment mechanism
	 */
	@Override
	public void init(AgentContainer anAgentGrid, XMLParser root){
		super.init(anAgentGrid, root);
		// kDet has units of: um-1.hr-1
		// this gives speed in um.hr-1
		kDet = root.getParamDbl("kDet");
		double value=root.getParamDbl("maxTh");
		maxTh=(Double.isNaN(value)? Double.POSITIVE_INFINITY:value);

		// world info
		length_X = anAgentGrid.domain.length_X;
		length_Y = anAgentGrid.domain.length_Y;
		length_Z = anAgentGrid.domain.length_Z;
		is3D = anAgentGrid.domain.is3D;
	}

	/**
	 *\brief Calculate and return the local detachment speed using this detachment method
	 *
	 * Calculate and return the local detachment speed using this detachment method
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aGroup	Located group for which the local detachment speed is being determined
	 * @return Double stating local detachment speed for this group
	 *
	 */
	@Override
	protected double getLocalDetachmentSpeed(LocatedGroup aGroup, Simulator aSim) {
		
		if (is3D) {
			// calculate distance^2 from center of the domain
			r = (aGroup.cc.x-(length_X*0.5))*(aGroup.cc.x-(length_X*0.5))
					+(aGroup.cc.y-(length_Y*0.5))*(aGroup.cc.y-(length_Y*0.5))
					+(aGroup.cc.z-(length_Z*0.5))*(aGroup.cc.z-(length_Z*0.5));
			
			// if max granule diameter is exceeded return max value.
			if (r>(0.5*maxTh)*(0.5*maxTh)*(0.5*maxTh)) return Double.MAX_VALUE;
		} else {
			r = (aGroup.cc.x-(length_X*0.5))*(aGroup.cc.x-(length_X*0.5))
					+(aGroup.cc.y-(length_Y*0.5))*(aGroup.cc.y-(length_Y*0.5));
			
			// if max granule diameter is exceeded return max value.
			if (r>(0.5*maxTh)*(0.5*maxTh)) return Double.MAX_VALUE;
		}
		
		//return detachment rate
		return kDet*r;
	}	

}
