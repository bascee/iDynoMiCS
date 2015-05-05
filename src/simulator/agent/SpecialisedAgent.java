///**
// * \package agent
// * \brief Package of utilities that create and manage agents in the simulation and their participation in relevant reactions
// * 
// * Package of utilities that create and manage agents in the simulation and their participation in relevant reactions. This package is 
// * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
// * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
// * the following URL  "http://www.cecill.info".
// */
//package simulator.agent;
//
//import simulator.SpatialGrid;
//
///**
// * \brief Extension of Agent class, adds location and parameter information
// * for an object of a particular species in the simulation.
// * 
// * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
// * for Infection Research (Germany)
// * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
// */
//public abstract class SpecialisedAgent extends Agent implements Cloneable 
//{
//	/**
//	 * Boolean noting whether this agent is still active in the simulation
//	 */
////	public Boolean           isDead = false;
//	
//	/**
//	 * Reason for agent's death. Added by Sonia Martins April 2010
//	 */
////	public String death;
//
//	/**
//	 * \brief Creates a SpecialisedAgent object and initialises the object in which associated parameters are stored
//	 * 
//	 * Creates a SpecialisedAgent object and initialises the object in which associated parameters are stored
//	 */
//	public SpecialisedAgent() 
//	{
//		// Call constructor of parent class
//		super();
//		_speciesParam = new SpeciesParam();
//	}
//
//	/**
//	 * \brief Clones this agent object, creating a new progeny of this agent. Ensures new clone inherits same parameters as parents
//	 * 
//	 * Clones this agent object, creating a new progeny of this agent. Ensures new clone inherits same parameters as parents
//	 * 
//	 * @throws CloneNotSupportedException	Exception should the class not implement Cloneable
//	 */
////	@Override
////	public Object clone() throws CloneNotSupportedException
////	{
////		Agent out = (Agent) super.clone();
////
////		// Copy the references (superficial copy)
////		out._species = this._species;
////		out._speciesParam = this._speciesParam;
////		return (Object) out;
////	}
//
//	/**
//	 * \brief Create a new agent with mutated parameters based on species
//	 * default values.
//	 */
//
////	public abstract void createNewAgent();
//
//
//	
//	/**
//	 * \brief Used in the calculation of delta move in Agent Container class.
//	 * 
//	 * Will investigate further why this class returns 0.
//	 * 
//	 * @return	0
//	 */
////	public Double move()
////	{
////		return 0.0;
////	}
//	
//	/**
//	 * \brief Models a mechanical interaction between two located agents.
//	 * 
//	 * Implemented by extending classes (LocatedAgent)
//	 * 
//	 * @param MUTUAL	Whether movement is shared between two agents or 
//	 * applied only to this one.
//	 * @return	The move to be applied once the shoving or pull calculations
//	 * have been performed.
//	 */
////	public Double interact(boolean MUTUAL)
////	{
////		return 0.0;
////	}
//
//	
////	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex) { 
////	}
////	
////	public void fitMassOnGrid(SpatialGrid aSpG) { 
////	}
////	
////	public void fitVolRateOnGrid(SpatialGrid aSpG) {
////	}
//
//}
