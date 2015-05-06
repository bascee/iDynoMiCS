///**
// * \package agent
// * \brief Package of utilities that create and manage agents in the simulation
// * and their participation in relevant reactions.
// * 
// * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
// * under French law and abides by the rules of distribution of free software.  
// * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
// * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL
// * "http://www.cecill.info".
// */
//package simulator.agent;
//
//import java.util.Iterator;
//
//import simulator.*;
//
///**
// * \brief Extends ActiveAgent by adding functionality to control agent grid
// * location, agent shoving, agent death and division, and agent movement.
// *  
// * During each global timestep, agent divisions and agent growth lead to many
// * cases where neighbouring agents will overlap. A relaxation algorithm is
// * used to find iteratively the new overlap-minimising steady state
// * configuration of agent locations at the end of each timestep. 
// * 
// * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
// * for Infection Research (Germany)
// * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
// * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology,
// * University of Birmingham (UK)
// * @author Rob Clegg (rjc096@bham.ac.uk), Centre for Systems Biology,
// * University of Birmingham (UK)
// */
//public abstract class LocatedAgent extends ActiveAgent implements Cloneable 
//{
//	/**
//	 * \brief Constructor used to generate progenitor and initialise an object
//	 * to store relevant parameters. 
//	 */
//	public LocatedAgent()
//	{
//		super();
//		_speciesParam = new LocatedParam();
//	}
//
//	/* ______________________ SHOVING ___________________________________ */
//
//	/* _______________ FILE OUTPUT _____________________ */
//
//	/* _______________ RADIUS, MASS AND VOLUME _____________________ */
//
//}