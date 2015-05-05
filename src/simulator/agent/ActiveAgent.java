///**
// * \package simulator.agent.zoo
// * \brief Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types
// * 
// * Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types. This package is 
// * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
// * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
// * the following URL  "http://www.cecill.info".
// */
//package simulator.agent;
//
//import java.util.ArrayList;
//
//import org.jdom.Element;
//
//import simulator.Simulator;
//import simulator.reaction.Reaction;
//import utils.ExtraMath;
//import utils.XMLParser;
//
///**
// * \brief Extension of Agent and SpecialisedAgent - adds reaction information to model agent involvement in reactions
// * 
// * Extension of Agent and SpecialisedAgent - adds reaction information to model agent involvement in reactions
// * 
// * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
// * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
// * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
// *
// */
//public abstract class ActiveAgent extends Agent {
//
//	// Parameters common to all agents of this class
//
//	// Parameters common (strict egality) to all agents of a Species
//
//	// massic growth rate of the agent (the sum of the growth rate of all of
//	// its compounds)
//
//	
//	
//	// Reaction parameters : (potentially)mutated from species parameters
//	
//	/**
//	 * \brief Creates an ActiveAgent object and initialises the object in
//	 * which associated parameters are stored.
//	 */
//	public ActiveAgent()
//	{
//		super();
//		_speciesParam = new ActiveParam();
//
//	}
//
//	/**
//	 * \brief Creates an agent of the specified species and notes the grid in
//	 * which this is assigned.
//	 * 
//	 * @param aSim	The simulation object used to simulate the conditions
//	 * specified in the protocol file.
//	 * @param xmlMarkUp	A species mark-up within the specified protocol file.
//	 */
//	@Override
//	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) 
//	{
//		// Initialisation common to all specialised agents.
//		super.initFromProtocolFile(aSim, xmlMarkUp);
//		
//		/* Create internal compounds_______________________________________ */
//		// Initialise tables for the compartments description.
//		int nParticle = aSim.particleDic.size();
//		int nReaction = aSim.reactionList.length;
//		int nSolute = aSim.soluteList.length;
//		int reacIndex;
//		
//		/*
//		 * Build the list of particles. Set the average mass of each particle
//		 * within the initial population. 
//		 */
//		particleMass = ExtraMath.newDoubleArray(nParticle);
//		int particleIndex;
//		for ( XMLParser aParser : xmlMarkUp.getChildrenParsers("particle") )
//		{
//			particleIndex = aSim.getParticleIndex(aParser.getName());
//			particleMass[particleIndex] = aParser.getParamMass("mass");
//		}
//		
//		deltaParticle = ExtraMath.newDoubleArray(particleMass.length);
//		
//		updateMass();
//		
//		/* Create description of reactions ________________________________ */
//		// Initialise the arrays.
//		allReactions = aSim.reactionList;
//		reactionKnown = new ArrayList<Integer>();
//		reactionActive = new ArrayList<Integer>();
//		growthRate = ExtraMath.newDoubleArray(nReaction);
//		soluteYield = ExtraMath.newDoubleArray(nReaction, nSolute);
//		/* 
//		 * Do not initialise reactionKinetic using ExtraMath.newDoubleArray()
//		 * as the number of j-elements in each i-array varies. Each i-array is
//		 * cloned from the reaction mark up, so no need to fill with zeros now.
//		 */
//		reactionKinetic = new Double[nReaction][];
//		
//		particleYield = ExtraMath.newDoubleArray(nReaction, nParticle);
//		
//		Reaction aReaction;
//		/*
//		 * Read the XML file. Note that this is not a parser (yet) and so we
//		 * should not use XMLParser.getName(), etc.
//		 */
//		for (Element aReacElem : xmlMarkUp.getChildrenElements("reaction"))
//		{
//			reacIndex = aSim.getReactionIndex(
//										aReacElem.getAttributeValue("name"));
//			aReaction = allReactions[reacIndex];
//			/*
//			 * Add the reaction to the list of known (and active) reactions.
//			 */
//			reactionKnown.add(reacIndex);
//			if (aReacElem.getAttributeValue("status").equals("active"))
//				reactionActive.add(reacIndex);
//			/* 
//			 * If reaction parameters have been redefined, load them; 
//			 * else load the parameters defined for the reaction.
//			 */
//			if ( aReacElem.getContentSize() == 0 )
//			{
//				soluteYield[reacIndex] = aReaction.getSoluteYield();
//				particleYield[reacIndex] = aReaction.getParticulateYield();
//				reactionKinetic[reacIndex] = aReaction.getKinetic();
//			}
//			else
//			{
//				/*
//				 * TODO Rob 15Apr2014: This seems a dangerous way of doing
//				 * things... having kinetics here that differ from those in
//				 * the simulator may play havoc with the diffusion-reaction
//				 * solver(s).
//				 */
//				aReaction.initFromAgent(this, aSim, new XMLParser(aReacElem));
//			}
//		}
//		/*
//		 * Now copy these value in the speciesParam structure.
//		 */
//		getActiveParam().soluteYield = soluteYield.clone();
//		getActiveParam().particleYield = particleYield.clone();
//		getActiveParam().reactionKinetic = reactionKinetic.clone();
//
//	}
//
//	/**
//	 * \brief Create an agent using information in a previous state or
//	 * initialisation file.
//	 * 
//	 * @param aSim	The simulation object used to simulate the conditions
//	 * specified in the protocol file.
//	 * @param singleAgentData	Data from the result or initialisation file
//	 * that is used to recreate this agent.
//	 */
//	@Override
//	public void initFromResultFile(Simulator aSim, String[] singleAgentData)
//	{
//		/*
//		 * This routine will read data from the end of the singleAgentData
//		 * array and then pass the remaining values onto the super class.
//		 * 
//		 * First find the position to start at by using length and number of
//		 * values read. 
//		 */
//		int nValsRead = 2 + particleMass.length;
//		int iDataStart = singleAgentData.length - nValsRead;
//		
//		// Read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT
//		// Particle masses:
//		for ( int iComp = 0; iComp < particleMass.length; iComp++ )
//		{
//			particleMass[iComp] = 
//						Double.parseDouble(singleAgentData[iDataStart+iComp]);
//		}
//		deltaParticle = ExtraMath.newDoubleArray(particleMass.length);
//		// Other values:
//		_netGrowthRate = Double.parseDouble(
//						singleAgentData[iDataStart+particleMass.length]);
//		_netVolumeRate = Double.parseDouble(
//						singleAgentData[iDataStart+particleMass.length+1]);
//		
//		// Now go up the hierarchy with the rest of the data.
//		String[] remainingSingleAgentData = new String[iDataStart];
//		for ( int i = 0; i < iDataStart; i++ )
//			remainingSingleAgentData[i] = singleAgentData[i];
//		super.initFromResultFile(aSim, remainingSingleAgentData);
//		
//		// Finally some creation-time calls.
//		updateSize();
//		registerBirth();		
//	}	
//	
//	/**
//	 * \brief Create a new agent with mutated parameters based on species
//	 * default values.
//	 * 
//	 * Agent is not located.
//	 */
//	@Override
//	public void createNewAgent()
//	{
//		try
//		{
//			Agent baby = (Agent) sendNewAgent();
//			// Register the baby in the pathway guilds.
//			baby.registerBirth();
//		}
//		catch (CloneNotSupportedException e)
//		{
//			System.out.println("At ActiveAgent: createNewAgent error " + e);
//		}
//	}
//
//	/**
//	 * \brief Registers a created agent into a respective container.
//	 *  
//	 * Each agent must be referenced by one such container. In this case, the 
//	 * species is registered into the agent grid.
//	 */
//	@Override
//	public void registerBirth()
//	{
//		super.registerBirth();
//		// Register the agent in the metabolic containers.
//		registerOnAllActiveReaction();
//	}
//	
//	/**
//	 * \brief Clones this agent object, creating a new progeny of this agent.
//	 * 
//	 * Ensures new clone inherits same parameters as parents.
//	 * 
//	 * @throws CloneNotSupportedException	Exception should the class not
//	 * implement Cloneable.
//	 */
//	@Override
//	@SuppressWarnings("unchecked")
//	public Object clone() throws CloneNotSupportedException
//	{
//		Agent out = (Agent) super.clone();
//		out.reactionActive = (ArrayList<Integer>) this.reactionActive.clone();
//		out.reactionKnown = (ArrayList<Integer>) this.reactionKnown.clone();
//		out.allReactions = this.allReactions.clone();
//		out.growthRate = ExtraMath.newDoubleArray(this.growthRate.length);
//		/*
//		 * No need to initialise out.soluteYield, out.reactionKinetic, or
//		 * out.particleYield using ExtraMath.newDoubleArray() as their
//		 * elements are all cloned from this agent.
//		 */
//		out.soluteYield = new Double[this.soluteYield.length][];
//		for (int iter = 0; iter < this.soluteYield.length; iter++)
//			out.soluteYield[iter] = this.soluteYield[iter].clone();
//		out.reactionKinetic = new Double[this.reactionKinetic.length][];
//		out.particleYield = new Double[this.particleYield.length][];
//		for ( int jReac : this.reactionKnown )
//		{
//			if ( this.reactionKinetic[jReac] == null )
//				out.reactionKinetic[jReac] = ExtraMath.newDoubleArray(1);
//			else
//			{
//				out.reactionKinetic[jReac] =
//										this.reactionKinetic[jReac].clone();
//			}
//			out.particleYield[jReac] = this.particleYield[jReac].clone();
//		}
//		out.particleMass = this.particleMass.clone();
//		return (Object) out;
//	}
//	
//	
//
//	/* ______________________ REACTION MANAGEMENT __________________________ */
//
//	/**
//	 * \brief Used in creation of results files - specifies the header of the columns of output information for this agent
//	 * 
//	 * Used in creation of results files - specifies the header of the columns of output information for this agent
//	 * 
//	 * @return	String specifying the header of each column of results associated with this agent
//	 */
//	@Override
//	public StringBuffer sendHeader()
//	{
////		StringBuffer tempString = super.sendHeader();
////		for (String particle : _species.currentSimulator.particleDic)
////			tempString.append(","+particle);
////		tempString.append(",growthRate,volumeRate");
////		return tempString;
//		return super.sendHeader();
//	}
//
//	/**
//	 * \brief Used in creation of results files - creates an output string of information generated on this particular agent
//	 * 
//	 * Used in creation of results files - creates an output string of information generated on this particular agent
//	 * 
//	 * @return	String containing results associated with this agent
//	 */
//	@Override
//	public StringBuffer writeOutput()
//	{
//		// write the data matching the header file
//		return super.writeOutput();
//
////		// Mass of different particles
////		for (int i = 0; i < particleMass.length; i++)
////			tempString.append(","+particleMass[i]);
////		
////		// Agent growth and volume rates
////		tempString.append(","+_netGrowthRate+","+_netVolumeRate);
////		return tempString;
//	}
//
//
//
//}
