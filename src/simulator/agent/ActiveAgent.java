/**
 * \package simulator.agent.zoo
 * \brief Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types
 * 
 * Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent;

import java.util.ArrayList;

import org.jdom.Element;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.SpatialGrid;
import simulator.reaction.Reaction;
import utils.ExtraMath;
import utils.XMLParser;

/**
 * \brief Extension of Agent and SpecialisedAgent - adds reaction information to model agent involvement in reactions
 * 
 * Extension of Agent and SpecialisedAgent - adds reaction information to model agent involvement in reactions
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 *
 */
public abstract class ActiveAgent extends SpecialisedAgent implements HasReaction {

	// Parameters common to all agents of this class

	// Parameters common (strict egality) to all agents of a Species

	// massic growth rate of the agent (the sum of the growth rate of all of
	// its compounds)

	/**
	 * Array of all the reactions this agent is involved in
	 */
	public Reaction[]            allReactions;
	
	/**
	 * Array of the reactions that are active
	 */
	protected ArrayList<Integer> reactionActive;
	
	/**
	 * Array of all reactions (even those not active on this agent)
	 */
	protected ArrayList<Integer> reactionKnown;

	/**
	 * Net growth rate of this agent due to reactions
	 */
	protected Double _netGrowthRate = 0.0;
	
	/**
	 * Net volume change rate of this agent due to reactions
	 */
	protected Double _netVolumeRate = 0.0;

	/**
	 * Growth rate of this agent due to reactions
	 */
	protected Double[] growthRate;
	
	/**
	 * The change to each particle due to reactions. Units fg.
	 */
	protected Double[] deltaParticle;
	
	// Reaction parameters : (potentially)mutated from species parameters
	
	/**
	 * Solute yield for reactions involving this agent
	 */
	public Double[][] soluteYield;
	
	/**
	 * Array of the kinetic factor parameters for the reactions this agent is involved in
	 */
	public Double[][] reactionKinetic;
	
	/**
	 * Particle yield for reactions involving this agent
	 */
	public Double[][] particleYield;

	/**
	 * Mass of the agent (table for all particles belonging to the agent)
	 */
	public Double[] particleMass;
	
	/**
	 * Sum of masses of all particles
	 */
	protected Double _totalMass;

	/**
	 * \brief Creates an ActiveAgent object and initialises the object in
	 * which associated parameters are stored.
	 */
	public ActiveAgent()
	{
		super();
		_speciesParam = new ActiveParam();

	}

	/**
	 * \brief Creates an agent of the specified species and notes the grid in
	 * which this is assigned.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param xmlMarkUp	A species mark-up within the specified protocol file.
	 */
	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) 
	{
		// Initialisation common to all specialised agents.
		super.initFromProtocolFile(aSim, xmlMarkUp);
		
		/* Create internal compounds_______________________________________ */
		// Initialise tables for the compartments description.
		int nParticle = aSim.particleDic.size();
		int nReaction = aSim.reactionList.length;
		int nSolute = aSim.soluteList.length;
		int reacIndex;
		
		/*
		 * Build the list of particles. Set the average mass of each particle
		 * within the initial population. 
		 */
		particleMass = ExtraMath.newDoubleArray(nParticle);
		int particleIndex;
		for ( XMLParser aParser : xmlMarkUp.getChildrenParsers("particle") )
		{
			particleIndex = aSim.getParticleIndex(aParser.getName());
			particleMass[particleIndex] = aParser.getParamMass("mass");
		}
		
		deltaParticle = ExtraMath.newDoubleArray(particleMass.length);
		
		updateMass();
		
		/* Create description of reactions ________________________________ */
		// Initialise the arrays.
		allReactions = aSim.reactionList;
		reactionKnown = new ArrayList<Integer>();
		reactionActive = new ArrayList<Integer>();
		growthRate = ExtraMath.newDoubleArray(nReaction);
		soluteYield = ExtraMath.newDoubleArray(nReaction, nSolute);
		/* 
		 * Do not initialise reactionKinetic using ExtraMath.newDoubleArray()
		 * as the number of j-elements in each i-array varies. Each i-array is
		 * cloned from the reaction mark up, so no need to fill with zeros now.
		 */
		reactionKinetic = new Double[nReaction][];
		
		particleYield = ExtraMath.newDoubleArray(nReaction, nParticle);
		
		Reaction aReaction;
		/*
		 * Read the XML file. Note that this is not a parser (yet) and so we
		 * should not use XMLParser.getName(), etc.
		 */
		for (Element aReacElem : xmlMarkUp.getChildrenElements("reaction"))
		{
			reacIndex = aSim.getReactionIndex(
										aReacElem.getAttributeValue("name"));
			aReaction = allReactions[reacIndex];
			/*
			 * Add the reaction to the list of known (and active) reactions.
			 */
			reactionKnown.add(reacIndex);
			if (aReacElem.getAttributeValue("status").equals("active"))
				reactionActive.add(reacIndex);
			/* 
			 * If reaction parameters have been redefined, load them; 
			 * else load the parameters defined for the reaction.
			 */
			if ( aReacElem.getContentSize() == 0 )
			{
				soluteYield[reacIndex] = aReaction.getSoluteYield();
				particleYield[reacIndex] = aReaction.getParticulateYield();
				reactionKinetic[reacIndex] = aReaction.getKinetic();
			}
			else
			{
				/*
				 * TODO Rob 15Apr2014: This seems a dangerous way of doing
				 * things... having kinetics here that differ from those in
				 * the simulator may play havoc with the diffusion-reaction
				 * solver(s).
				 */
				aReaction.initFromAgent(this, aSim, new XMLParser(aReacElem));
			}
		}
		/*
		 * Now copy these value in the speciesParam structure.
		 */
		getSpeciesParam().soluteYield = soluteYield.clone();
		getSpeciesParam().particleYield = particleYield.clone();
		getSpeciesParam().reactionKinetic = reactionKinetic.clone();

	}

	/**
	 * \brief Create an agent using information in a previous state or
	 * initialisation file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param singleAgentData	Data from the result or initialisation file
	 * that is used to recreate this agent.
	 */
	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData)
	{
		/*
		 * This routine will read data from the end of the singleAgentData
		 * array and then pass the remaining values onto the super class.
		 * 
		 * First find the position to start at by using length and number of
		 * values read. 
		 */
		int nValsRead = 2 + particleMass.length;
		int iDataStart = singleAgentData.length - nValsRead;
		
		// Read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT
		// Particle masses:
		for ( int iComp = 0; iComp < particleMass.length; iComp++ )
		{
			particleMass[iComp] = 
						Double.parseDouble(singleAgentData[iDataStart+iComp]);
		}
		deltaParticle = ExtraMath.newDoubleArray(particleMass.length);
		// Other values:
		_netGrowthRate = Double.parseDouble(
						singleAgentData[iDataStart+particleMass.length]);
		_netVolumeRate = Double.parseDouble(
						singleAgentData[iDataStart+particleMass.length+1]);
		
		// Now go up the hierarchy with the rest of the data.
		String[] remainingSingleAgentData = new String[iDataStart];
		for ( int i = 0; i < iDataStart; i++ )
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
		
		// Finally some creation-time calls.
		updateSize();
		registerBirth();		
	}	
	
	/**
	 * \brief Create a new agent with mutated parameters based on species
	 * default values.
	 * 
	 * Agent is not located.
	 */
	@Override
	public void createNewAgent()
	{
		try
		{
			ActiveAgent baby = (ActiveAgent) sendNewAgent();
			// Register the baby in the pathway guilds.
			baby.registerBirth();
		}
		catch (CloneNotSupportedException e)
		{
			System.out.println("At ActiveAgent: createNewAgent error " + e);
		}
	}

	/**
	 * \brief Registers a created agent into a respective container.
	 *  
	 * Each agent must be referenced by one such container. In this case, the 
	 * species is registered into the agent grid.
	 */
	@Override
	public void registerBirth()
	{
		super.registerBirth();
		// Register the agent in the metabolic containers.
		registerOnAllActiveReaction();
	}
	
	/**
	 * \brief Clones this agent object, creating a new progeny of this agent.
	 * 
	 * Ensures new clone inherits same parameters as parents.
	 * 
	 * @throws CloneNotSupportedException	Exception should the class not
	 * implement Cloneable.
	 */
	@Override
	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException
	{
		ActiveAgent out = (ActiveAgent) super.clone();
		out.reactionActive = (ArrayList<Integer>) this.reactionActive.clone();
		out.reactionKnown = (ArrayList<Integer>) this.reactionKnown.clone();
		out.allReactions = this.allReactions.clone();
		out.growthRate = ExtraMath.newDoubleArray(this.growthRate.length);
		/*
		 * No need to initialise out.soluteYield, out.reactionKinetic, or
		 * out.particleYield using ExtraMath.newDoubleArray() as their
		 * elements are all cloned from this agent.
		 */
		out.soluteYield = new Double[this.soluteYield.length][];
		for (int iter = 0; iter < this.soluteYield.length; iter++)
			out.soluteYield[iter] = this.soluteYield[iter].clone();
		out.reactionKinetic = new Double[this.reactionKinetic.length][];
		out.particleYield = new Double[this.particleYield.length][];
		for ( int jReac : this.reactionKnown )
		{
			if ( this.reactionKinetic[jReac] == null )
				out.reactionKinetic[jReac] = ExtraMath.newDoubleArray(1);
			else
			{
				out.reactionKinetic[jReac] =
										this.reactionKinetic[jReac].clone();
			}
			out.particleYield[jReac] = this.particleYield[jReac].clone();
		}
		out.particleMass = this.particleMass.clone();
		return (Object) out;
	}
	
	/**
	 * \brief Notifies the simulation that this agent has become too small and
	 * is then counted as dead.
	 * 
	 * Decreases the population of this species.
	 * 
	 * @param isStarving	Boolean noting whether the agent currently has
	 * access to any resources.
	 */
	@Override
	public void die(Boolean isStarving)
	{
		super.die(isStarving);
		// Unregister from the metabolic guilds.
		unregisterFromAllActiveReactions();
	}

	/**
	 * \brief Called at each time step of the simulation to compute cell
	 * growth, update size, and monitor cell death and division.
	 * 
	 * Also determines whether the agent has reached the size at which it must
	 * divide.
	 */
	@Override
	protected void internalStep()
	{
		grow();
		updateSize();
	}
	
	/**
	 * \brief Put growth rates into effect by changing the particle masses.
	 * 
	 * We adjust the particle masses after calculating all the deltaParticle
	 * values so that the reactions occur simultaneously.
	 */
	public void grow()
	{
		updateGrowthRates();
		for (int i = 0; i < particleMass.length; i++)
			particleMass[i] += deltaParticle[i];
	}
	
	/**
	 * \brief Perform agent growth by calling all active reaction pathways.
	 * 
	 */
	protected void updateGrowthRates()
	{
		Double deltaMass = 0.0;
		for (int i = 0; i < particleMass.length; i++)
			deltaParticle[i] = 0.0;
		
		int reacIndex = 0;
		// TODO RC 20 Jan 2014: Shouldn't this be the agentTimeStep?
		Double tStep = SimTimer.getCurrentTimeStep();
		Double catMass = 0.0; // Catalyst mass
		Double catYield = 0.0;
		_netGrowthRate = 0.0;
		_netVolumeRate = 0.0;
		
		// Compute mass growth rate of each active reaction
		for (int iReac = 0; iReac<reactionActive.size(); iReac++)
		{
			// Compute the growth rate
			reacIndex = reactionActive.get(iReac);
			catMass = particleMass[allReactions[reacIndex]._catalystIndex];
			// get the growth rate in [fgX.hr-1]
			growthRate[reacIndex] = allReactions[reacIndex].computeSpecGrowthRate(this);
			
			for (int i = 0; i<particleYield[reacIndex].length; i++)
			{
				if (allReactions[reacIndex].autocatalytic)
				{
					// Exponential growth/decay
					catYield = particleYield[reacIndex][allReactions[reacIndex]._catalystIndex];
					deltaParticle[i] += catMass * (particleYield[reacIndex][i]/catYield)
									    	* Math.expm1(catYield * growthRate[reacIndex]*tStep);
					deltaMass = deltaParticle[i]/tStep;
				}
				else
				{
					// Constant growth/decay
					deltaMass = catMass * particleYield[reacIndex][i]*growthRate[reacIndex];
					deltaParticle[i] += deltaMass*tStep;
				}
				_netGrowthRate += deltaMass;
				_netVolumeRate += deltaMass/getSpeciesParam().particleDensity[i];
			}
		}
		
	}
	
	/**
	 * \brief Update size of agent to take growth into account
	 * 
	 * Update mass of agent to take growth into account
	 * 
	 */
	public void updateSize() {
		updateMass();
	}

	/**
	 * \brief Update mass of agent, summing the particle mass
	 * 
	 * Update mass of agent, summing the particle mass
	 */
	public void updateMass()
	{
		_totalMass = ExtraMath.sum(particleMass);
	}

	/* ______________________ REACTION MANAGEMENT __________________________ */

	/**
	 * \brief Add the reaction to the list of known reactions
	 * 
	 * Add the reaction to the list of known reactions
	 * 
	 * @param aReaction	The reaction to add to the list
	 * @param useDefaultParam	Whether to use default reaction parameters or bespoke parameters have been specified in the protocol file
	 */
	@Override
	public void addReaction(Reaction aReaction, Boolean useDefaultParam)
	{
		// Add the reaction to the list of known reaction
		reactionKnown.add(aReaction.reactionIndex);

		// Test if specific parameters exist for this reaction
		int index = aReaction.reactionIndex;
		boolean test = getSpeciesParam().soluteYield[index]==null;

		if (useDefaultParam||test) {
			// Use parameters defined in the reaction object
			reactionKinetic[index] = aReaction.getKinetic();
			soluteYield[index] = aReaction.getSoluteYield();
			particleYield[index] = aReaction.getParticulateYield();
		} else {

			// Use parameters defined in the speciesParam structure
			reactionKinetic[index] = getSpeciesParam().reactionKinetic[index];
			soluteYield[index] = getSpeciesParam().soluteYield[index];
			particleYield[index] = getSpeciesParam().particleYield[index];
		}
	}

	/**
	 * \brief Adds an active reaction to the list of known reactions and switches the reaction on
	 * 
	 * Adds an active reaction to the list of known reactions and switches the reaction on
	 * 
	 * @param aReaction	The reaction to add to the list
	 * @param useDefaultParam	Whether to use default reaction parameters or bespoke parameters have been specified in the protocol file
	 * 
	 */
	@Override
	public void addActiveReaction(Reaction aReaction, Boolean useDefaultParam)
	{
		addReaction(aReaction, useDefaultParam);
		switchOnReaction(aReaction);
	}

	/**
	 * \brief Remove a reaction from the list of known reactions
	 * 
	 * Remove a reaction from the list of known reactions
	 * 
	 * @param aPathway	The reaction to remove from the list
	 */
	@Override
	public void removeReaction(Reaction aPathway) {
		switchOffreaction(aPathway);
		reactionKnown.remove(aPathway);
	}

	/**
	 * \brief Switches off a reaction by removing it from the active reaction array
	 * 
	 * Switches off a reaction by removing it from the active reaction array. BVM 27.11.08: added the '.reactionIndex' calls to the 
	 * two lines below in order to get this function to work correctly
	 * 
	 * @param aPathway	The reaction to switch off
	 */ 
	@Override
	public void switchOffreaction(Reaction aPathway) {
		if (reactionActive.contains(aPathway.reactionIndex)) {
			// need to remove using indexOf because the remove(object) version thinks
			// the int being passed in is the index to remove rather than the object to remove
			reactionActive.remove(reactionActive.indexOf(aPathway.reactionIndex));
			aPathway.removeAgent(this);
		}
	}

	/**
	 * \brief Switches on a reaction by adding it to the active reaction array
	 * 
	 * Switches off a reaction by adding it to the active reaction array. BVM 27.11.08: added the if statement to prevent adding a 
	 * reaction if it is already present
	 * 
	 * @param aReaction	The reaction to switch on
	 */ 
	@Override
	public void switchOnReaction(Reaction aReaction) {
		//		System.out.println("Turn it on? "+aReaction.reactionName);
		if (!reactionActive.contains(aReaction.reactionIndex)) {
			//			System.out.println("Turn on: "+aReaction.reactionName);
			reactionActive.add(aReaction.reactionIndex);
			aReaction.addAgent(this);
		}
	}

	/**
	 * \brief Register the agent on each guild of its activated pathways. Called by makeKid
	 * 
	 * Register the agent on each guild of its activated pathways. Called by makeKid
	 */
	public void registerOnAllActiveReaction() {
		for (int iReac = 0; iReac<reactionActive.size(); iReac++) {
			allReactions[reactionActive.get(iReac)].addAgent(this);
		}
	}

	/**
	 * \brief Called by the die method, to unregister the agent from all active reactions
	 * 
	 * Called by the die method, to unregister the agent from all active reactions
	 */
	public void unregisterFromAllActiveReactions() {
		for (int iReac = 0; iReac<reactionActive.size(); iReac++) {
			allReactions[reactionActive.get(iReac)].removeAgent(this);
		}
	}

	/**
	 * \brief Add the reacting concentration of an agent to the received grid
	 * 
	 * Add the reacting concentration of an agent to the received grid
	 * 
	 * @param aSpG	Spatial grid used to sum catalysing mass
	 * @param catalystIndex	Index of the compartment of the cell supporting the reaction
	 */
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex) {
	}

	/**
	 * \brief Add the reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * Add the total reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * @param aRateGrid	Spatial grid used to store total reaction rate
	 * @param reactionIndex	Index of this declared reaction in the simulation dictionary
	 */
	public void fitReacRateOnGrid(SpatialGrid aRateGrid, int reactionIndex) {
	}

	/**
	 * \brief Add the total concentration of an agent on received grid
	 * 
	 * Add the total concentration of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum total mass
	 */
	public void fitMassOnGrid(SpatialGrid aSpG) {
	}

	/**
	 * \brief Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * @return	String specifying the header of each column of results associated with this agent
	 */
	@Override
	public StringBuffer sendHeader()
	{
		StringBuffer tempString = super.sendHeader();
		for (String particle : _species.currentSimulator.particleDic)
			tempString.append(","+particle);
		tempString.append(",growthRate,volumeRate");
		return tempString;
	}

	/**
	 * \brief Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * @return	String containing results associated with this agent
	 */
	@Override
	public StringBuffer writeOutput()
	{
		// write the data matching the header file
		StringBuffer tempString = super.writeOutput();

		// Mass of different particles
		for (int i = 0; i < particleMass.length; i++)
			tempString.append(","+particleMass[i]);
		
		// Agent growth and volume rates
		tempString.append(","+_netGrowthRate+","+_netVolumeRate);
		return tempString;
	}

	/**
	 * \brief Return total mass of this agent
	 *
	 * Return total mass of this agent
	 *
	 * @return Double stating the total mass of this agent
	 */
	public Double getTotalMass() {
		return _totalMass;
	}

	/**
	 * \brief Return particle mass of this agent
	 *
	 * Return particle mass of this agent
	 *
	 * @return Double stating the particle mass of this agent
	 */
	public Double getParticleMass(int particleIndex)
	{
		return particleMass[particleIndex];
	}

	/**
	 * \brief Return net growth rate of this agent
	 *
	 * Return net growth rate of this agent
	 *
	 * @return Double stating the net growth rate of this agent
	 */
	public Double getNetGrowth()
	{
		return _netGrowthRate;
	}

	/**
	 * \brief Return volume growth rate of this agent
	 *
	 * Return volume growth rate of this agent
	 *
	 * @return Double stating the volume growth rate of this agent
	 */
	public Double getVolGrowth()
	{
		return _netVolumeRate;
	}

	/**
	 * \brief Set net growth rate of this agent
	 *
	 * Set net growth rate of this agent
	 *
	 * @param value	Double stating the net growth rate of this agent
	 */
	public void setNetGrowth(Double value)
	{
		_netGrowthRate = value;
	}

	/**
	 * \brief Return solute yield for a particular reaction
	 * 
	 * Return solute yield for a particular reaction
	 * 
	 * @param indexReaction	Index to a reaction in the simulation dictionary
	 * @return	Double array of the solute yields for that reaction
	 */
	public Double[] getSoluteYield(int indexReaction)
	{
		return soluteYield[indexReaction];
	}

	/**
	 * \brief Return the reaction kinetic parameters for a particular reaction
	 * 
	 * Return rection kinetic parameters for a particular reaction
	 * 
	 * @param indexReaction	Index to a reaction in the simulation dictionary
	 * @return	Double array of the rection kinetic parameters for that reaction
	 */
	public Double[] getReactionKinetic(int indexReaction)
	{
		return reactionKinetic[indexReaction];
	}

	/**
	 * \brief Return the parameters associated with this agent (speciesParam)
	 * 
	 * Return the parameters associated with this agent (speciesParam)
	 * 
	 * @return ActiveParam object containing the species associated with this agent
	 */
	@Override
	public ActiveParam getSpeciesParam()
	{
		return (ActiveParam) _speciesParam;
	}



}
