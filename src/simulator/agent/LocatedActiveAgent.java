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

import java.awt.Color;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.jdom.Element;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.SpatialGrid;
import simulator.geometry.ContinuousVector;
import simulator.geometry.Domain;
import simulator.geometry.boundaryConditions.AllBC;
import simulator.geometry.boundaryConditions.BoundaryCyclic;
import simulator.reaction.Reaction;
import utils.ExtraMath;
import utils.LogFile;
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
public abstract class LocatedActiveAgent extends Agent {

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
	private Double _netVolumeRate = 0.0;

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
	private Double _totalMass;

	/**
	 * Radius of this agent.
	 */
	protected Double _radius = 0.0;

	/**
	 * Cell radius including any capsules.
	 */
	protected Double _totalRadius = 0.0;

	/**
	 * Radius at which this agent will divide.
	 */
	private Double _myDivRadius = 0.0;

	/**
	 * Radius at which agent death is triggered.
	 */
	private Double _myDeathRadius = 0.0;

	/**
	 * Volume of this agent.
	 */
	protected Double _volume = 0.0;

	/**
	 * Cell volume including any capsules.
	 */
	protected Double _totalVolume = 0.0;

	/**
	 * Agent position - continuous coordinates.
	 */
	protected ContinuousVector _location = new ContinuousVector();

	/**
	 * ContinuousVector noting the move that will be applied to the agents position.
	 */
	protected ContinuousVector _movement = new ContinuousVector();

	/**
	 * Direction in which this cell divides.
	 */
	protected ContinuousVector _divisionDirection = new ContinuousVector();

	/**
	 * List of neighbouring agents in this agent's vicinity.
	 */
	protected LinkedList<Agent> _myNeighbors = new LinkedList<Agent>();

	/**
	 * Index of the agent position on the vectorized grid.
	 */
	protected int _agentGridIndex;

	/**
	 * Boolean noting whether this agent is interacting with a surface (true)
	 * or not (false).
	 */
	protected Boolean _isAttached = false;

	/**
	 * Detachment priority
	 */
	private Double detPriority = 0.0;

	/**
	 * Stores the simulation time since the last division check
	 */
	public Double _timeSinceLastDivisionCheck = Double.MAX_VALUE;

	/**
	 * Distance based probability from a given neighbour (used in HGT).
	 */
	public Double _distProb = 0.0;

	/**
	 * Cumulative probability as to whether the plasmid will be transferred.
	 */
	private Double _distCumProb = 0.0;

	/**
	 * \brief Creates an ActiveAgent object and initialises the object in
	 * which associated parameters are stored.
	 */
	public LocatedActiveAgent()
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
		getLocatedParam().soluteYield = soluteYield.clone();
		getLocatedParam().particleYield = particleYield.clone();
		getLocatedParam().reactionKinetic = reactionKinetic.clone();

		setMyDivRadius(getRandomisedDivRadius());
		setMyDeathRadius(getRandomisedDeathRadius());
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
		
		int nValsRead = 5;
		int iDataStart = singleAgentData.length - nValsRead;
		/*
		 * This is necessary for the case when agents in a biofilm
		 * simulation are transferred into a chemostat.
		 */
		if ( Simulator.isChemostat )
			_location.reset();
		else
		{
			Double newAgentX, newAgentY, newAgentZ;
			newAgentX = Double.parseDouble(singleAgentData[iDataStart]);
			newAgentY = Double.parseDouble(singleAgentData[iDataStart+1]);
			if ( _agentGrid.is3D )
				newAgentZ = Double.parseDouble(singleAgentData[iDataStart+2]);
			else
				newAgentZ = 0.0;
			_location.set(newAgentX, newAgentY, newAgentZ);
		}
		/*
		 * Agent size.
		 */
		_radius      = Double.parseDouble(singleAgentData[iDataStart+3]);
		_totalRadius = Double.parseDouble(singleAgentData[iDataStart+4]);
		/*
		 * These are randomly generated.
		 */
		setMyDivRadius(getRandomisedDivRadius());
		setMyDeathRadius(getRandomisedDeathRadius());
		/*
		 * Now go up the hierarchy with the rest of the data.
		 */
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];
		
		//-->active
		/*
		 * This routine will read data from the end of the singleAgentData
		 * array and then pass the remaining values onto the super class.
		 * 
		 * First find the position to start at by using length and number of
		 * values read. 
		 */
		nValsRead = 2 + particleMass.length;
		iDataStart = singleAgentData.length - nValsRead;
		
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
		setNetVolumeRate(Double.parseDouble(
						singleAgentData[iDataStart+particleMass.length+1]));
		
		// Now go up the hierarchy with the rest of the data.
		remainingSingleAgentData = new String[iDataStart];
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
			Agent baby = sendNewAgent();
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
		LocatedActiveAgent out = (LocatedActiveAgent) super.clone();
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
		
		//-->located
		out._location = (ContinuousVector) this._location.clone();
		out._movement = (ContinuousVector) this._movement.clone();
		out._divisionDirection = (ContinuousVector)
											this._divisionDirection.clone();
		out._myNeighbors = (LinkedList<Agent>) this._myNeighbors.clone();
		out._agentGridIndex = this._agentGridIndex;
		return out;
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
//	@Override
//	protected void internalStep()
//	{
//		grow();
//		updateSize();
//	}
	
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
	@Override
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
		setNetVolumeRate(0.0);
		
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
				setNetVolumeRate(getNetVolumeRate() + deltaMass/getLocatedParam().particleDensity[i]);
			}
		}
		
	}
	
	/**
	 * \brief Update size of agent to take growth into account
	 * 
	 * Update mass of agent to take growth into account
	 * 
	 */
//	public void updateSize() {
//		updateMass();
//	}

	/**
	 * \brief Update mass of agent, summing the particle mass
	 * 
	 * Update mass of agent, summing the particle mass
	 */
	public void updateMass()
	{
		setTotalMass(ExtraMath.sum(particleMass));
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
	public void addReaction(Reaction aReaction, Boolean useDefaultParam)
	{
		// Add the reaction to the list of known reaction
		reactionKnown.add(aReaction.reactionIndex);

		// Test if specific parameters exist for this reaction
		int index = aReaction.reactionIndex;
		boolean test = getLocatedParam().soluteYield[index]==null;

		if (useDefaultParam||test) {
			// Use parameters defined in the reaction object
			reactionKinetic[index] = aReaction.getKinetic();
			soluteYield[index] = aReaction.getSoluteYield();
			particleYield[index] = aReaction.getParticulateYield();
		} else {

			// Use parameters defined in the speciesParam structure
			reactionKinetic[index] = getLocatedParam().reactionKinetic[index];
			soluteYield[index] = getLocatedParam().soluteYield[index];
			particleYield[index] = getLocatedParam().particleYield[index];
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
		
		// location info and radius
		tempString.append(",locationX,locationY,locationZ,radius,totalRadius");
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
		tempString.append(","+_netGrowthRate+","+getNetVolumeRate());
		
		// location info and radius
		tempString.append(","+_location.x+","+_location.y+","+_location.z+",");
		tempString.append(_radius+","+_totalRadius);
		return tempString;
	}

	/**
	 * \brief Return total mass of this agent
	 *
	 * Return total mass of this agent
	 *
	 * @return Double stating the total mass of this agent
	 */
	@Override
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
	@Override
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
		return getNetVolumeRate();
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

	@Override
	protected void setTotalMass(Double _totalMass) {
		this._totalMass = _totalMass;
	}

	/**
	 * \brief Creates a daughter Located Agent by cloning this agent and
	 * parameter objects.
	 * 
	 * @throws CloneNotSupportedException Thrown if the agent cannot be cloned.
	 */
//	@Override
//	@SuppressWarnings("unchecked")
//	public Object clone() throws CloneNotSupportedException {
//		LocatedAgent out = (LocatedAgent) super.clone();
//		out._location = (ContinuousVector) this._location.clone();
//		out._movement = (ContinuousVector) this._movement.clone();
//		out._divisionDirection = (ContinuousVector)
//											this._divisionDirection.clone();
//		out._myNeighbors = (LinkedList<Agent>) this._myNeighbors.clone();
//		out._agentGridIndex = this._agentGridIndex;
//		return out;
//	}

	/**
	 * \brief Create a new agent in a specified position.
	 * 
	 * @param position	Vector stating where this agent should be located.
	 */
	@Override
	public void createNewAgent(ContinuousVector position) {
		try 
		{
			// Get a clone of the progenitor.
			Agent baby = sendNewAgent();
			baby.giveName();
			baby.updateSize();
			
			this.setMyDivRadius(getRandomisedDivRadius());
			baby.setMyDivRadius(getRandomisedDivRadius());
			baby.setMyDeathRadius(getRandomisedDeathRadius());
			
			// Just to avoid to be in the carrier.
			// TODO Rob 13Mar2015: Is this correct?
			position.x += this._totalRadius;
			
			baby.setLocation(position);
			baby.registerBirth();
		} 
		catch (CloneNotSupportedException e) 
		{
			LogFile.writeError(e, "LocatedAgent.createNewAgent()");
		}
	}

	/**
	 * \brief Creates an agent of the specified species and notes the grid in
	 * which this is assigned.
	 *
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param xmlMarkUp	A species mark-up within the specified protocol file.
	 */
//	@Override
//	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) {	
//		super.initFromProtocolFile(aSim, xmlMarkUp);
//		_myDivRadius = getDivRadius();
//		_myDeathRadius = getDeathRadius();
//	}

	/**
	 * \brief Create an agent using information in a previous state or
	 * initialization file.
	 *
	 * Reads in data from the end of the singleAgentData array and then passes
	 * the remaining values onto the super class.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param singleAgentData	Data from the result or initialisation file
	 * that is used to recreate this agent.
	 */
//	@Override
//	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
//		/*
//		 * Find the position to start at.
//		 */
//		int nValsRead = 5;
//		int iDataStart = singleAgentData.length - nValsRead;
//		/*
//		 * This is necessary for the case when agents in a biofilm
//		 * simulation are transferred into a chemostat.
//		 */
//		if ( Simulator.isChemostat )
//			_location.reset();
//		else
//		{
//			Double newAgentX, newAgentY, newAgentZ;
//			newAgentX = Double.parseDouble(singleAgentData[iDataStart]);
//			newAgentY = Double.parseDouble(singleAgentData[iDataStart+1]);
//			if ( _agentGrid.is3D )
//				newAgentZ = Double.parseDouble(singleAgentData[iDataStart+2]);
//			else
//				newAgentZ = 0.0;
//			_location.set(newAgentX, newAgentY, newAgentZ);
//		}
//		/*
//		 * Agent size.
//		 */
//		_radius      = Double.parseDouble(singleAgentData[iDataStart+3]);
//		_totalRadius = Double.parseDouble(singleAgentData[iDataStart+4]);
//		/*
//		 * These are randomly generated.
//		 */
//		_myDivRadius = getDivRadius();
//		_myDeathRadius = getDeathRadius();
//		/*
//		 * Now go up the hierarchy with the rest of the data.
//		 */
//		String[] remainingSingleAgentData = new String[iDataStart];
//		for (int i=0; i<iDataStart; i++)
//			remainingSingleAgentData[i] = singleAgentData[i];
//		super.initFromResultFile(aSim, remainingSingleAgentData);
//	}

	/**
	 * \brief Called at each time step of the simulation to compute cell
	 * growth, update size, and monitor cell death and division.
	 * 
	 * Also determines whether the agent has reached the size at which it must
	 * divide.
	 */
	@Override
	protected void internalStep() {
		/*
		 * Compute mass growth over all compartments.
		 */
		grow();
		/*
		 * Apply this mass growth of all compounds on global radius and mass.
		 */
		updateSize();
		/*
		 * Divide if you have to.
		 */
		if ( willDivide() )
			divide();
		/*
		 * Die if you have to.
		 */
		if ( willDie() )
			die(true);
	}

	/**
	 * \brief Update the radius of the agent from the current mass (and then
	 * the volume) of the agent (EPS included).
	 */
	@Override
	public void updateSize() {
		/* 
		 * Update the totalMass field (sum of the particles masses).
		 */
		updateMass();
		/*
		 * Check the mass is positive.
		 */
		if ( getTotalMass() < 0.0 )
			LogFile.writeLog("Warning: negative mass on agent "+sendName());
		/*
		 * Sum of (particles masses / particles density).
		 */
		updateVolume();
		/*
		 * Compute radius according to the volume.
		 */
		updateRadius();
		/*
		 * Check if by chance the agent is close enough to a support to be
		 * attached.
		 */
		if ( ! Simulator.isChemostat )
			updateAttachment();
	}

	/**
	 * \brief Captures cell division by making a clone of this agent using the
	 * makeKid method.
	 */
	@Override
	public void divide() {
		try
		{
			makeKid();
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "LocatedAgent.divide()");
		}
	}

	/**
	 * \brief Determines whether or not a cell has reached the radius where
	 * cell division can be triggered.
	 * 
	 * @return	Boolean stating whether cell division should be triggered
	 * (true) or not (false).
	 */
	@Override
	public boolean willDivide() {
		/*
		 * This ensures that the checks for when to divide don't occur too
		 * often; at most they will occur at the rate of AGENTTIMESTEP.
		 */
		_timeSinceLastDivisionCheck += SimTimer.getCurrentTimeStep();
		if ( _timeSinceLastDivisionCheck < _agentGrid.getAgentTimeStep() )
			return false;
		_timeSinceLastDivisionCheck = 0.0;
		/*
		 * At this point we will actually check whether to divide.
		 */
		return getRadius(false) > this.getMyDivRadius();
	}

	/**
	 * \brief Determines whether or not a cell has reached the radius where
	 * cell death can be triggered.
	 * 
	 * @return	Boolean stating whether cell death should be triggered (true)
	 * or not (false).
	 */
	@Override
	public boolean willDie() {
		return (getTotalMass() < 0.0) || (getRadius(false) <= getMyDeathRadius());
	}

	/**
	 * \brief With it determined that cell division will occur, create a new
	 * agent from the existing one.
	 * 
	 * @throws CloneNotSupportedException Thrown if the agent cannot be cloned.
	 */
	@Override
	public void makeKid() throws CloneNotSupportedException {
		/*
		 * Create the new instance.
		 */
		Agent baby = sendNewAgent();
		/*
		 * These are all generated randomly.
		 */
		this.setMyDivRadius(getRandomisedDivRadius());
		baby.setMyDivRadius(getRandomisedDivRadius());
		baby.setMyDeathRadius(getRandomisedDeathRadius());
		/*
		 * Update the lineage.
		 */
		recordGenealogy(baby);
		/*
		 * Share mass of all compounds between two daughter cells and compute
		 * new size.
		 */
		divideCompounds(baby, getBabyMassFrac());
		/*
		 * In a chemostat, the daughter cells remain with the coordinates of
		 * their progenitor. Otherwise, compute movement to apply to both
		 * cells and apply it.
		 */
		if ( ! Simulator.isChemostat )
		{
			setDivisionDirection(getInteractDistance(baby)/2);
			baby.subtractMovement(_divisionDirection);
			this.addMovement(_divisionDirection);
		}
		/*
		 * Now register the agent inside the guilds and the agent grid.
		 */
		baby.registerBirth();
		baby.setNetVolumeRate(0.0);
	}

	/**
	 * \brief On agent division, divides the mass between the old and new
	 * agent, at a specified fraction.
	 * 
	 * @param baby	The new agent, which is inheriting mass.
	 * @param babyMassFrac	The fraction of this agents mass that should be
	 * transferred to the new agent.
	 */
	public void divideCompounds(Agent baby, Double babyMassFrac) {
		/*
		 * Apply babyMassFrac.
		 */
		for (int i = 0; i<particleMass.length; i++)
		{
			baby.multiplyParticleMass(babyMassFrac,i);
			this.multiplyParticleMass(1-babyMassFrac,i);
		}
		/*
		 * Update radius, mass, volumes and growth rates.
		 */
		updateSize();
		baby.updateSize();
		updateGrowthRates();
		baby.updateGrowthRates();
	}

	/**
	 * \brief On agent division, transfers biomass and EPS between the old and
	 * new agent, at a specified ratio.
	 * 
	 * @param agent	The new agent, which is inheriting mass.
	 * @param babyMassFrac	The ratio of the biomass/EPS that should be 
	 * transferred to the new agent.
	 */
	public void transferCompounds(Agent agent, Double babyMassFrac) {
		Double[] massToTransfer = new Double[particleMass.length];
		for (int i = 0; i<particleMass.length; i++)
			massToTransfer[i] = this.particleMass[i] * babyMassFrac;
		agent.addToParticleMasses(massToTransfer);
		this.subtractFromParticleMasses(massToTransfer);
		/*
		 * Update radius, mass and volumes.
		 */
		updateSize();
		agent.updateSize();
	}

	/**
	 * \brief Set the movement vector that states where to put a newly-created
	 * particle.
	 * 
	 * @param distance	Distance between the this agent and the new agent.
	 */
	public void setDivisionDirection(Double distance) {
		Double phi, theta;
		phi = ExtraMath.getUniRandAngle();
		theta = ExtraMath.getUniRandAngle();
		_divisionDirection.x = distance * Math.sin(phi) * Math.cos(theta);
		_divisionDirection.y = distance * Math.sin(phi) * Math.sin(theta);
		if ( _agentGrid.is3D )
			_divisionDirection.z = distance * Math.cos(phi);
		else
			_divisionDirection.z = 0.0;
	}

	/**
	 * \brief Models a mechanical interaction between two located agents.
	 * 
	 * Implemented by extending classes (LocatedAgent).
	 * 
	 * @param MUTUAL	Whether movement is shared between two agents or
	 * applied only to this one. Set in the protocol file.
	 * @return	The move to be applied once the shoving or pull calculations
	 * have been performed.
	 */
	@Override
	public Double interact(boolean MUTUAL) {
		move();
		/*
		 * Rebuild your neighbourhood.
		 */
		getPotentialShovers(getInteractDistance(),_location,_radius);
		for ( Agent neighbour : _myNeighbors )
			addPushMovement(neighbour, MUTUAL);
		_myNeighbors.clear();
		return move();
	}

	/**
	 * \brief Mutual shoving : The movement by shoving of an agent is calculated based on the cell overlap and added to the agents movement vector.
	 * 
	 * Mutual shoving : The movement by shoving of an agent is calculated based on the cell overlap and added to the agents movement vector. 
	 * Both agents are moved of half the overlapping distance in opposite directions.
	 * 
	 * @param aNeighbour	 Reference to the potentially shoving neighbour
	 * @param isMutual	Whether movement is shared between two agents or applied only to this one
	 * @param gain	Double noting change in position
	 * @return Boolean stating whether shoving is detected (true) or not (false)
	 */
	public void addPushMovement(Agent neighbour, boolean isMutual) {
		/*
		 * Cannot push oneself!
		 */
		if ( neighbour == this )
			return;
		/*
		 * Find the vector from your neighbour's cell centre to your cell
		 * centre.
		 */
		ContinuousVector diff = computeDifferenceVector(neighbour);
		/*
		 * Compute effective cell-cell distance.
		 */
		Double delta = diff.norm() - getInteractDistance(neighbour);
		/*
		 * Apply the shoving calculated. If it's mutual, apply half to each.
		 */
		if ( delta < 0.0 )
		{
			diff.normalizeVector(delta);
			if ( isMutual )
			{
				diff.times(0.5);
				neighbour.addMovement(diff);
			}
			this.subtractMovement(diff);
		}
	}

	/**
	 * \brief Pulling : The movement of agents by a shrinking biofilm. Move calculated and added to the agents movement vector.
	 * 
	 * The movement of agents by a shrinking biofilm. Move calculated and added to the agents movement vector. 
	 * 
	 * TODO Not currently used... consider deleting?
	 * 
	 * @param aNeighbor	 Reference to the potentially shoving neighbour
	 * @param isMutual	Whether movement is shared between two agents or applied only to this one
	 * @param gain	Double noting change in position
	 * @return Boolean stating whether pulling is detected (true) or not (false)
	 */
	public void addSpringMovement(Agent aNeighbor, boolean isMutual) {
		/*
		 * Cannot push oneself!
		 */
		if ( aNeighbor == this )
			return;
		/*
		 * Build the escape vector and find the distance between you and your
		 * neighbour. 
		 */
		ContinuousVector diff = computeDifferenceVector(aNeighbor);
		/*
		 * Compute effective cell-cell distance. This part differs from 
		 * addPushMovement() in that the 
		 */
		Double delta = diff.norm() - getInteractDistance(aNeighbor);
		if (delta > _totalRadius)
			return;
		// TODO Rob 13Mar2015: where does this 5 come from?
		if ( delta > 0.0 )
			delta *= Math.exp(-delta * 5 / _totalRadius);
		/*
		 * Apply the shoving calculated. If it's mutual, apply half to each.
		 */
		diff.normalizeVector(delta);
		if ( isMutual )
		{
			diff.times(0.5);
			aNeighbor.addMovement(diff);
		} 
		this.subtractMovement(diff);
	}

	/**
	 * \brief Computes the shortest vector between this agent and a position
	 * given as a ContinuousVector. Assumes cyclic boundaries.
	 * 
	 * If the vector is all zero's, returns a vector of random direction and
	 * length = 0.01 * radius.
	 * 
	 * TODO Can we do this without assuming cyclic boundaries? I.e. actually
	 * check..
	 * 
	 * @param neighbour	ContinuousVector of position to calculate distance to.
	 * @return The shortest movement vector to go from a to b, taking into
	 * account the cyclic boundary.
	 * @see addOverlapMovement
	 * @see addPullMovement works in 2 and 3D
	 */
	public ContinuousVector computeDifferenceVector(ContinuousVector position) {
		Double gridLength;
		ContinuousVector diff = new ContinuousVector();
		diff.sendDiff(_location, position);
		/*
		 * Check periodicity in X.
		 */
		gridLength = _species.domain.length_X;
		if ( Math.abs(diff.x) > 0.5 * gridLength )
			diff.x -= Math.signum(diff.x) * gridLength;
		/*
		 * Check periodicity in Y.
		 */
		gridLength = _species.domain.length_Y;
		if ( Math.abs(diff.y) > 0.5 * gridLength )
			diff.y -= Math.signum(diff.y) * gridLength;
		/*
		 * Check periodicity in Z.
		 */
		if (_agentGrid.is3D)
		{
			gridLength = _species.domain.length_Z;
			if (Math.abs(diff.z) > 0.5 * gridLength)
				diff.z -= Math.signum(diff.z) * gridLength;
		}
		/*
		 * If this is a zero vector, give it random direction and a norm of
		 * 0.01 * radius.
		 */
		if ( diff.isZero() )
		{
			diff.alea(_agentGrid.is3D);
			diff.normalizeVector(0.01*_radius);
		}
		return diff;
	}

	/**
	 * 
	 * @param aLoc
	 * @return
	 */
	public ContinuousVector computeDifferenceVector(Agent aLoc) {
		return computeDifferenceVector(aLoc.getLocation());
	}

	/**
	 * \brief Find neighbouring agents in a range around you.
	 * 
	 * @param radius	The distance to search around the agent location.
	 */
	public void getPotentialShovers(Double radius) {
		_agentGrid.getPotentialShovers(_agentGridIndex, radius, _myNeighbors);
	}

	public void getPotentialShovers(Double interactDistance, ContinuousVector location, Double radius) {
		List<Agent> tempAgentList = _agentGrid.neighborhoodSearch(getSearchCoord(radius+interactDistance),(radius+interactDistance)*2);
			
		
		for (Agent a: tempAgentList)
			_myNeighbors.add(a);
		
	}

	/**
	 * \brief Pick a random neighbour from the _myNeigbors collection.
	 * 
	 * @return	A randomly picked neighbour (LocatedAgent object) from the
	 * list of neighbours.
	 */
	public Agent pickNeighbor() {
		if (_myNeighbors.isEmpty())
			return null;
		return _myNeighbors.get(ExtraMath.getUniRandInt(_myNeighbors.size()));
	}

	/**
	 * \brief Find siblings of this agent in the immediate surroundings.
	 * 
	 * @param indexSpecies	The index used to reference this species in the
	 * simulation dictionary.
	 */
	public void findCloseSiblings(int indexSpecies) {
		Double shoveDist;
		Agent aNb;
		/*
		 * Find and count neighbours.
		 */
		getPotentialShovers(getInteractDistance(),_location,_radius);
		int nNb = _myNeighbors.size();
		/*
		 * Loop through them, only re-appending them to the neighbour list
		 * if they are: (1) different to this agent, (2) the same species as 
		 * this agent, and (3) close enough to this agent.  
		 */
		for ( int iNb = 0; iNb < nNb; iNb++ )
		{
			aNb = _myNeighbors.removeFirst();
			if ( aNb == this || indexSpecies != aNb.speciesIndex)
				continue;
			shoveDist = 2 * (getShoveRadius() + aNb.getShoveRadius());
			if ( getDistance(aNb) <= shoveDist )
				_myNeighbors.addLast(aNb);
		}
	}

	/**
	 * \brief With the agent move calculated, apply this movement, taking care
	 * to respect boundary conditions.
	 * 
	 * @return Distance moved relative to total radius.
	 */
	@Override
	public Double move() {
		/*
		 * Check the movement is valid.
		 */
		if ( ! _movement.isValid() )
		{
			LogFile.writeLog("Incorrect movement coordinates");
			_movement.reset();
		}
		/*
		 * Check we're not trying to move in the Z direction in 2D.
		 */
		if ( !(_agentGrid.is3D) && !(_movement.z.equals(0.0)) )
		{
			_movement.z = 0.0;
			_movement.reset();
			LogFile.writeLog("Agent tried to move in Z direction!");
		}
		/*
		 * No movement planned, finish here.
		 */
		if (_movement.isZero())
			return 0.0;

		/*
		 * Now apply the movement.
		 */
		setLocation(getVerifiedMovement(_movement));
		_agentGrid.registerMove(this);
		/*
		 * Calculate how far we've traveled relative to the total radius.
		 */
		Double delta = _movement.norm();
		_movement.reset();
		return delta/_totalRadius;
	}

	/**
	 * \brief Used by the move method to determine if an agent's move crosses
	 * any of the domain's boundaries.
	 */
	public ContinuousVector getVerifiedMovement(ContinuousVector movement) {
		// Search a boundary which will be crossed
		ContinuousVector newLoc = new ContinuousVector(_location);
		newLoc.add(movement);
		return getVerifiedLocation(newLoc);
	}
	
	public ContinuousVector getVerifiedLocation(ContinuousVector location) {
		AllBC aBoundary = getDomain().testCrossedBoundary(this,location);
		int nDim = (_agentGrid.is3D ? 3 : 2);
		Boolean test = ( aBoundary != null );
		int counter = 0;
		/*
		 * Test all boundaries and apply corrections according to crossed
		 * boundaries.
		 */
		while (test)
		{
			counter++;
			aBoundary.applyBoundary(this, location);
			aBoundary = getDomain().testCrossedBoundary(this,location);
			test = (aBoundary != null) || (counter > nDim);
			// TODO Rob 16Mar2015: Not sure why iDynoMiCS is failing at cyclic
			// boundaries.
//			test = test && !(aBoundary instanceof BoundaryCyclic);
			if (counter > nDim)
			{
				LogFile.writeLogAlways(
						"Problem in LocatedAgent.checkBoundaries(): "+
						"\n\tLocatedAgent at "+_location.toString()+
						", trying to move to "+location.toString()+
						", with radius "+_totalRadius+" on boundary "+
						aBoundary.getSide()+" ("+aBoundary.getClass()+")"+
						"\n\tcounter ("+counter+") > nDim ("+nDim+")");
				return _location;
			}
		}
		return location;
	}

	/**
	 * \brief Add the reacting concentration of an agent to the received grid
	 * 
	 * Add the reacting concentration of an agent to the received grid
	 * 
	 * @param aSpG	Spatial grid used to sum catalysing mass
	 * @param catalystIndex	Index of the compartment of the cell supporting the reaction
	 */
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex) //FIXME: Bas: method is misleading this method does not fit a mass but a 
	// 'concentration' on a grid (mass / grid volume)
	{
		if (isDead)
			return;
	
		Double value = particleMass[catalystIndex]/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		aSpG.addValueAt(value, _location);
	}

	/**
	 * \brief Add the total concentration of an agent on received grid
	 * 
	 * Add the total concentration of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum total mass
	 */
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG) //FIXME: Bas: method is misleading this method does not fit a mass but a 
	// 'concentration' on a grid (mass / grid volume)
	{
		if (isDead)
			return;
	
		Double value = getTotalMass()/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		aSpG.addValueAt(value, _location);
	}

	/**
	 * \brief Add the total volume rate of an agent on received grid
	 * 
	 * Add the total volume rate of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum volume
	 */
	@Override
	public void fitVolRateOnGrid(SpatialGrid aSpG) {
		Double value = getNetVolumeRate()/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		try
		{
			aSpG.addValueAt(value, _location);
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			LogFile.writeLogAlways("Could not put LocatedAgent mass on grid");
			LogFile.writeLogAlways("Problem with location "
													+_location.toString());
			System.exit(-1);
		}
	}

	/**
	 * \brief Add the reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * Add the total reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * @param aRateGrid	Spatial grid used to store total reaction rate
	 * @param reactionIndex	Index of this declared reaction in the simulation dictionary
	 */
	@Override
	public void fitReacRateOnGrid(SpatialGrid aRateGrid, int reactionIndex) {
		if (isDead)
			return;
		
		// growthRate is in [fgX.hr-1] so convert to concentration:
		// [fgX.um-3.hr-1 = gX.L-1.hr-1]
		Double value = growthRate[reactionIndex]/aRateGrid.getVoxelVolume();
	
		if ( ! Double.isFinite(value) )
			value = 0.0;
	
		aRateGrid.addValueAt(value, _location);
	}

	/**
	 * \brief Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * @return	String specifying the header of each column of results associated with this agent
	 */
//	@Override
//	public StringBuffer sendHeader() {
//		// return the header file for this agent's values after sending those for super
//		StringBuffer tempString = super.sendHeader();
//		
//		// location info and radius
//		tempString.append(",locationX,locationY,locationZ,radius,totalRadius");
//		
//		return tempString;
//	}

	/**
	 * \brief Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * @return	String containing results associated with this agent
	 */
//	@Override
//	public StringBuffer writeOutput() {
//		// write the data matching the header file
//		StringBuffer tempString = super.writeOutput();
//		
//		// location info and radius
//		tempString.append(","+_location.x+","+_location.y+","+_location.z+",");
//		tempString.append(_radius+","+_totalRadius);
//		
//		return tempString;
//	}

	/**
	 * \brief Compute the volume on the basis of the mass and density of different compounds defined in the cell
	 * 
	 * Compute the volume on the basis of the mass and density of different compounds defined in the cell
	 */
	public void updateVolume() {
		_volume = 0.0;
		for (int i = 0; i<particleMass.length; i++) {
			_volume += particleMass[i]/getLocatedParam().particleDensity[i];
		}
		_totalVolume = _volume;
	}

	/**
	 * \brief Compute the radius on the basis of the volume The radius evolution is stored in deltaRadius (used for shrinking)
	 * 
	 * Compute the radius on the basis of the volume The radius evolution is stored in deltaRadius (used for shrinking)
	 */
	public void updateRadius() {
	
		//sonia:chemostat 22.02.2010
		if(Simulator.isChemostat || _species.domain.is3D){
			_radius = ExtraMath.radiusOfASphere(_volume);
			_totalRadius = ExtraMath.radiusOfASphere(_totalVolume);
		}else{
			_radius = ExtraMath.radiusOfACylinder(_volume,
					_species.domain.length_Z);
			_totalRadius = ExtraMath.radiusOfACylinder(_totalVolume,
					_species.domain.length_Z);
		}
	}

	/**
	 * \brief Update the attachment, checking if this agent is close enough to
	 * any boundaries.
	 * 
	 * TODO Rob 13Mar2015: Where does this 3 come from?!
	 * 
	 * @return	Boundary that has been crossed.
	 */
	public void updateAttachment() {
		for (AllBC aBoundary : getDomain().getAllBoundaries())
			if ( aBoundary.isSupport() &&
						aBoundary.getDistance(_location) <= 3 * _totalRadius )
			{
				_isAttached = true;
				return;
			}
	}

	/**
	 * \brief Add movement to the ContinuousVector storing the agents move.
	 * 
	 * @param aMove	ContinuousVector to add to the movement vector.
	 */
	@Override
	public void addMovement(ContinuousVector aMove) {
		this._movement.add(aMove);
	}
	
	@Override
	public void subtractMovement(ContinuousVector aMove) {
		this._movement.subtract(aMove);
	}

	/**
	 * \brief Return the set of parameters associated with this agent
	 * (LocatedParam object).
	 * 
	 * @return LocatedParam object of parameters associated with this agent.
	 */
	@Override
	public LocatedParam getLocatedParam() {
		return (LocatedParam) _speciesParam;
	}
	
	@Override
	public SpeciesParam getSpeciesParam() {
		return _speciesParam;
	}

	/**
	 * \brief Return the volume of this agent, with or without the capsule.
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the volume of this agent.
	 */
	@Override
	public Double getVolume(boolean withCapsule) {
		return withCapsule ? _totalVolume : _volume;
	}

	/**
	 * \brief Return the radius of this agent, with or without the capsule.
	 * 
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the radius of this agent
	 */
	@Override
	public Double getRadius(boolean withCapsule) {
		return (withCapsule ? _totalRadius : _radius);
	}

	/**
	 * \brief Return the mass of this agent, with or without the capsule.
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the mass of this agent.
	 */
	public Double getMass(Boolean withCapsule) {
		return (withCapsule ? getTotalMass() : getTotalMass());
	}

	/**
	 * \brief Report whether this cell has any EPS.
	 * 
	 * @return	Boolean noting whether this cell has any EPS.
	 */
	@Override
	public Boolean hasEPS() {
		return false;
	}

	/**
	 * \brief Return the shove factor to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove factor that will be applied.
	 */
	public Double getShoveFactor() {
		return ((LocatedParam) _speciesParam).shoveFactor;
	}

	/**
	 * \brief Return the shove radius to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove radius that will be applied.
	 */
	@Override
	public Double getShoveRadius() {
		return _totalRadius * getShoveFactor();
	}

	/**
	 * \brief Return the shove limit to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove limit that will be applied.
	 */
	public Double getShoveLimit() {
		return ((LocatedParam) _speciesParam).shoveLimit;
	}

	/**
	 * \brief Return the shoving interaction distance to be used in shoving
	 * for this species of agent.
	 * 
	 * @return	Double specifying the shoving interaction distance that will
	 * be applied.
	 */
	public Double getInteractDistance() {
		return getInteractDistance(this);
	}

	/**
	 * \brief Return the shoving interaction distance to be used in shoving
	 * against a specified agent.
	 * 
	 * @return	Double specifying the shoving interaction distance that will
	 * be applied.
	 */
	public Double getInteractDistance(Agent neighbour) {
		return getShoveRadius() + neighbour.getShoveRadius() + getShoveLimit();
	}

	/**
	 * \brief Return the fraction of mass that is transferred to the new agent
	 * on cell division.
	 * 
	 * @return	Double stating the fraction of mass that is transferred to the
	 * new agent on cell division.
	 */
	public Double getBabyMassFrac() {
		return ExtraMath.deviateFromCV(getLocatedParam().babyMassFrac,
											getLocatedParam().babyMassFracCV);
	}

	/**
	 * \brief Return the agent radius at which cell division is triggered.
	 * 
	 * @return	Double stating the agent radius at which cell division is
	 * triggered.
	 */
	public Double getRandomisedDivRadius() {
		return ExtraMath.deviateFromCV(getLocatedParam().divRadius,
											getLocatedParam().divRadiusCV);
	}

	/**
	 * \brief Return the agent radius at which cell death is triggered
	 * 
	 * @return	Double stating the agent radius at which cell death is triggered
	 */
	public Double getRandomisedDeathRadius() {
		return ExtraMath.deviateFromCV(getLocatedParam().deathRadius,
											getLocatedParam().deathRadiusCV);
	}

	/**
	 * \brief Report if this agent is attached to a surface.
	 * 
	 * @return Boolean noting whether the agent is attached to a surface.
	 */
	public Boolean isAttached() {
		return _isAttached;
	}

	/**
	 * \brief Return the active fraction of this agent.
	 * 
	 * @return	Double value stating the active fraction of this agent.
	 */
	@Override
	public Double getActiveFrac() {
		return 1.0;
	}

	/**
	 * \brief Return the color assigned to this agent in POV-Ray output.
	 * 
	 * TODO Rob 13Mar2015: Consider deleting as part of move away from POV-Ray.
	 * 
	 * @return	Colour assigned to this agent as specified in the protocol
	 * file.
	 */
	@Override
	public Color getColor() {
		return _species.color;
	}

	/**
	 * \brief Return the colour assigned to any capsules contained in this
	 * agent in POV-Ray output.
	 * 
	 * @return	Colour assigned to this agent capsules as specified in the
	 * protocol file.
	 */
	@Override
	public Color getColorCapsule() {
		return Color.green;
	}

	/**
	 * \brief Return the location of this agent.
	 * 
	 * @return	ContinuousVector stating the location of this agent.
	 */
	@Override
	public ContinuousVector getLocation() {
		return _location;
	}

	/**
	 * \brief Return the distance from this agent to a ContinuousVector.
	 * 
	 * @param position	ContinuousVector to find distance to.
	 * @return distance between this agent and cV (assuming cyclic boundaries).
	 */
	public Double getDistance(ContinuousVector position) {
		return computeDifferenceVector(position).norm();
	}

	/**
	 * \brief Return the distance from this agent to another.
	 * 
	 * @param aLoc	LocatedAgent to find distance to.
	 * @return Distance from this agent to that given (assuming cyclic
	 * boundaries).
	 */
	public Double getDistance(Agent aLoc) {
		return getDistance(aLoc.getLocation());
	}

	@Override
	public Double[] getParticleMass() {
		return particleMass;
	}

	@Override
	public void addParticleMass(Double mass, int particleIndex) {
			particleMass[particleIndex] += mass;
	}

	@Override
	public void multiplyParticleMass(Double multiplier, int particleIndex) {
			particleMass[particleIndex] *= multiplier;
	}

	@Override
	public void addToParticleMasses(Double[] mass) {
		for (int i=0; i < particleMass.length; i++)
			particleMass[i] += mass[i];
	}

	@Override
	public void subtractFromParticleMasses(Double[] mass) {
		for (int i=0; i < particleMass.length; i++)
			particleMass[i] -= mass[i];
	}

	@Override
	public void setDetPriority(Double priority) {
		this.detPriority = priority;
	}

	protected Double getMyDivRadius() {
		return _myDivRadius;
	}

	@Override
	protected void setMyDivRadius(Double _myDivRadius) {
		this._myDivRadius = _myDivRadius;
	}

	@Override
	public void addDetPriority(Double priority) {
		setDetPriority(getDetPriority() + priority);
	}

	@Override
	public void multiplyDetPriority(Double multiplier) {
		setDetPriority(getDetPriority() * multiplier);
	}

	@Override
	public void setDistProb(Double prob) {
		this._distProb = prob;
	}

	@Override
	public Double getDistProb() {
		return _distProb;
	}

	@Override
	public Double getDistCumProb() {
		return _distCumProb;
	}

	@Override
	public void setDistCumProb(Double _distCumProb) {
		this._distCumProb = _distCumProb;
	}

	@Override
	public Double getDetPriority() {
		return detPriority;
	}

	/**
	 * \brief Set the location of this agent to the supplied continuous vector.
	 * 
	 * @param cc	Location which this agent should be assigned to.
	 */
	@Override
	public void setLocation(ContinuousVector cc) {
		// In a chemostat set the location of the newborns to zero.
		if ( Simulator.isChemostat )
			_location.reset();
		else
			_location.set(cc);
	}

	/**
	 * \brief Return the continuous vector that states this agents move.
	 * 
	 * @return Continuous vector that states this agents move.
	 */
	@Override
	public ContinuousVector getMovement() {
		return _movement;
	}

	/**
	 * \brief Return the index of the grid on which this agent is placed.
	 * 
	 * @return Integer grid index of where this agent is placed.
	 */
	@Override
	public int getGridIndex() {
		return _agentGridIndex;
	}

	/**
	 * \brief Return the LocatedGroup of agents that are present in the
	 * location where this agent is placed.
	 * 
	 * @return	LocatedGroup containing all agents present in the same grid
	 * space as this agent.
	 */
	public LocatedGroup getGridElement() {
		return _agentGrid.getShovingGrid()[_agentGridIndex];
	}

	/**
	 * \brief Move this agent to another grid index.
	 * 
	 * @param aGridIndex Grid index in which this agent should now be placed.
	 */
	@Override
	public void setGridIndex(int aGridIndex) {
		_agentGridIndex = aGridIndex;
	}

	/**
	 * \brief Return the domain where this agent is contained.
	 * 
	 * @return The domain where this agent is contained (Domain object).
	 */
	public Domain getDomain() {
		return _species.domain;
	}

	@Override
	public double[] getBoundingBoxCoord() {
		int dim = 2;
		if (_agentGrid.is3D)
			dim = 3;
		double[] coord = new double[dim];
		for (int i = 0; i < dim; i++) {
			coord[i] = ((_location.get()[i] )-_radius);
		}
		return coord;
		 
	}

	@Override
	public double[] getSearchCoord(Double range) {
		int dim = 2;
		if (_agentGrid.is3D)
			dim = 3;
		double[] coord = new double[dim];
		for (int i = 0; i < dim; i++) {
			coord[i] = ((_location.get()[i] )-range);
		}
		return coord;
		 
	}

	@Override
	public double[] getBoundingBoxDimensions() {
		int dim = 2;
		if (_agentGrid.is3D)
			dim = 3;
		double[] dimensions = new double[dim];
		for (int i = 0; i < dim; i++) {
			dimensions[i] = (_radius*2);
		}
		return dimensions;
	}

	protected Double getMyDeathRadius() {
		return _myDeathRadius;
	}

	@Override
	protected void setMyDeathRadius(Double _myDeathRadius) {
		this._myDeathRadius = _myDeathRadius;
	}

	protected Double getNetVolumeRate() {
		return _netVolumeRate;
	}

	@Override
	protected void setNetVolumeRate(Double _netVolumeRate) {
		this._netVolumeRate = _netVolumeRate;
	}

}
