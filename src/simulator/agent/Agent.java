/**
 * \package agent
 * \brief Package of utilities that create and manage agents in the simulation
 * and their participation in relevant reactions.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent;

import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;

import org.jdom.Element;

import idyno.SimTimer;
import simulator.AgentContainer;
import simulator.Simulator;
import simulator.SpatialGrid;
import simulator.reaction.Reaction;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Major class of iDynoMiCS - defines the agents that are involved in
 * an iDynoMiCS simulation.
 * 
 * Extended by a number of agent types.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany).
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 */
public abstract class Agent implements Cloneable
{
	/* Parameters common to all agents of this class ________________________ */
	/* Temporary variables stored in static fields __________________________ */
	/* Parameters common (strict equality) to all agents of a Species _________ */
	/* Parameters mutated from species parameters ___________________________ */

	/**
	 * Integer noting the last simulation timestep when this agent was stepped
	 */
	protected int _lastStep;

	/**
	 * The number of generations between the progenitor and the current agent
	 */
	protected int _generation = 0;
	
	/**
	 * Integer for the binary reading of the 0 and 1 coding the lineage. When
	 * a cells divides, one daughter has the index value 1, the other the
	 * index value 0, then this index is added on the left of the lineage
	 * description.
	 */
	protected BigInteger _genealogy  = BigInteger.ZERO;
	
	/**
	 * Integer noting the family which this agent belongs to
	 */
	protected int        _family     = 0;
	
	/**
	 * Integer noting the next family that any newly created agent will belong
	 * to.
	 */
	protected static int nextFamily  = 0;
	
	/**
	 * Time at which this agent was created.
	 */
	protected Double _birthday;

	public boolean isDead;

	public String death;

	/**
	 * Type of species that this agent is representing
	 */
	protected Species _species;

	/**
	 * Integer index to that species within the simulation dictionary
	 */
	public int speciesIndex;

	/**
	 * Set of parameters associated with this specialised agent
	 */
	protected SpeciesParam _speciesParam;

	/**
	 * Grid in which this agent is contained
	 */
	protected AgentContainer _agentGrid;

	/**
	 * Array of all the reactions this agent is involved in
	 */
	public Reaction[] allReactions;

	/**
	 * Array of the reactions that are active
	 */
	public ArrayList<Integer> reactionActive;

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
	 * \brief Initialise an agent object, setting its time of creation and
	 * thus the time it was last stepped.
	 */
	public Agent()
	{
		_birthday = SimTimer.getCurrentTime();
		_lastStep = SimTimer.getCurrentIter()-1;
	}
	
	/**
	 * \brief Initialise the agent from the protocol file.
	 * 
	 * Implemented by classes that extend this class.
	 * 
	 * @param aSimulator	The simulation object used to simulate the
	 * conditions specified in the protocol file.
	 * @param aSpeciesRoot	A Species mark-up within the specified protocol
	 * file.
	 */
//	public void initFromProtocolFile(Simulator aSimulator,
//													XMLParser aSpeciesRoot)
//	{
//		
//	}
	
	/**
	 * \brief Create an agent using information in a previous state or
	 * initialisation file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param singleAgentData	Data from the result or initialisation file 
	 * that is used to recreate this agent.
	 */
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) 
	{

		
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


		// read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT
		_family     = Integer.parseInt(singleAgentData[0]);
		_genealogy  = new BigInteger(singleAgentData[1]);
		_generation = Integer.parseInt(singleAgentData[2]);
		_birthday   = Double.parseDouble(singleAgentData[3]);
		
		// Finally some creation-time calls.
		updateSize();
		registerBirth();	
	}
	
	/**
	 * \brief Creates a new agent from an existing one, and registers this new
	 * agent in the simulation.
	 * 
	 * @throws CloneNotSupportedException	Exception should the class not
	 * implement Cloneable.
	 */
	public void makeKid() throws CloneNotSupportedException 
	{

		this.clone();
		// Now register the agent in the appropriate container
		registerBirth();
	}

	/**
	 * \brief Clones this agent object, creating a new progeny of this agent.
	 * 
	 * @throws CloneNotSupportedException	Exception should the class not
	 * implement Cloneable.
	 */
	@Override
	public Object clone() throws CloneNotSupportedException
	{
		Agent out = (Agent) super.clone();

		// Copy the references (superficial copy)
		out._species = this._species;
		out._speciesParam = this._speciesParam;

		//from active
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
	 * \brief Registers a created agent into a respective container.
	 * 
	 * Each agent must be referenced by one such container. Implemented by
	 * classes that extend Agent.
	 */
//	public abstract void registerBirth();
	
	/**
	 * \brief Perform the next timestep of the simulation for this agent.
	 * 
	 * _lastStep is implemented to note that the agent has been stepped.
	 * Implemented fully by agent types that extend Agent.
	 */
	public void step()
	{
		_lastStep = SimTimer.getCurrentIter();
		internalStep();
	}
	
	/**
	 * \brief Called at each time step (under the control of the method Step
	 * of the class Agent to avoid multiple calls).
	 * 
	 * Implemented by classes that extend Agent.
	 */
//	protected abstract void internalStep();
	
	/**
	 * \brief Specifies the header of the columns of output information for
	 * this agent.
	 * 
	 * Used in creation of results files.
	 * 
	 * @return	String specifying the header of each column of results
	 * associated with this agent.
	 */
	public StringBuffer sendHeader()
	{
		StringBuffer tempString = new StringBuffer("family,genealogy,generation,birthday");
		for (String particle : _species.currentSimulator.particleDic)
			tempString.append(","+particle);
		tempString.append(",growthRate,volumeRate");
		return tempString;
	}
	
	/**
	 * \brief Creates an output string of information generated on this
	 * particular agent.
	 * 
	 * Used in creation of results files. Data written matches the headers in
	 * sendHeader().
	 * 
	 * @return	String containing results associated with this agent.
	 */
	public StringBuffer writeOutput()
	{
		StringBuffer tempString = new StringBuffer(
						_family+","+_genealogy+","+_generation+","+_birthday);
		
		// Mass of different particles
		for (int i = 0; i < particleMass.length; i++)
			tempString.append(","+particleMass[i]);
		
		// Agent growth and volume rates
		tempString.append(","+_netGrowthRate+","+_netVolumeRate);
		return tempString;
	}
	
	/**
	 * \brief Called when creating an agent : updates _generation and
	 * _genealogy field.
	 * 
	 * @param baby The newly created agent that is the next generation of this
	 * agent.
	 */
	protected void recordGenealogy(Agent baby) 
	{
		// Rob 18/1/11: Shuffled around slightly to include odd numbers
		//baby._genealogy = _genealogy+ExtraMath.exp2long(this._generation);
		// Rob 11/2/15: Changed to BigInteger
		baby._genealogy = new BigInteger("2");
		baby._genealogy.pow(this._generation);
		baby._genealogy.add(this._genealogy);
		
		this._generation++;
		baby._generation = this._generation;
		
		// We want to know if the genealogy goes negative.
		if ( baby._genealogy.signum() < 0 )
		{
			LogFile.writeLog("Warning: baby's genealogy has gone negative:");
			LogFile.writeLog("family "+baby._family+", genealogy "+
							baby._genealogy+", generation "+baby._generation);
		}
		
		// Rob 21/1/11: changed so that only the baby is given a new birthday
		// this._birthday = SimTimer.getCurrentTime();
		// baby._birthday = this._birthday;
		baby._birthday = SimTimer.getCurrentTime();
	}

	/**
	 * \brief Returns a string containing the family name and genealogy of
	 * this agent.
	 * 
	 * @return	String containing the family name and genealogy of this agent.
	 */
	public String sendName()
	{
		return _family+"-"+_genealogy;
	}
	
	/**
	 * \brief Set the family for this agent, based on the next family.
	 */
	public void giveName() 
	{
		_family = ++nextFamily;
	}
	
	/**
	 * \brief Return the simulation time at which this agent was created.
	 * 
	 * @return	Double noting the simulation time at which this agent was
	 * created.
	 */
	public Double getBirthday()
	{
		return this._birthday;
	}

	public Double move() {
		return null;
		// TODO Auto-generated method stub
		
	}

	public Double interact(boolean mUTUAL) {
		// TODO Auto-generated method stub
		return null;
	}

	public void die(boolean b) {

		// Unregister from the metabolic guilds.
		unregisterFromAllActiveReactions();
		
	}


	public void fitVolRateOnGrid(SpatialGrid biomassGrid) {
		// TODO Auto-generated method stub
		
	}


	/**
	 * \brief Returns the object containing a set of parameters associated with a particular agent (species)
	 * 
	 * Returns the object containing a set of parameters associated with a particular agent (species)
	 */
	public SpeciesParam getSpeciesParam() {
		return _speciesParam;
	}

	/**
	 * \brief Returns the species object that is represented by this agent.
	 * 
	 * @return	Object of the Species class that this agent is representing.
	 */
	public Species getSpecies() {
		return _species;
	}

	/**
	 * \brief Returns the species object that is represented by this agent.
	 * 
	 * @return	Object of the Species class that this agent is representing.
	 */
	public int getSpeciesIndex() {
		return _species.speciesIndex;
	}

	/**
	 * \brief Set the progenitor Specialised agent to a specified species.
	 * 
	 * @param aSpecies	A species object to use as the progenitor.
	 */
	public void setSpecies(Species aSpecies) {
		_species = aSpecies;
		speciesIndex = aSpecies.speciesIndex;
	}

	/**
	 * \brief Return the name of the species represented by this agent.
	 * 
	 * @return	Name of the species represented.
	 */
	public String getName() {
		return _species.speciesName;
	}

	/**
	 * \brief Creates an agent of the specified species and notes the grid in which this is assigned
	 *
	 * Creates an agent of the specified species and notes the grid in which this is assigned
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A species mark-up within the specified protocol file
	 */

	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot) {
		try 
		{
			// Create the agent object
			_agentGrid = aSim.agentGrid;
		} 
		catch (Exception e) 
		{
			LogFile.writeLog("Creating "+this.getSpecies().speciesName);
			System.exit(-1);
		}
		
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
		for ( XMLParser aParser : aSpeciesRoot.getChildrenParsers("particle") )
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
		for (Element aReacElem : aSpeciesRoot.getChildrenElements("reaction"))
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
		getActiveParam().soluteYield = soluteYield.clone();
		getActiveParam().particleYield = particleYield.clone();
		getActiveParam().reactionKinetic = reactionKinetic.clone();

	}

	/**
	 * \brief Registers a created agent into a respective container.
	 *  
	 * Each agent must be referenced by one such container. In this case, the 
	 * species is registered into the agent grid.
	 */
	public void registerBirth() {
		_agentGrid = _species.currentSimulator.agentGrid;
		_agentGrid.registerBirth(this);
		_species.notifyBirth();
		// Register the agent in the metabolic containers.
				registerOnAllActiveReaction();
	}

	/**
	 * \brief Notifies the simulation that this agent has become too small and
	 * is then counted as dead.
	 * 
	 * @param isStarving	Boolean noting whether the agent currently has
	 * access to any resources.
	 */
	public void die(Boolean isStarving) {
		/*
		 * If you are too small, you must die!
		 * Decrease the population of your species
		 */
		_species.notifyDeath();
		isDead = true;
		_agentGrid.registerDeath(this);
		unregisterFromAllActiveReactions();
	}

	/**
	 * \brief Used in POV-Ray output to display this species - writes a colour definition to the passed-in file
	 * 
	 * Used in POV-Ray output to display this species. This writes a color definition to the passed-in file. Meant for later use in 
	 * macros. Note that this routine is put here and not in Species to allow derived agents to use different colors for different 
	 * states; EpiBac is one example, with different colors for donor, recipient, and transconjugant states
	 * 
	 * @param fr	POV-Ray output file where the colour definition should be applied
	 */
	public void writePOVColorDefinition(FileWriter fr) throws IOException {
		fr.write("#declare "+_species.speciesName+" = color rgb < ");
		fr.write((_species.color.getRed()) / 255.0 + " , ");
		fr.write((_species.color.getGreen()) / 255.0 + " , ");
		fr.write((_species.color.getBlue()) / 255.0 + " >");
		fr.write(";\n");
	}

	/**
	 * \brief Obtain another instance of the same species (totally
	 * independent).
	 * 
	 * Implemented by classes that extend this class.
	 */
	public abstract Agent sendNewAgent() throws CloneNotSupportedException;

	public void createNewAgent()
	{
		try
		{
			Agent baby = (Agent) sendNewAgent();
			// Register the baby in the pathway guilds.
			baby.registerBirth();
		}
		catch (CloneNotSupportedException e)
		{
			System.out.println("At ActiveAgent: createNewAgent error " + e);
		}
	}

//	/**
//	 * \brief Notifies the simulation that this agent has become too small and
//	 * is then counted as dead.
//	 * 
//	 * Decreases the population of this species.
//	 * 
//	 * @param isStarving	Boolean noting whether the agent currently has
//	 * access to any resources.
//	 */
//	@Override
//	public void die(Boolean isStarving) {
//		super.die(isStarving);
//		// Unregister from the metabolic guilds.
//		unregisterFromAllActiveReactions();
//	}

	/**
	 * \brief Called at each time step of the simulation to compute cell
	 * growth, update size, and monitor cell death and division.
	 * 
	 * Also determines whether the agent has reached the size at which it must
	 * divide.
	 */
	protected void internalStep() {
	//		grow();
	//		updateSize();
		}

	/**
	 * \brief Put growth rates into effect by changing the particle masses.
	 * 
	 * We adjust the particle masses after calculating all the deltaParticle
	 * values so that the reactions occur simultaneously.
	 */
	public void grow() {
		updateGrowthRates();
		for (int i = 0; i < particleMass.length; i++)
			particleMass[i] += deltaParticle[i];
	}

	/**
	 * \brief Perform agent growth by calling all active reaction pathways.
	 * 
	 */
	protected void updateGrowthRates() {
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
				_netVolumeRate += deltaMass/getActiveParam().particleDensity[i];
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
	public void updateMass() {
		_totalMass = ExtraMath.sum(particleMass);
	}

	/**
	 * \brief Add the reaction to the list of known reactions
	 * 
	 * Add the reaction to the list of known reactions
	 * 
	 * @param aReaction	The reaction to add to the list
	 * @param useDefaultParam	Whether to use default reaction parameters or bespoke parameters have been specified in the protocol file
	 */
	public void addReaction(Reaction aReaction, Boolean useDefaultParam) {
		// Add the reaction to the list of known reaction
		reactionKnown.add(aReaction.reactionIndex);
	
		// Test if specific parameters exist for this reaction
		int index = aReaction.reactionIndex;
		boolean test = getActiveParam().soluteYield[index]==null;
	
		if (useDefaultParam||test) {
			// Use parameters defined in the reaction object
			reactionKinetic[index] = aReaction.getKinetic();
			soluteYield[index] = aReaction.getSoluteYield();
			particleYield[index] = aReaction.getParticulateYield();
		} else {
	
			// Use parameters defined in the speciesParam structure
			reactionKinetic[index] = getActiveParam().reactionKinetic[index];
			soluteYield[index] = getActiveParam().soluteYield[index];
			particleYield[index] = getActiveParam().particleYield[index];
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
	public void addActiveReaction(Reaction aReaction, Boolean useDefaultParam) {
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
	public Double getParticleMass(int particleIndex) {
		return particleMass[particleIndex];
	}

	/**
	 * \brief Return net growth rate of this agent
	 *
	 * Return net growth rate of this agent
	 *
	 * @return Double stating the net growth rate of this agent
	 */
	public Double getNetGrowth() {
		return _netGrowthRate;
	}

	/**
	 * \brief Return volume growth rate of this agent
	 *
	 * Return volume growth rate of this agent
	 *
	 * @return Double stating the volume growth rate of this agent
	 */
	public Double getVolGrowth() {
		return _netVolumeRate;
	}

	/**
	 * \brief Set net growth rate of this agent
	 *
	 * Set net growth rate of this agent
	 *
	 * @param value	Double stating the net growth rate of this agent
	 */
	public void setNetGrowth(Double value) {
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
	public Double[] getSoluteYield(int indexReaction) {
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
	public Double[] getReactionKinetic(int indexReaction) {
		return reactionKinetic[indexReaction];
	}

	/**
	 * \brief Return the parameters associated with this agent (speciesParam)
	 * 
	 * Return the parameters associated with this agent (speciesParam)
	 * 
	 * @return ActiveParam object containing the species associated with this agent
	 */
	public ActiveParam getActiveParam() {
		return (ActiveParam) _speciesParam;
	}
}
