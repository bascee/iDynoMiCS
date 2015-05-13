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

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;

import idyno.SimTimer;
import simulator.AgentContainer;
import simulator.Simulator;
import simulator.SpatialGrid;
import simulator.geometry.ContinuousVector;
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
	 * \brief Initialise an agent object, setting its time of creation and
	 * thus the time it was last stepped.
	 */
	public Agent()
	{
		_birthday = SimTimer.getCurrentTime();
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
		// read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT
		_family     = Integer.parseInt(singleAgentData[0]);
		_genealogy  = new BigInteger(singleAgentData[1]);
		_generation = Integer.parseInt(singleAgentData[2]);
		_birthday   = Double.parseDouble(singleAgentData[3]);
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
		return out;
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
		internalStep();
	}
	
	/**
	 * \brief Called at each time step (under the control of the method Step
	 * of the class Agent to avoid multiple calls).
	 * 
	 * Implemented by classes that extend Agent.
	 */
	protected abstract void internalStep();
	
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
		return new StringBuffer("family,genealogy,generation,birthday");
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
		return new StringBuffer(
						_family+","+_genealogy+","+_generation+","+_birthday);
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

	public void interact(boolean mUTUAL) {
		// TODO Auto-generated method stub

	}

	public void die(boolean b) {
		// TODO Auto-generated method stub
		
	}

	public void fitMassOnGrid(SpatialGrid biomassGrid) {
		// TODO Auto-generated method stub
		
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

	public void createNewAgent() {
		// TODO Auto-generated method stub
	}

	public Double getTotalMass() {
		// TODO Auto-generated method stub
		return null;
	}

	public Double getVolume(boolean b) {
		// TODO Auto-generated method stub
		return null;
	}

	public void addMovement(ContinuousVector move) {
		// TODO Auto-generated method stub
		
	}

	public ContinuousVector getLocation() {
		// TODO Auto-generated method stub
		return null;
	}

	public Double getShoveRadius() {
		// TODO Auto-generated method stub
		return null;
	}

	public void updateSize() {
		// TODO Auto-generated method stub
		
	}

	public boolean willDie() {
		// TODO Auto-generated method stub
		return false;
	}
	
	public void setDetPriority(Double priority) {
	}

	public void addDetPriority(Double detFunction) {
		// TODO Auto-generated method stub
	}

	public void multiplyDetPriority(Double multiplier) {
		// TODO Auto-generated method stub
		
	}

	public Double getRadius(boolean b) {
		// TODO Auto-generated method stub
		return null;
	}

	public void setDistProb(Double prob) {
		// TODO Auto-generated method stub
	}

	public Double getDistProb() {
		// TODO Auto-generated method stub
		return null;
	}

	public void setDistCumProb(Double prob) {
		// TODO Auto-generated method stub
		
	}

	public Double getDistCumProb() {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean willDivide() {
		// TODO Auto-generated method stub
		return false;
	}

	public void setParticleMass(Double[] particleMass) {
		// TODO Auto-generated method stub
		
	}

	public Double[] getParticleMass() {
		// TODO Auto-generated method stub
		return null;
	}

	public void addToParticleMasses(Double[] mass) {
		// TODO Auto-generated method stub
		
	}

	public void subtractFromParticleMasses(Double[] mass) {
		// TODO Auto-generated method stub
		
	}

	public void addParticleMass(Double mass, int particleIndex) {
		// TODO Auto-generated method stub
		
	}

	public void multiplyParticleMass(Double multiplier, int particleIndex) {
		// TODO Auto-generated method stub
		
	}

	public Double getDetPriority() {
		// TODO Auto-generated method stub
		return null;
	}

	protected void setTotalMass(Double _totalMass) {
		// TODO Auto-generated method stub
	}

	public double[] getBoundingBoxCoord() {
		// TODO Auto-generated method stub
		return null;
	}

	public double[] getSearchCoord(Double range) {
		// TODO Auto-generated method stub
		return null;
	}

	public double[] getBoundingBoxDimensions() {
		// TODO Auto-generated method stub
		return null;
	}

	public void setGridIndex(int aGridIndex) {
		// TODO Auto-generated method stub
		
	}

	public int getGridIndex() {
		// TODO Auto-generated method stub
		return 0;
	}

	public Double getNetGrowth() {
		// TODO Auto-generated method stub
		return null;
	}

	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex) {
		// TODO Auto-generated method stub
		
	}

	public void fitReacRateOnGrid(SpatialGrid aRateGrid, int reactionIndex) {
		// TODO Auto-generated method stub
		
	}

	public ContinuousVector getMovement() {
		// TODO Auto-generated method stub
		return null;
	}

	protected void updateGrowthRates() {
		// TODO Auto-generated method stub
		
	}

	public void createNewAgent(ContinuousVector position) {
		// TODO Auto-generated method stub
		
	}

	public Color getColor() {
		// TODO Auto-generated method stub
		return null;
	}

	public Double getActiveFrac() {
		// TODO Auto-generated method stub
		return null;
	}

	public Boolean hasEPS() {
		// TODO Auto-generated method stub
		return null;
	}

	public Color getColorCapsule() {
		// TODO Auto-generated method stub
		return null;
	}

	protected void setMyDivRadius(Double _myDivRadius) {
		// TODO Auto-generated method stub
		
	}

	protected void setMyDeathRadius(Double _myDeathRadius) {
		// TODO Auto-generated method stub
		
	}

	public void subtractMovement(ContinuousVector _divisionDirection) {
		// TODO Auto-generated method stub
		
	}

	protected void setNetVolumeRate(Double _netVolumeRate) {
		// TODO Auto-generated method stub
		
	}

	public void divide() {
		// TODO Auto-generated method stub
		
	}

	public void setLocation(ContinuousVector position) {
		// TODO Auto-generated method stub
		
	}

	public LocatedParam getLocatedParam() {
		// TODO Auto-generated method stub
		return null;
	}

	public ContinuousVector getVerifiedLocationFromMovement(ContinuousVector continuousVector) {
		// TODO Auto-generated method stub
		return null;
		
	}

	public void interact(Double dt, boolean MUTUAL) {
		// TODO Auto-generated method stub
		
	}

	public void environment() {
		// TODO Auto-generated method stub
		
	}

	protected ContinuousVector getVelocity() {
		// TODO Auto-generated method stub
		return null;
	}

	protected ContinuousVector getForce() {
		// TODO Auto-generated method stub
		return null;
	}

	protected void setForce(ContinuousVector _force) {
		// TODO Auto-generated method stub
		
	}

	public void updateMovement(Double dt) {
		// TODO Auto-generated method stub
		
	}
}
