package simulator.agent;

import utils.Vect;

public class Point {
    static int UNIQUE_ID = 0;
    int uid = ++UNIQUE_ID;
	Double[] position;
	Double[] velocity;
	Double[] force;
	
	public Point(int nDim) 
	{
		setPosition(Vect.zeros(nDim));
		velocity = Vect.zeros(nDim);
		force = Vect.zeros(nDim);
	}
	
	public Point(int nDim, double domain) 
	{
		setPosition(Vect.randomDirection(nDim,domain));
		velocity = Vect.zeros(nDim);
		force = Vect.zeros(nDim);
	}
	
	public Point(Double[] p) 
	{
		this.setPosition(p);
		this.velocity = Vect.zeros(p.length);
		this.force = Vect.zeros(p.length);
	}

	public int identifier() 
	{
        return uid;
    }
	
	/**
	 * \brief performs one euler step for the mechanical relaxation.
	 * The velocity is expressed as v = (sum forces) / (3 Pi diameter viscosity)
	 * Currently the viscosity of water is assumed.
	 * @param vSquare
	 * 			Highest squared velocity in the system
	 * @param dt
	 * 			Current timestep of the mechanical relaxation
	 * @param radius
	 * 			Radius of the Point
	 * @return vSquare, if the squared velocity of this point is higher vSquare
	 * is updated.
	 */
	public Double euStep(Double vSquare, double dt, Double radius) 
	{
		position = Vect.sum(position, Vect.product(velocity,dt));
		velocity = Vect.product(force, 1.0/(radius*0.01885));
		if (Vect.normSquare(velocity) > vSquare)
			vSquare = Vect.normSquare(velocity);
		Vect.reset(force);
		return vSquare;
	}
	
	float[] coord(Double radius) 
	{
		float[] coord = new float[position.length];
		for (int i = 0; i < position.length; i++) 
			coord[i] = (float) (position[i]-radius);
		return coord;
	}
	
	float[] dimensions(Double radius) 
	{
		float[] dimensions = new float[position.length];
		for (int i = 0; i < position.length; i++) 
			dimensions[i] = (float) (radius*2.0);
		return dimensions;
	}
	
	float[] upper(Double radius) 
	{
		float[] coord = new float[position.length];
		for (int i = 0; i < position.length; i++) 
			coord[i] = (float) (position[i]+radius);
		return coord;
	}
	
	public int nDim() 
	{
		return position.length;
	}

	public Double[] getPosition() 
	{
		return this.position;
	}

	void setPosition(Double[] position) 
	{
		this.position = position;
	}

}
