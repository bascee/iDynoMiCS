package simulator.diffusionSolver.multigrid;

import simulator.geometry.IsComputationDomain;
import simulator.SoluteGrid;
import utils.ExtraMath;

/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *_____________________________________________________
 * Implements static utility functions for used in multigrid method.
 * 
 */
 
/**
 * @since June 2006
 * @version 1.0
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * 
 */


public class SinglegridPressure
{
	/**
	 * 
	 */
	public String soluteName;
	
	/**
	 * 
	 */
	public SoluteGrid realGrid;
	
	/**
	 * 
	 */
	protected Double _referenceSystemSide;
	
	/**
	 * 
	 */
	protected Double _diffusivity;
	
	/**
	 * 
	 */
	protected IsComputationDomain _domain;

	/**
	 * 
	 */
	protected Double sBulkMax, sBulk;
	
	/**
	 * 
	 */
	protected SoluteGrid _relDiff, _bLayer;
	
	/**
	 * 
	 */
	public SoluteGrid _conc, _reac;
	
	/**
	 * 
	 */
	protected SoluteGrid _rhs, _itemp, _itau;
	
	/**
	 * 
	 */
	public Double truncationError;
	
	/**
	 * Not used (?)
	 */
	//private static final Double[][][] _diff    = new double[3][3][3];
	
	/**
	 * 
	 */
	private static Double[][][] u;
	
	/**
	 * Not used (?)
	 */
	//private static Double[][][] rd;
	
	/**
	 * 
	 */
	private static Double[][][] bl;
	
	/**
	 * 
	 */
	private static int _i, _j, _k;
	
	/**
	 * 
	 */
	public static final double BLTHRESH = 0.1;
	
	/**
	 * 
	 */
	@SuppressWarnings("unused")
	private static int maxOrder;
	
	/**
	 * 
	 */
	private static int _nI, _nJ, _nK;

	/* ____________ ______________________ */
	public SinglegridPressure(SoluteGrid aSolute, SoluteGrid bLayer, Double sBulk)
	{
		realGrid = aSolute;
		realGrid.diffusivity = 1.0;
		soluteName = realGrid.gridName;
		
		_nI = realGrid.getGridSizeI();
		_nJ = realGrid.getGridSizeJ();
		_nK = realGrid.getGridSizeK();
		_domain = realGrid.getDomain();
		
		setReferenceSide();
		
		this.sBulkMax = sBulk;
		this.sBulk = sBulk;
		
		_bLayer = bLayer;
		
		//for (int iGrid = 0; iGrid<maxOrder; iGrid++) {
		int iGrid=0;
		_i = (_nI-1)/ExtraMath.exp2(iGrid)+1;
		_j = (_nJ-1)/ExtraMath.exp2(iGrid)+1;
		_k = (_nK-1)/ExtraMath.exp2(iGrid)+1;
		double r = _referenceSystemSide/referenceIndex(_nI,_nJ,_nK);
		// Padding is automatically generated by the constructor
		_conc = new SoluteGrid(_i, _j, _k, r, aSolute);
		_rhs = new SoluteGrid(_i, _j, _k, r, aSolute);
		_reac = new SoluteGrid(_i, _j, _k, r, aSolute);
	}

	/* _______________ ______________________________________ */

	public double relax()
	{
		int nI = _conc.getGridSizeI();
		int nJ = _conc.getGridSizeJ();
		int nK = _conc.getGridSizeK();

		double h = _referenceSystemSide/referenceIndex(_nI,_nJ,_nK);
		double h2i = 0.5f/(h*h);
		// red-black relaxation
		// iterate through system
		// isw, jsw and ksw alternate between values 1 and 2
		u = _conc.grid;
		bl = _bLayer.grid;
		double lop, dlop, res;
		double totalError = 0;
		// bvm 22.12.09: now allows red-black for 2d AND 3d
		int ksw = 1;
		int isw, jsw;
		for (int pass = 1; pass<=2; pass++, ksw = 3-ksw)
		{
			jsw = ksw;
			for (_k = 1; _k<=nK; _k++, jsw = 3-jsw)
			{
				isw = jsw;
				for (_j = 1; _j <= nJ; _j++, isw = 3-isw)
				{
					for (_i = isw; _i <= nI; _i += 2)
					{
						if (bl[_i][_j][_k]>=BLTHRESH)
						{
							// Case: Inside boundary layer
							// Equations must be solved here
							// compute diffusivity values
							// and that of surrounding neighbors
							//fillDiff();
							// compute L operator
							lop = computeLop(h2i);
							// compute derivative of L operator
							dlop = computeDiffLop(h2i);
							// compute residual
							res = (lop-_rhs.grid[_i][_j][_k])/dlop;
							totalError += lop;
							// update concentration (test for NaN)
							//LogFile.writeLog("NaN generated in multigrid solver "+"while computing rate for "+soluteName);
							u[_i][_j][_k] -= res;
						}
					}
				}
			}
			// refresh the padding elements to enforce
			// boundary conditions for all solutes
			_conc.refreshBoundary();
		}
		return totalError;
	}

	//@SuppressWarnings("unused")
	//private void fillDiff() {
	//	_diff[0][1][1] = realGrid.diffusivity*rd[_i-1][_j][_k];
	//	_diff[2][1][1] = realGrid.diffusivity*rd[_i+1][_j][_k];
	//	_diff[1][0][1] = realGrid.diffusivity*rd[_i][_j-1][_k];
	//	_diff[1][2][1] = realGrid.diffusivity*rd[_i][_j+1][_k];
	//	_diff[1][1][0] = realGrid.diffusivity*rd[_i][_j][_k-1];
	//	_diff[1][1][2] = realGrid.diffusivity*rd[_i][_j][_k+1];
	//	_diff[1][1][1] = realGrid.diffusivity*rd[_i][_j][_k];
	//}

	private double computeLop(double h2i)
	{
		return ((2)*(u[_i+1][_j][_k]-u[_i][_j][_k])+(2)*(u[_i-1][_j][_k]-u[_i][_j][_k])+(2)
		        *(u[_i][_j+1][_k]-u[_i][_j][_k])+(2)*(u[_i][_j-1][_k]-u[_i][_j][_k])+(2)
		        *(u[_i][_j][_k+1]-u[_i][_j][_k])+(2)*(u[_i][_j][_k-1]-u[_i][_j][_k]))
		        *h2i+_reac.grid[_i][_j][_k];
	}
	
	/** 
	 * TODO do properly (it's kinda just coppied from MultigridSolute)
	 * 
	 * @param h2i
	 * @return
	 */
	private double computeDiffLop(double h2i)
	{
		return -h2i*(12)+0;
	}
	
	/* _________________________ TOOLBOX ____________________________ */
	
	/**
	 * 
	 */
	public void resetMultigridCopies()
	{
		setSoluteGridToBulk();
		_reac.resetToZero();;
		_rhs.resetToZero();;
	}
	
	/**
	 * Set all grids elements to the value defined for Bulk. For elements
	 * located in the convective part (i.e. outside the BLayer, we take the
	 * value defined in the BulkBoundary Class)
	 */
	public void setSoluteGridToBulk()
	{
		for (_i = 1; _i <= _conc.getGridSizeI(); _i++)
			for (_j = 1; _j <= _conc.getGridSizeJ(); _j++)
				for (_k = 1; _k <= _conc.getGridSizeK(); _k++)
					if (_bLayer.grid[_i][_j][_k] <= BLTHRESH)
					{
						// outside the boundary layer (will not be solved)
						_conc.grid[_i][_j][_k] = sBulk;
					}
					else
					{
						// inside the biofilm (value is not really important
						// now)
						_conc.grid[_i][_j][_k] = sBulkMax;
					}
	}
	
	/**
	 * Determine order of the finest grid.
	 * 
	 */
	public void setReferenceSide()
	{
		_referenceSystemSide = (double) Math.min(_nI, _nJ);
		if (_nK > 1)
			_referenceSystemSide = Math.min(_referenceSystemSide, _nK);
		maxOrder = ExtraMath.log2(_referenceSystemSide).intValue();
		_referenceSystemSide -= 1;
		_referenceSystemSide *= realGrid.getResolution();
	}
	
	/**
	 * \brief This is meant to return the correct index value following
	 * the logic of setReferenceSide() above.
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return
	 */
	private double referenceIndex(int i, int j, int k)
	{
		if (_nK > 1)
			return Math.min(i, Math.min(j, k)) - 1;
		return Math.min(i, j) - 1;
	}
}
