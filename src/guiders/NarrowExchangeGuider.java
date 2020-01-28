package guiders;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public abstract interface NarrowExchangeGuider {

	

	
	
	/***
	 * @return A cumulative probability vector [p1,p2,p3,...] of sampling the added trees [1,2,3,...]
	 */
	public double[] getProbabilityVector();
	
	/***
	 * Clear the cache of tree neighbours
	 */
	public void reset(int numNeighbours);
	
	/***
	 * Add a new tree neighbour
	 * @param tree
	 */
	public void addNeighbour(Tree tree);
	
	
	

}
