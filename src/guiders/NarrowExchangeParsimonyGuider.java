package guiders;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.parsimony.FitchParsimony;

public class NarrowExchangeParsimonyGuider extends BEASTObject implements NarrowExchangeGuider {

	public final Input<Alignment> dataInput = new Input<>("data", "Multiple sequence alignment for computing parsimony scores", Input.Validate.REQUIRED);
	
	
	
	double[] scores = new double[] { };
	int currentNeighbourNum = 0;
	
	FitchParsimony fitch;
	Alignment alignment;
	Tree treeClone;
	
	@Override
	public void initAndValidate() {
		alignment = dataInput.get();
		fitch = new FitchParsimony(alignment, false);
	}

	@Override
	public double[] getProbabilityVector() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void reset(int numNeighbours) {
		if (scores.length != numNeighbours) scores = new double[numNeighbours];
		this.currentNeighbourNum = 0;
		
		for (int i = 0; i < numNeighbours; i ++) scores[i] = 0.0;
		
	}

	@Override
	public void addNeighbour(Tree tree) {
		fitch.reset();
		scores[currentNeighbourNum] = fitch.getScore(tree);
		currentNeighbourNum++;
	}
	
	



}
