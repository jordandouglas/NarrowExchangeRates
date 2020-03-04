package branchratemodel;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

public class ClockModelAveraging extends BranchRateModel.Base {

	
	final public Input<IntegerParameter> modelIndicator = new Input<>("indicator",  "Clock model indicator variable", Input.Validate.REQUIRED);
	final public Input<List<BranchRateModel>> modelsInput = new Input<>("model", "Clock models", new ArrayList<>());
	
	private IntegerParameter indicators;
	private List<BranchRateModel> models;

	@Override
	public void initAndValidate() {
		
		this.indicators = modelIndicator.get();
		this.models = modelsInput.get();
		
		
		this.indicators.setDimension(1);
		this.indicators.setBounds(0, this.models.size()-1);
		//this.indicators.setValue(Randomizer.nextInt(this.models.size()));
	}
	
	
	@Override
	public double getRateForBranch(Node node) {
		
		BranchRateModel clockModel = this.models.get(this.indicators.getValue());
		return clockModel.getRateForBranch(node);
	}


	
	
}
