package edu.scripps.yates.censustmt2msstatstmt.singletons;

import java.util.ArrayList;

public class PTMList extends ArrayList<Float> implements Clearable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -7641065478839736484L;
	private static PTMList instance;
	public static final double MASS_PRECISION = 0.01;

	public static PTMList getInstance() {
		if (instance == null) {
			instance = new PTMList();
		}
		return instance;
	}

	private PTMList() {

	}

	@Override
	public void clearData() {
		this.clear();
	}
}
