package edu.scripps.yates.censustmt2msstatstmt.filters;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.censustmt2msstatstmt.util.DataUtil;

public class PurityFilter {
	private final float minimumPurity;

	public PurityFilter(float minimumPurity) {
		this.minimumPurity = minimumPurity;
	}

	public boolean passFilter(QuantifiedPSMInterface psm) {
		return DataUtil.isTMTPurityValid(psm, this.minimumPurity);
	}
}
