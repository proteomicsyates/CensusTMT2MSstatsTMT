package edu.scripps.yates.censustmt2msstatstmt;

public enum PSMAggregationType {
	SUM, AVG, HIGHEST;

	public static String printValidValues() {
		final StringBuilder sb = new StringBuilder();
		for (final PSMAggregationType p : values()) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(p.name());
		}
		return sb.toString();
	}
}
