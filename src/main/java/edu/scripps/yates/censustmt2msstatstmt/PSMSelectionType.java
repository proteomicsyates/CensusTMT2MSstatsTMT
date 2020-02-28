package edu.scripps.yates.censustmt2msstatstmt;

public enum PSMSelectionType {
	SUM, AVERAGE, HIGHEST;

	public static String printValidValues() {
		final StringBuilder sb = new StringBuilder();
		for (final PSMSelectionType p : values()) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(p.name());
		}
		return sb.toString();
	}
}
