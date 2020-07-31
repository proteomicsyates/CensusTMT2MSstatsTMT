package edu.scripps.yates.censustmt2msstatstmt;

public enum SPC_FILTER_TYPE {
	SEQUENCE, ION, PROTEIN_SITE;

	public static String getValuesString(String separator) {
		final StringBuilder sb = new StringBuilder();
		for (final SPC_FILTER_TYPE spc : values()) {
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			sb.append(spc.name());
		}
		return sb.toString();
	}

	public static SPC_FILTER_TYPE getValueOf(String name) {
		for (final SPC_FILTER_TYPE ret : values()) {
			if (ret.name().equalsIgnoreCase(name)) {
				return ret;
			}
		}
		return null;
	}
}
