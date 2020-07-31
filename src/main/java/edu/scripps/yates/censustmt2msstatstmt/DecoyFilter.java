package edu.scripps.yates.censustmt2msstatstmt;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;

public class DecoyFilter {
	private final Pattern decoyPattern;

	public DecoyFilter(String decoyPaternString) {
		decoyPattern = Pattern.compile(decoyPaternString);
	}

	public boolean isDecoy(String acc) {

		final Matcher matcher = decoyPattern.matcher(acc);
		if (matcher.find()) {
			return true;
		}
		return false;
	}

	public boolean isDecoy(QuantifiedPSMInterface psm) {
		final List<String> accs = psm.getQuantifiedProteins().stream().map(p -> p.getAccession())
				.collect(Collectors.toList());
		for (final String acc : accs) {
			if (!isDecoy(acc)) {
				return false;
			}
		}
		return true;
	}
}
