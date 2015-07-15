enable_if_configured("bundled:cdd") or return;
enable_if_configured("bundled:jreality") or return;
prefer_now 'polytope::cdd';

jreality(load("1")->BB_VISUAL, File=>diff_with("1j"))
	and
jreality(load("2")->BB_VISUAL, File=>diff_with("2j"))
	and
jreality(load("3")->BB_VISUAL, File=>diff_with("3j"))
	and
jreality(load("4")->BB_VISUAL, File=>diff_with("4j"));


