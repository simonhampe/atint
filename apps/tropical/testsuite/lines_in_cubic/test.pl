compare_output {
	my $L = lines_in_cubic(load_data("poly"));
	print $L->N_ISOLATED,"\n";
	print $L->N_FAMILIES,"\n";
} 'out'
