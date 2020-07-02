#! /usr/bin/perl -w

if (@ARGV < 2) {
	die "usage: ntfreq.pl <filename> x\nwhere x = 1 (short output) or 2 (long output)\n";
} elsif (($ARGV[1] != 1) && ($ARGV[1] != 2) ) {
	die "usage: ntfreq.pl <filename> x\nwhere x = 1 (short output) or 2 (long output)\n";
}


%seqs = readFasta($ARGV[0]);

# for each sequence, calculate
#       n = number of nucleotides
#       nx = freq of nucleotide x
#       nxy = freq of dinucleotide xy


foreach $id (keys %seqs) {
		
        $seq = $seqs{$id};
        $n = length($seq);

        # two hash tables are used to collect counts for nucleotide/dinucleotide frequencies
        # Each of the nucleotides and dinucleotides have to be initialized to zero
	# in order to have aligned columns in the resulting output
 
	%nx = ();
	$nx{"a"}++; $nx{"a"}--; 
	$nx{"c"}++; $nx{"c"}--; 
	$nx{"g"}++; $nx{"g"}--; 
	$nx{"t"}++; $nx{"t"}--; 
	$nx{"n"}++; $nx{"n"}--; 

	%nxy = (); 
	$nxy{"aa"}++; $nxy{"aa"}--;
	$nxy{"ac"}++; $nxy{"ac"}--;
	$nxy{"ag"}++; $nxy{"ag"}--;
	$nxy{"at"}++; $nxy{"at"}--;
	$nxy{"an"}++; $nxy{"an"}--;
	$nxy{"ca"}++; $nxy{"ca"}--;
	$nxy{"cc"}++; $nxy{"cc"}--;
	$nxy{"cg"}++; $nxy{"cg"}--;
	$nxy{"ct"}++; $nxy{"ct"}--;
	$nxy{"cn"}++; $nxy{"cn"}--;
	$nxy{"ga"}++; $nxy{"ga"}--;
	$nxy{"gc"}++; $nxy{"gc"}--;
	$nxy{"gg"}++; $nxy{"gg"}--;
	$nxy{"gt"}++; $nxy{"gt"}--;
	$nxy{"gn"}++; $nxy{"gn"}--;
	$nxy{"ta"}++; $nxy{"ta"}--;
	$nxy{"tc"}++; $nxy{"tc"}--;
	$nxy{"tg"}++; $nxy{"tg"}--;
	$nxy{"tt"}++; $nxy{"tt"}--;
	$nxy{"tn"}++; $nxy{"tn"}--;
	$nxy{"na"}++; $nxy{"na"}--;
	$nxy{"nc"}++; $nxy{"nc"}--;
	$nxy{"ng"}++; $nxy{"ng"}--;
	$nxy{"nt"}++; $nxy{"nt"}--;
	$nxy{"nn"}++; $nxy{"nn"}--;
  

	# loop through each of the nucleotides in the sequence
        for ($i=0; $i<$n; $i++) {

		# if a given nucleotide is found then its number is incremented in the hash table
                $nx{lc(substr($seq, $i, 1))}++; 
 		
		# skip the very last base
                if ($i < $n-1) {     

			# if a given dinucleotide is found then its number is incremented in the hash table
                        $nxy{lc(substr($seq, $i, 2))}++;  
                }
        }
     
	# short output
	if ($ARGV[1] == 1){

		# print out counts 
       		print "$id:";
        	
		# calculate and print out the CpGOE
		# CpGOE = N(cg) / (N(c) * N(g)) 
		if ($nxy{"cg"} > 0) {
			$freq = ($nxy{"cg"}*$n)/(($nx{"c"})*($nx{"g"}));
        		print "$freq\n";
		} else {
			print "N/A\n";
		}


	# long output
	} elsif ($ARGV[1] == 2)  {

		# print out counts 
       		 print "$id:  ~N (# of nucleotides) ~ $n". "~";
        	
	
		# print out the number of each nucleotide in the sequence
        	foreach $nuc (sort keys %nx) {
                	print "$nuc" . "~$nx{$nuc}" . "~";
        	}
	
		# print out the number of each dinucleotide in the sequence  #added sort keys to se if it would help output
        	foreach $dinuc (sort keys %nxy) {
                	print "$dinuc" . "~$nxy{$dinuc}" . "~";
        	}
	
		# calculate and print out the CpGOE
		# CpGOE = N(cg) / (N(c) * N(g)) 
		if ($nxy{"cg"} > 0) {
			$freq = ($nxy{"cg"}*$n)/(($nx{"c"})*($nx{"g"}));
        		print "CpGOE~$freq\n";
		} else {
			print "CpGOE~N/A\n";
		}	
	}
}



# ------
# subroutine readFasta
#       input : filename
#       output : hashtable containing all the sequences

sub readFasta {
        my %sequences;    # make this a local variable
        my $filename = $_[0];

        open (FILE, $filename) or die "Cannot open $filename";

        while (<FILE>) {
                chomp;
                if (/^>s*(.+)/) {   # if it's the name
                        $name = $1;
                } elsif ($_ ne "") {
                        $sequences{$name} .= $_;
                }
        }
        return %sequences;
}
