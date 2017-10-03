#! /usr/bin/perl

#	do the asymmetric omega ratio test in the triplet pattern
#	create_codemlctl.pl should be ready and added to the ENV or its value should be OK (see below)
#	Author: Sishuo Wang from the University of British Columbia (sishuowang@hotmail.ca)
#	the name of the target sequences should be used by specifying "-t|-target 123"


########################################################################
use strict;
use 5.010;
use Getopt::Long;


########################################################################
my ($input, @target, $quiet_swi, @non_target, $help, $mlc_file, $is_show_rates);
my ($seqtype, $model, $codon);
my ($treefile);
my (@lnL, %rates);
our $counter = 0;

my $create_codemlctl = "/home/sswang/tools/self_bao_cun/PAML/create_codemlctl.pl";


#########################################################################
GetOptions(
	'i|input=s'     =>	\$input,	# paml format
	't|target=s'    =>	\@target,
	'q|quiet!'      =>	\$quiet_swi,
	#'non_target|non_t=s'	=>	\@non_target,
    'seqtype=s'     =>  \$seqtype,
    'model=s'       =>  \$model,
    'codon|codon_table=s'       =>  \$codon,
    'show_rates|show_rate!' =>  \$is_show_rates,
	'help|h!'		=>	\$help,
) || die "illegal param!\n";

&show_help() if $help;
$treefile = $input . ".tree";
die "target number should be 2!\n" if @target != 2; 


#########################################################################
my $seq_href = &read_input_seq($input);

my ($target_title_href, $non_target_title_href) = &generate_target_non_target([keys %$seq_href], [@target]);

$counter++;
&triplet_test([keys %$target_title_href, keys %$non_target_title_href]);
$counter++;
&triplet_test([keys %$target_title_href], [keys %$non_target_title_href]);

&chi2_calculate(\@lnL);

&show_rates(\%rates, $target_title_href) if $is_show_rates;


#########################################################################
sub show_rates{
    my ($rates_href, $target_title_href) = @_;
    foreach my $key (keys %$rates_href){
        next if not exists $target_title_href->{$key};
        my $omega = $rates_href->{$key};
        print $key."\t".$omega."\n";
    }
}


sub chi2_calculate{
	my ($lnL_aref) = @_;
	my @lnL = @$lnL_aref;
	my $abs_diff = 2 * abs($lnL[0] - $lnL[1]);
	#print "abs:\t$abs_diff\n";
	#system "chi2 1 $abs_diff";
	open (IN, "chi2 1 $abs_diff |");
	my @in=<IN>;
	@in = grep {$_=~/\S/} @in;
	my ($prob) = $in[-1] =~ /(\S+)\s*$/;
	print $prob."\n";
}


sub triplet_test{
    my $clock=0;
	my ($target_title_aref, $non_target_title_aref) = @_;
	given (@$non_target_title_aref){
		when ($_){
            if (($seqtype == 2 or $seqtype == 3) and $counter == 1){
                $clock = 1;
            }
			&create_triplet_tree_for_paml([@$target_title_aref], [@$non_target_title_aref]);
		}
		when (! $_){
			&create_triplet_tree_for_paml([@$target_title_aref]);
		}
	}
	$mlc_file = 'mlc' if not defined $mlc_file;
	&run_codeml($seqtype, $model, $clock, $codon);
	&get_lnL($mlc_file, $target_title_aref);
}


sub generate_target_non_target{
	my (%seq_title, %target_title, %non_target_title);	# seq titles are desposited here
	my ($seq_name_aref, $target_aref) = @_;
	map {$seq_title{$_}=1} @$seq_name_aref;
	map {$target_title{$_} = 1} @$target_aref;
	map {$non_target_title{$_} = 1} grep {! exists $target_title{$_}} keys %seq_title;
	foreach(keys %target_title){
		print "target $_ is not found in sequence file!\n" and die if not exists $seq_title{$_};
	}
	return (\%target_title, \%non_target_title);
}


sub create_triplet_tree_for_paml{
	my (%seq_title, %target_title, %non_target_title);	# seq titles are desposited here
	my ($target_title_aref, $non_target_title_aref) = @_;
	#print "@$target_title_aref\n";
	open (my $OUT, '>', "$input".'.tree') || die;
	select $OUT;
	print "(";
	map {print "$_,"} @$non_target_title_aref;
    print $target_title_aref->[2] . ',' if scalar @$target_title_aref == 3;
    print "(";
    given (scalar @$target_title_aref){
        when (2)	{print join (',', map {my $a=$_+1; $target_title_aref->[$_] . "#" . $a} 0..1)}
        when (3)	{print join (',', map {$target_title_aref->[$_] . "#1"} 0..1)}
    }
    print ")";
	print ")";
	print ";";
	print "\n";
	close $OUT; select STDOUT;
}


sub get_lnL{
	my ($mlc_lnL, $target_title_aref) = @_;
	open (my $IN, '<', "$mlc_lnL") || die "mlc_lnL $mlc_lnL cannot be opened!\n";
	while(<$IN>){
		next if /^$/;
		chomp;
		if (/^lnL [^-]+ (\-\d+(?:\.\d+)?)/x){
			push @lnL, $1;
		}
        if ($seqtype == 1){
            if (/w ratios as labels for TreeView:/){
                my $line = <$IN>;
                chomp($line);
                #(AT1G06900 #0.1008 , (Bra015536_LF #0.2516 , Bra030665_MF1 #0.1888 ) #0.1008 );
                foreach my $title (@$target_title_aref){
                    if ($line =~ /($title) [#](\d+\.\d+)/){
                        $rates{$1}=$2;
                    }
                }
                last;
            }
        }
        elsif ($seqtype == 2 or $seqtype == 3){
            if (/^tree length/){
                my $line = <$IN>;
                my $line = <$IN>;
                my $line = <$IN>;
                my $line = <$IN>;
                chomp($line);
                #(AT1G06900: 0.05954, (Bra015536_LF: 0.05929, Bra030665_MF1: 0.05929): 0.00025);
                foreach my $title (@$target_title_aref){
                    if ($line =~ /($title): (\d+\.\d+)/){
                        $rates{$1}=$2;
                    }
                }
                last;
            }

        }
	}
}


sub read_input_seq
{
    my ($input_seq)=@_;
    open(IN,$input_seq);
    my $line1=<IN>;
    chomp ($line1);
    my %seq;
    my ($seq_num,$num_char) = $line1 =~ m/(\d+)\s+(\d+)/; # i.e.  2  531
    my ($count_line,$add,$num_line_per_seq) = (0,0,0);

    #print $seq_num."\t".$num_char."\n";
    $num_char%60 ? ($add=1) : ($add=0);
    $num_line_per_seq=int($num_char/60)+1+$add;
    #print $num_line_per_seq."\n";

    my $seq_name;
    while(my $line=<IN>){
        chomp($line);
        ++$count_line;
        ($count_line % $num_line_per_seq == 1) ? ($seq_name=$line) : ($seq{$seq_name}.=$line);
    }
    close IN;
    return (\%seq);
}


sub run_codeml{
	no strict;
    my ($seqtype, $model, $clock, $codon) = @_;
    $model = 2 if not defined $model;
	print "codeml starts to be run ......\n" if not $quiet_swi;
	my ($output_codemlctl_dir, $NSsites);
	$output_codemlctl_dir="./";
	$NSsites = 0;
	`\"$create_codemlctl\" -seqfile=\"$input\" -treefile=\"$treefile\" -outfile=\"$outfile\" -model=\"$model\" -NSsites=\"$NSsites\" -ncatG=\"$ncatG\" -seqtype=\"$seqtype\" -clock=\"$clock\" -output_codemlctl_dir=\"$output_codemlctl_dir\" -codon=\"$codon\"`;
	`codeml`;
}


#####################################################################
sub show_help{
	print "\tUSAGE:\tperl $0 <-input infile> <-t 1> <-t 2>\n\n";
	exit;
}


