#!/usr/bin/perl -w
use strict;
use File::Path;
use List::Util qw[min max];
use CGI;
use CGI::Carp qw(fatalsToBrowser);
$CGI::POST_MAX = 1024 * 1024;

my $cgi = CGI->new();



# print header
print "Content-type: text/html\n\n";

print <<"HTML_header";

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>ZEOGS: Zebrafish Expression Ontology</title>
<link rel="stylesheet" type="text/css" href="../view.css" media="all">
<script type="text/javascript" src="../view.js"></script>

</head>
<body id="main_body" >
	
	<img id="top" src="../ZEOGS_banner.png" alt="">
		
	<div id="form_container">
	
	<table summary="Header table for links to different website pages" class="links">
	<thead class="links"><tr>
	<th scope="col" class="links"><a href="http://zeogs.x10.bz/index.html">ZEOGS form</a></th>
	<th scope="col" class="links"><a href="http://zeogs.x10.bz/algorithm.html">Algorithm scheme</a></th>
	<th scope="col" class="links"><a href="http://zeogs.x10.bz/help.html">Help</a></th>
	<th scope="col" class="links"><a href="http://zeogs.x10.bz/tutorial.html">Tutorial</a></th>
	<th scope="col" class="links"><a href="http://zeogs.x10.bz/links.html">Links</a></th>
	<th scope="col" class="links"><a href="http://zeogs.x10.bz/contact.html">Contact</a></th>
	</tr>
	</thead>
	</table>
	
HTML_header



###################################
# STEP 1 - Obtaining the input data
###################################

# Values for option parameters need to be obtained
my %params = $cgi->Vars();
my $input = $params{'gene_list'}; 
my $format  = $params{'format'}; # possible formats are: ZFIN ID, ZFIN gene symbol
my $stage   = $params{'stage'};  # stages are numerical values from 1 to 44
my $path = '../uploads';
my $filehandle = $cgi->upload('upload_file');
my $filename = $params{'upload_file'};
my $data_source = $params{'data_source'};
my $sorting_option = $params{'sorting'};

# The first task of the program is to parse the input from the form
# OR
# Get the data from the indicated file

# First, gene ID data need to be obtained either from the text area or a file
# array for storing these IDs is initialized, a specific order of elements is not important
my @geneids_input = ();

# treat this input as something potentially unreliable and in need of additional examination and processing
if ($input ne ''){	# the first alternative is to enter all the data into the text area
	@geneids_input  = split(/\n/, $input); 
	
	# Do additional examination and processing of gene identifiers
	# 1) remove any superfluous characters from each of the gene identifiers
	
	foreach (@geneids_input){
		chomp($_);
		$_ = trim($_);
	}
} else { # the second alternative is to upload the data from the file 
	
	$filename = (split(/[\\\/]/, $filename))[-1];
	$filename =~ s/[^A-Za-z0-9_\.\-]//g;
	
	# the upload procedure should output the contents of the file in an array
	@geneids_input = upload_file($filename,$filehandle,$path);
	
	foreach (@geneids_input){
		chomp($_);
		$_ = trim($_);
	}
	
	}
	

	# 2) verify that each of the entered gene identifiers satisfies the requirements for certain gene identifier formats
	# Here it only makes sense to verify Ensembl, Entrez and ZFIN gene IDs since for them one can define suitable 
	# regular expressions. For ZFIN gene symbols one could only define very fuzzy regular expressions, so this program leaves
	# it up to the user to ensure the correctness of this type of gene IDs
	 

	my @geneids = ();	# an array to store verified gene IDS
	
# do the verification of ZFIN gene identifier formas using regular expressions


		if ($format eq 'zfin_id') {
			foreach (@geneids_input) {
				if ($_ =~ /ZDB-GENEP?-[0-9]{6}-[0-9]{1,4}/) {push(@geneids, $_);}
			}
		} else {
			@geneids = @geneids_input; 
		}


	
	
############################################################
# STEP 2 - Converting the input gene IDs to ZFIN IDs
############################################################

# In this step, the gene IDs from the previous step are taken and converted to Ensembl gene IDs according
# to the format input information

# initialize an array for storing ZFIN gene IDs resulting from conversion of current gene IDs
my @ZFIN_IDs = ();
my %conversion_table = (); # a hash for storing the correspondence between the input gene IDs and ZFIN gene IDs 

# if loops are used to determine which conversion table should be read in or none at all
if ($format eq 'zfin_id') {
	
	# no conversion is necessary
	@ZFIN_IDs = @geneids;

} elsif ($format eq 'zfin_symbol') {
	
	open ZFINSYMB, "<../zeogs_data/ZFIN_symbol-ZFIN_IDs_map.txt" or die $!;
		
		# initialize the variables necessary for holding gene IDs
		my $zfin_symb;
		my $zfin_id;
		
		# read the data from the indicated file and populate the conversion table hash 
		while(<ZFINSYMB>) {
			
			# obtain values for 
			($zfin_symb, $zfin_id) = split(/\t/, $_);
			# clean up these variables
			chomp($zfin_symb); chomp($zfin_id);
			
			if(exists $conversion_table{$zfin_symb}){
				# add a new Ensembl gene ID to the array reference stored under the key of $ensembl
				push(@{$conversion_table{$zfin_symb}}, $zfin_id);
			
			} elsif ($zfin_id =~ /ZDB-GENEP?-[0-9]{6}-[0-9]{1,4}/) # zfinID should match its correct pattern
			{
				# create a new array reference under the key of $ensembl
				$conversion_table{$zfin_symb} = [$zfin_id];
			}	
		}
		close ZFINSYMB; 
	
		# use values of gene IDs in the array to obtain the corresponding ZFIN IDs from the conversion table hash
		# and store them in the ZFIN_IDs hash
		
		# a hash for ZFIN_IDs needs to be defined to make sure that in the end all of the converted gene IDs are unique
		
		my %ZFIN_IDs = (); 
		
		foreach my $id (@geneids) {
			if (exists $conversion_table{$id}){
				# run through the array corresponding to this ID and collect unique ZFIN gene IDs
				foreach (@{$conversion_table{$id}}) {
					unless(exists $ZFIN_IDs{$_}) {
						$ZFIN_IDs{$_} = 1;
					}
				}
			}
		}
		# fill the ZFIN gene ID array with ZFIN gene IDs stored in the corresponding hash
		@ZFIN_IDs = keys %ZFIN_IDs;
	} 

#######################################################
# STEP 2.2 - Reading the total gene list selection file 
#######################################################
 
# an array for storing uploaded gene IDs is defined
my $file_total_genes = $params{'file_genes_tested'};
my $fh_total_genes = $cgi->upload('file_genes_tested');
my @total_genes_tested = ();

# treat this input as something potentially unreliable and in need of additional examination and processing
if ($file_total_genes ne ''){	# check that the user entered some kind of filename
		
	$file_total_genes = (split(/[\\\/]/, $file_total_genes))[-1];
	$file_total_genes =~ s/[^A-Za-z0-9_\.\-]//g;
	
	# the upload procedure should output the contents of the file in an array
	@total_genes_tested = upload_file($file_total_genes, $fh_total_genes, $path); 
	
	foreach (@total_genes_tested){
		chomp($_);
		$_ = trim($_);
	}
}
	
#######################################################################
# STEP 3 - Obtaining the expression term information from a stored file
#######################################################################

##################################################################################
# STEP 4 - Filtering expression terms based on the stage-specific expression terms 
##################################################################################

# Algorithm
#
# 1. Read in expression terms from the main file into a hash with ZFIN IDs as keys and all the expression terms
#    into an array, references of which are values of the keys  
# 2. Get a slice of this hash for the available ZFIN IDs.
# 3. Put ZFIN IDs without expression term annotation into a special array.



# define a hash for storing expression terms
# and variables for temporary storing and processing of the information from the main data set
my %all_genes_patterns = (); # ZFIN-ID-centric hash
my %patterns_all_genes = (); # expression pattern-centric hash
my $ZFIN_ID;
my $expr_term;


# Do filtering of expression terms while populating the expression term hashes
# Algorithm for filtering:
# 0. Choose the data source for collecting the expression data using an if-loop
# 1. Generate stage-specific files containing all of the expression terms at each particular stage of development.
# 2. Depending on the stage variable read the corresponding file and populate a hash with stage-specific expression terms
# 3. Check each of the expression terms in the hash for sample genes on its presence in the stage-specific 
#    expression terms' hash. If it is not present in the stage-specific term hash, do not include such a term into a hash for expression terms

# In addition, the program will check whether each of the ZFIN IDs is present in the list of user-defined total genes 
# For this, it is first necessary to generate a hash with all of the user-defined ZFIN IDs. 
# Generating such a hash requires first converting IDs of the total gene list to ZFIN IDs

# initialize a hash for storing all of the tested ZFIN IDs
my %total_ZFIN_IDs = ();
my $number_total_genes = scalar @total_genes_tested;

# a condition that the number of genes is larger than 0, that is, the user used the option of entering total genes tested
if ($number_total_genes > 0){

# if loops are used to determine which conversion table should be read in or none at all
if ($format eq 'zfin_id') {
	
	# no conversion is necessary, simply assign to each ZFIN ID a value of 1 in the hash
	foreach(@total_genes_tested) {
		$total_ZFIN_IDs{$_} = 1;
	}

} elsif ($format eq 'zfin_symbol') {
	
	open ZFINSYMB, "<../zeogs_data/ZFIN_symbol-ZFIN_IDs_map.txt" or die $!;
		
		# initialize the variables necessary for holding gene IDs
		my $zfin_symb;
		my $zfin_id;
		
		# read the data from the indicated file and populate the conversion table hash 
		while(<ZFINSYMB>) {
			
			# obtain values for 
			($zfin_symb, $zfin_id) = split(/\t/, $_);
			# clean up these variables
			chomp($zfin_symb); chomp($zfin_id);
			
			if(exists $conversion_table{$zfin_symb}){
				# add a new Ensembl gene ID to the array reference stored under the key of $ensembl
				push(@{$conversion_table{$zfin_symb}}, $zfin_id);
			
			} elsif ($zfin_id =~ /ZDB-GENEP?-[0-9]{6}-[0-9]{1,4}/) # zfinID should match its correct pattern
			{
				# create a new array reference under the key of $ensembl
				$conversion_table{$zfin_symb} = [$zfin_id];
			}	
		}
		close ZFINSYMB; 
	
		# use values of gene IDs in the array to obtain the corresponding ZFIN IDs from the conversion table hash
		# and store them in the ZFIN_IDs hash
			
		foreach my $id (@total_genes_tested) {
			if (exists $conversion_table{$id}){
				# run through the array corresponding to this ID and collect unique ZFIN gene IDs
				foreach (@{$conversion_table{$id}}) {
					unless(exists $total_ZFIN_IDs{$_}) {
						$total_ZFIN_IDs{$_} = 1;
					}
				}
			}
		}
		
	} 
}


# 0. Choose the data source for collecting the expression data using an if-loop
# The first option is 'gene_summary' meaning that all the structures are taken from the gene summary pages

# initialize a hash for storing all of the stage-specific expression terms
my %expr_terms_stage = ();

if ($data_source eq 'gene_summary') {

# 1. Generate stage-specific files containing all of the expression terms at each particular stage of development.
# initialize a hash for storing all of the stage-specific expression terms
my %expr_terms_stage = ();

# read the file containing all of the stage-specific terms, which can be accessed using the stage variable value
open (STAGES, "<../zeogs_data/stages/stage$stage.txt") or die ("Cannot open this file");

# 2. Depending on the stage variable read the corresponding file and populate a hash with stage-specific expression terms
while(<STAGES>){
	chomp($_);
	$expr_terms_stage{$_}= 1;
}

close STAGES;

# open the file with all of the expression information available from ZFIN
open (READ, "<../zeogs_data/ZFIN_IDs_genes_expr_terms.txt") or die ("Cannot open this file");


# Read in expression terms and ZFIN IDs from the main file and assign expression terms to keys of the hash
# and ZFIN IDs to arrays, whose references are values of the hash
# In addition, generate a hash containing ZFIN IDs as keys of the hash and lists of expression terms as values in this hash
# For both hashes, check that each expression term exists at a given stage. 

# Go through all of the lines in the file, split them into a ZFIN ID and an expression term
while (<READ>) {
	# split each line of the text file into its tab-separated components
	($ZFIN_ID, $expr_term) = split('\t', $_);
	chomp($ZFIN_ID);
	chomp($expr_term);
	trim($expr_term);
	
	# checking that the expression term is not 'unspecified' and not 'whole organism'
	# The term 'whole organism' is problematic because for many genes it is not well-defined
	# temporarily and many genes are expressed in the whole organism only very early in development 
	# checking the existence of a hash entry with a particular ZFIN ID
	
	# generate a hash where each ZFIN ID has a list of all of its expression terms
	if(($expr_term ne 'unspecified') && ($expr_terms_stage{$expr_term}) && ($expr_term !~ /whole organism/)) {
		if( exists $all_genes_patterns{$ZFIN_ID}) {
			push(@{$all_genes_patterns{$ZFIN_ID}}, $expr_term); 
		} else {
			$all_genes_patterns{$ZFIN_ID} = [$expr_term];
		}
	}
	
	# generate a hash where each expression term has a list of all of ZFIN IDs for genes expressed in this expression term
	if(($expr_term ne 'unspecified') && ($expr_terms_stage{$expr_term}) && ($expr_term !~ /whole organism/)) {
		if( exists $patterns_all_genes{$expr_term}) {
			push(@{$patterns_all_genes{$expr_term}}, $ZFIN_ID); 
		} else {
			$patterns_all_genes{$expr_term} = [$ZFIN_ID];
		}
	}
  }
 close(READ);

} elsif($data_source eq 'wt_expression'){
	
# open the file with all of the expression information available from ZFIN
open (READ, "<../zeogs_data/wildtype_expression_staged.txt") or die ("Cannot open this file");

# Go through all of the lines in the file, split them into a ZFIN ID and an expression term
while (<READ>) {
	# split each line of the text file into its tab-separated components
	my ($ZFIN_ID, $stage_info, $expr_term) = split('\t', $_);
	chomp($expr_term);
	
	# checking that the expression term is not 'unspecified' and not 'whole organism'
	# The term 'whole organism' is problematic because for many genes it is not well-defined
	# temporarily and many genes are expressed in the whole organism only very early in development 
	# checking the existence of a hash entry with a particular ZFIN ID
	
	# for each line test whether the stage information in the file is equal to the selected stage option
	if(($stage_info == $stage) &&($expr_term ne 'unspecified') && ($expr_term !~ /whole organism/)) {
	
	# generate a hash where each ZFIN ID has a list of all of its expression terms
		if( exists $all_genes_patterns{$ZFIN_ID}) {
			push(@{$all_genes_patterns{$ZFIN_ID}}, $expr_term); 
		} else {
			$all_genes_patterns{$ZFIN_ID} = [$expr_term];
		}
	
	# generate a hash where each expression term has a list of all of ZFIN IDs for genes expressed in this expression term
		if( exists $patterns_all_genes{$expr_term}) {
			push(@{$patterns_all_genes{$expr_term}}, $ZFIN_ID); 
		} else {
			$patterns_all_genes{$expr_term} = [$ZFIN_ID];
		}
	}
	
		
}
close(READ);
}

# The next step is to adjust the total available information about gene expression patterns of zebrafish genes
# to those genes that were actually tested


if ($number_total_genes > 0) { # a list of total genes was input to the program

	# initialize a hash for storing expression terms for ZFIN IDs from the previous step of the script
	my %total_genes_patterns = ();
	my %patterns_total_genes = ();

	# go through all of the ZFIN IDs and make a slice of the original expression term hash
	foreach (keys %total_ZFIN_IDs){
		if (exists $all_genes_patterns{$_}){
			$total_genes_patterns{$_} = $all_genes_patterns{$_};
		} 
	}

	# generate a hash with expression terms as keys and lists of ZFIN IDs as values from the user input
	# for the total tested genes

	foreach my $zid (keys %total_genes_patterns) {
		foreach (@{$total_genes_patterns{$zid}}){
			if( exists $patterns_total_genes{$_}) {
				push(@{$patterns_total_genes{$_}}, $zid); 
			} else {
				$patterns_total_genes{$_} = [$zid];
			}
		}
	} 
	
	# update the hashes containing information on expression of all known genes with the 
	# information obtained from the list of total tested genes
	
	%all_genes_patterns = %total_genes_patterns;
	%patterns_all_genes = %patterns_total_genes;
	
}




# 2. Get a slice of this hash for the available ZFIN IDs.
# 3. Put ZFIN IDs without expression term annotation into a special array.

# initialize a hash for storing expression terms for ZFIN IDs from the previous step of the script
my %sample_genes_patterns = ();
my %patterns_sample_genes = ();
# initialize a hash for storing ZFIN IDs having no expression terms assigned
my @sample_ZFINIDs_no_expr_terms = ();

# go through all of the ZFIN IDs and make a slice of the original expression term hash
# as well as put ZFIN without expression terms into a special array
foreach (@ZFIN_IDs){
	if (exists $all_genes_patterns{$_}){
		$sample_genes_patterns{$_} = $all_genes_patterns{$_};
	} else {
		push(@sample_ZFINIDs_no_expr_terms, $_);
	}
}

# generate a hash with expression terms as keys and lists of ZFIN IDs as values from the user input

foreach my $zid (keys %sample_genes_patterns) {
	foreach (@{$sample_genes_patterns{$zid}}){
		if( exists $patterns_sample_genes{$_}) {
			push(@{$patterns_sample_genes{$_}}, $zid); 
		} else {
			$patterns_sample_genes{$_} = [$zid];
		}
	}
} 

# Obtain the total number of genes with expression terms both for the sample and "all genes" hashes
# These total counts can then be used for statistical testing

my $total_sample = scalar keys %sample_genes_patterns;
my $total_all_genes    = scalar keys %all_genes_patterns;


#####################################################
# Step 5 - Generating a data structure for the output
#####################################################

# initialize a new hash for storing all the necessary information for output on the webpage
# The hash will contain the following fields in its array reference: 
# 1) a link to ZFIN Anatomy ontology page for the expression term 
# 2) count of genes having this expression term in the sample
# 3) a reference to the gene list  
# 4) a gene count for all of the genes with expression terms
# 5) P-value for the hypergeometric test

# and two optional fields if the user has decided that a multiple testing correction should be done:

# 6) a corrected P-value
# 7) a 'YES' or 'NO' value depending on whether the original P-value has exceeded the Benjamini-Hochberg procedure's threshold

my %expression_terms_table = ();

# obtain an URL for each of the ZFIN anatomical terms

# initialize a hash for storing URLs of 
my %expression_terms_URLs = ();

# open a file containing identifier information for all of the zebrafish
# expression terms, read in the identifiers and names and make a hash, where a name of an expression term
# is a key and its URL is its value.
# URL has the following template: http://zfin.org/action/anatomy/anatomy-view/

open (TERMS, "<../zeogs_data/anatomy_item.txt") or die ("Cannot open this file");

while(<TERMS>){
	# read all of the elements from each line of the file
	my @items = split(/\t/, $_);

	# store expression term names and their corresponding addresses in a hash
	$expression_terms_URLs{$items[1]} = "http://zfin.org/action/anatomy/anatomy-view/$items[0]";
}
close TERMS;

# Now the final data structure should be populated based on the information collected so far
# and according to the structure outlined above 

foreach my $expression_term (keys %patterns_sample_genes) {
	$expression_terms_table{$expression_term} = [$expression_terms_URLs{$expression_term}, 
	scalar @{$patterns_sample_genes{$expression_term}}, \@{$patterns_sample_genes{$expression_term}},
	scalar @{$patterns_all_genes{$expression_term}}];
}


##################################################################
# Step 6 - Performing a hypergeometric statistical 
#          test on counts in sample data versus the all genes data 
##################################################################


# the simple hypergeometric test is the default behaviour of the script
# meaning that there are no if loops to specify what exactly should be done

# Here, $n is the total number of genes with an expression terms among all genes,
# $m is the difference between the total number of genes and $n
# N is the total number of genes in the sample
# $i is the number of genes with an expression term among sample genes
# hypercdf_pvalue($n, $m, $N, $i)

my $total_sample = scalar keys %sample_genes_patterns;
my $total_all_genes    = scalar keys %all_genes_patterns;
my $total_anatomy_terms    = scalar keys %patterns_all_genes;

# In the P-value calculation below, the the following values are used
# $n = $expression_terms_table{$expression_term}[3]
# $m = $total_all_genes -  $expression_terms_table{$expression_term}[3]
# $N = $total_sample
# $i = $expression_terms_table{$expression_term}[1]

foreach my $expression_term (keys %expression_terms_table) {
	$expression_terms_table{$expression_term}[4] = hypercdf_pvalue($expression_terms_table{$expression_term}[3], $total_all_genes -  $expression_terms_table{$expression_term}[3], $total_sample, $expression_terms_table{$expression_term}[1]);
}



#################################
# Benjamini-Hochberg correction
#################################



# The algorithm will be commented stepwise
# The code for calculating the adjusted P-value threshold for the current P-values is external, 
# the webpage for it being given in the subroutine itself  

# 1. Sort the %expression_terms_table hash using P-value fields of its array values

# initialize an array for ordered P-values
my @Pvalues = ();
  
foreach ( sort {$expression_terms_table{$a}[4] <=> $expression_terms_table{$b}[4]} keys %expression_terms_table) {
# 2. Obtain sorted P-values in a separate array
	push(@Pvalues, $expression_terms_table{$_}[4]);
}


# 3. Sort the hash entries again using the P-value criterium and do the following:
#    * calculate a corrected P-value = total/i*Pi and store it in the 6th field of the array reference for each hash entry
#    * compare the initial P-value with the calculated BH P-value threshold and write a value of 'YES' or 'NO' depending on
#    whether this value is smaller than the threshold or not. The value of 'YES' means that the anatomical term is significantly
#    over-represented after the multiple testing correction

# get the total number of expression terms for the purpose of multiple testing correction
my $number_expr_terms = scalar keys %expression_terms_table;

# initialize a rank counter for calculating the corrected P-value 
my $rank = $number_expr_terms;


foreach my $term ( sort {$expression_terms_table{$b}[4] <=> $expression_terms_table{$a}[4]} keys %expression_terms_table) {
		# calculate the corrected P-value
		$expression_terms_table{$term}[5] = ($number_expr_terms/$rank)*$expression_terms_table{$term}[4]; 
		# The P-values are probabilities. Therefore, they should be equal 1 at most
		if ($expression_terms_table{$term}[5] >= 1) {$expression_terms_table{$term}[5] = 1;}
		
		# increment the rank value 
		$rank--;
		
		# Do P-value comparison with the P-value threshold and write the result to the array
		if ($expression_terms_table{$term}[5] >= 0.1) {
			$expression_terms_table{$term}[6] = 'NO';
		} else { $expression_terms_table{$term}[6] = 'YES'; }
}

###########################################################################
# Step 7 - Sorting the data based on the options from the user
###########################################################################

# initialize an array for storing anatomical terms in the sorted order
my @expression_terms_sorted = ();

# perform the sorting and store the results in the array initialized above
if ($sorting_option eq 'count'){ 
	
	# sorting based on the frequency of expression terms among sample genes 
	@expression_terms_sorted = sort {$expression_terms_table{$b}[1] <=> $expression_terms_table{$a}[1]} keys %expression_terms_table;

} elsif($sorting_option eq 'Pvalue') {
	# sorting based on the original P-value from the hypergeometric test 
	@expression_terms_sorted = sort {$expression_terms_table{$a}[4] <=> $expression_terms_table{$b}[4]} keys %expression_terms_table;

} elsif($sorting_option eq 'Fraction_total') {
	# sorting based on the frequency of expression terms among sample genes 
	@expression_terms_sorted = sort {$expression_terms_table{$b}[1]/$expression_terms_table{$b}[3] <=> $expression_terms_table{$a}[1]/$expression_terms_table{$a}[3]} keys %expression_terms_table;

}


#######################################################################################
# Step 8 - Generating an output page and writing files to a disk with links on the page
#######################################################################################

# generating a random folder name and storing the newly generated data in that folder
my $rand_folder_name = &generate_random_string(10);

# create a directory with the $rand_folder name
&makeDir("../tmp/results/$rand_folder_name");

# writing the raw data to files for later use for download links
open (OUTPUT, ">../tmp/results/$rand_folder_name/table_ranked_anatomy_terms.txt") or die ("Cannot open this file");

printf OUTPUT "%50s", 'anatomical term';
print  OUTPUT "\t";
printf OUTPUT "%18s", 'Sample gene count';
print  OUTPUT "\t";
printf OUTPUT "%16s", 'Total gene count';
print  OUTPUT "\t"; 
printf OUTPUT "%7s", 'P-value';
print  OUTPUT "\t";
printf OUTPUT "%18s", 'corrected P-value'; 
print  OUTPUT "\t";
printf OUTPUT "%15s", 'Below threshold?';
print  OUTPUT "\n";

foreach (@expression_terms_sorted) {
	printf OUTPUT "%50s", $_; 
	print  OUTPUT "\t";
	printf OUTPUT "%18u", $expression_terms_table{$_}[1]; 
	print  OUTPUT "\t";
	printf OUTPUT "%16u", $expression_terms_table{$_}[3];
	print  OUTPUT "\t";
	
	# output P-values in one of two different formats depending on how small they are
	if($expression_terms_table{$_}[4] < 0.00001) {
		printf OUTPUT "%7.5e", $expression_terms_table{$_}[4];
	} else {
		printf OUTPUT "%7.5f", $expression_terms_table{$_}[4];
	}
	print OUTPUT "\t";
	printf OUTPUT "%18.5f", $expression_terms_table{$_}[5];
	print  OUTPUT "\t";
	printf OUTPUT "%15s", $expression_terms_table{$_}[6];
	print  OUTPUT "\n";
}
close OUTPUT;

# before writing the genes for each anatomical term to the file 
# we should obtain gene symbols and short descriptions for using genes' ZFIN IDs

# first, one needs to initialize a hash, where ZFIN IDs are keys and symbols as well as descriptions are stored in an array referenced by the value of the hash 

my %ZFINID_symb_desc = ();

open (GENEINFO, "<../zeogs_data/ZFIN_IDs-ZFIN_symbols.txt") or die ("Cannot open this file");

while(<GENEINFO>){
	# read all of the elements from each line of the file
	my @gene_infos = split(/\t/, $_);
		# remove any extraneous material from each text element
		chomp($gene_infos[0]);
		$gene_infos[0] = trim($gene_infos[0]);
		chomp($gene_infos[1]);
		$gene_infos[1] = trim($gene_infos[1]);
		chomp($gene_infos[2]);
		$gene_infos[2] = trim($gene_infos[2]);
		
		# populate the hash of interest
		# Here, the first element is ZFIN ID, the second is gene symbol and the third is the gene description
		$ZFINID_symb_desc{$gene_infos[0]} = [$gene_infos[1], $gene_infos[2]];
}
close GENEINFO;



# write lists of genes for each anatomical term encountered in the sample
open (OUTGENES, ">../tmp/results/$rand_folder_name/genes_anatomy_terms.txt") or die ("Cannot open this file");

foreach (@expression_terms_sorted) {
	print OUTGENES "\n$_\t\t\t\n\n";
	foreach my $g(@{$expression_terms_table{$_}[2]}) {print OUTGENES "\t$g\t$ZFINID_symb_desc{$g}[0]\t$ZFINID_symb_desc{$g}[1]\n";}
	
}
close OUTGENES;

# In addition, a list of genes without expression terms needs to be output
open (NOTERMS, ">../tmp/results/$rand_folder_name/genes_no_terms.txt") or die ("Cannot open this file");

print NOTERMS "The following genes currently do not have any associated anatomical terms:\n\n";

foreach my $zid (@sample_ZFINIDs_no_expr_terms) {
	print NOTERMS "$zid\t$ZFINID_symb_desc{$zid}[0]\t$ZFINID_symb_desc{$zid}[1]\n";
}

close NOTERMS;



# printing the all the results to an html page
# <!--The section containing links to the generated files available for download-->
print "<h2>The output of ZEOGS for your input data</h2>";

print "<h3>All of the information resulting from your input is available in the following files:</h3>";


print "<h4>&nbsp;&nbsp;&nbsp;&nbsp;<a href=\"http://zeogs.x10.bz/tmp/results/$rand_folder_name/table_ranked_anatomy_terms.txt\">File 1: Table of ranked anatomical terms</a></h4>";

print "<h4>&nbsp;&nbsp;&nbsp;&nbsp;<a href=\"http://zeogs.x10.bz/tmp/results/$rand_folder_name/genes_anatomy_terms.txt\">File 2: Tables of gene identifiers for each anatomical term</a></h4>";

print "<h4>&nbsp;&nbsp;&nbsp;&nbsp;<a href=\"http://zeogs.x10.bz/tmp/results/$rand_folder_name/genes_no_terms.txt\">File 3: A table of gene identifiers without known anatomical terms</a></h4>";


# printing expression terms from the %patterns_samples_genes hash in the form of a table
print "<table class=\"result\">";
print "<thead class=\"result\"><tr><th> Anatomical term </th><th> Sample gene count </th><th> Total gene count </th><th> P-value </th><th> corrected P-value </th><th> Below threshold? </th></tr></thead>";

foreach (@expression_terms_sorted) {
	if ($expression_terms_table{$_}[5] <= 0.1) {
	
		print "<tr id=\"high\"><td><a href=\"$expression_terms_table{$_}[0]\">$_</a></td><td>$expression_terms_table{$_}[1]</td><td>$expression_terms_table{$_}[3]</td><td>";
	} elsif ( $expression_terms_table{$_}[4] <= 0.1) {
		print "<tr id=\"signif\"><td><a href=\"$expression_terms_table{$_}[0]\">$_</a></td><td>$expression_terms_table{$_}[1]</td><td>$expression_terms_table{$_}[3]</td><td>";
	
	} else {
		print "<tr id=\"nonsignif\"><td><a href=\"$expression_terms_table{$_}[0]\">$_</a></td><td>$expression_terms_table{$_}[1]</td><td>$expression_terms_table{$_}[3]</td><td>";
	
	}

	
	# output P-values in one of two different formats depending on how small they are
	if($expression_terms_table{$_}[4] < 0.00001) {
		printf "%12.5e", $expression_terms_table{$_}[4];
	} else {
		printf "%12.5f", $expression_terms_table{$_}[4];
	}
	
	print "</td><td>";
	printf "%12.5f", $expression_terms_table{$_}[5]; 
	print "</td><td>$expression_terms_table{$_}[6]</td></tr>";
}

# FINISH THE TABLE
print "</table>";

print "<p>The number of genes with anatomical terms in the sample is $total_sample</p>";
print "<p>The number of genes with anatomical terms for all genes is $total_all_genes</p>";
print "<p>The number of anatomical terms encountered is $number_expr_terms</p>";
print "<p>The number of anatomical terms encountered for ALL GENES is $total_anatomy_terms</p>";


print <<"HTML_footer";

<br />
<br />
<br />
		
<p id="footer"> Copyright notice: ZEOGS website has been developed by Sergey Prykhozhij at Max-Planck Institute for molecular genetics, Berlin, Germany, and at IWK Health Centre and Dalhousie University in 2012. The last update of the data sets used was in April 2015.</p>

HTML_footer

delete_folders();

##################################
# Subroutines
##################################

sub hypercdf_pvalue {
	# obtain necessary data for the calculation
	# Here, $n is the total number of genes with expression term among all genes,
	# $m is the difference between the total number of genes and $n
	# N is the total number of genes in the sample
	# $i is the number of genes with an expression term among sample genes
	
	my ($n, $m, $N, $i) = @_;
	my $hypercdf = 0;
	
	# Do the calculation in a loop
	for (my $iref=$i; $iref <= min($N,$n); $iref++) {
		$hypercdf += hypergeom($n,$m,$N,$iref);
	}
	
	# return the result of the subroutine 
	return $hypercdf;
	
}

sub logfact {
	   return gammln(shift(@_) + 1);
	}

sub hypergeom {
	   my ($n, $m, $N, $i) = @_;
	   my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
	   my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact($N-$i)+logfact($m+$n);
	   return exp($loghyp1 - $loghyp2);
}

sub gammln {
	  my $xx = shift;
	  my @cof = (76.180091729471457, -86.505320329416776,
				 24.014098240830911, -1.231739572450155,
				 0.12086509738661, -5.395239384953e-06);
	  my $y = my $x = $xx;
	  my $tmp = $x + 5.5;
	  $tmp -= ($x + 0.5) * log($tmp);
	  my $ser = 1.0000000001900149;
	  foreach my $j (0 .. 5) {
		 $ser += $cof[$j] / ++$y;
	  }
	  -$tmp + log(2.5066282746310007*$ser/$x);
}



sub Error {
	print "Couldn't open temporary file: $!";
	exit;
}

sub trim {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub delete_folders {

	my $OneDay = 60*60*24; #60secs x 60mins * 24hrs
	$OneDay = time - $OneDay;
	my $path="../tmp/results"; 

	opendir(FP,"$path"); 
	my @del=readdir(FP); 
	closedir(FP); 
	
	my $count_dir;
	$count_dir = scalar @del;
	
	if ($count_dir > 0) {
	
	foreach my $i (@del) { 
		my $fp = "$path/$i"; 
		
		if (-M $fp > 1) {
			
			# first obtain a list of the files inside a directory
			opendir (DIR, "$fp");
			my @FILES = readdir(DIR);
			closedir (DIR);
			## DELETE THE .TXT FILES THAT ARE OLDER THAN 1 DAY

			foreach my $FILE (@FILES) {
				  unlink("$fp/$FILE");
			}
		}
		
		# now that the directory is empty one can apply the rmdir function to remove it
		# the condition for that is that there were any folders present inside the 'root' directory
		rmdir($fp);
		
		}
	} 
}



sub upload_file{
	# the information necessary for the upload operation to happen
	my ($filename,$filehandle,$path) = @_;
	# initialize an array to store the identifiers from the file 
	my @geneids = ();
	# where the files will be stored 
	my $target = $path.'/'.$filename;
	
	# checking whether the file already exists and if it doesn't the upload will occur
	if(-e $target){
		print "The target file already exists.\n";
		exit(0);
	}
	else{
		binmode $filehandle;
		open(TARGET,">$target") or die $!;
		binmode TARGET;
		my ($buffer);
	while(read $filehandle,$buffer,1024){
		print TARGET $buffer;
	}
	# store all the gene IDs from the file into the array
	
	close TARGET;
	
	# open file
	open(READ, "<$target") or die("Unable to open file");
 
	# read file into an array
	@geneids = <READ>;
 
	# close file 
	close(READ);
	
	# delete file that has been uploaded
	unlink($target);
	
	# return an array containing the gene IDs input by the user
	return @geneids;
 }
}



# This function generates random strings of a given length
###########################################################
# Written by Guy Malachi http://guymal.com
# 18 August, 2002
###########################################################
sub generate_random_string
{
	my $length_of_randomstring=shift;# the length of the random string to generate

	my @chars=('a'..'z','A'..'Z','0'..'9','_');
	my $random_string;
	foreach (1..$length_of_randomstring) 
	{
		# rand @chars will generate a random 
		# number between 0 and scalar @chars
		$random_string.=$chars[rand @chars];
	}
	return $random_string;
}



# This subroutine has been written by Robert Ketter
# and posted on a discussion forum on Webdeveloper.com
# it was adapted for the purposes of my script

sub makeDir {
  my $perm = 755;
  my $new_dir = shift;
  
  if (-e "$new_dir"){ 
	problem("Directory (<b> $new_dir </b>) already exists.\n") 
  } # Checks for existing directories with the same name

  mkdir ($new_dir,$perm) || problem("Error making Directory (<b> $new_dir </b>)\n");


}# end sub makedir

sub problem{
  # I always make my own sub to report problems. I find it alot easier to understand
  # an error if I explain what happened in my own words.
  print "Content-type: text/html\n\n";
  print "$_[0]\n";
  exit;
  }# end sub problem